library(rattaca)
library(argparse)
library(viridis)
library(scales)

################################################
# ---------------  ARGUMENTS  --------------- ##
################################################


parser <- argument_parser(
    argument_def('rattaca', type = 'character', default_val = 'rattaca',
        help = 'arbitrary positional argument included to enable argument parsing'),
    argument_def('--pred_type', type = 'character', required =  TRUE, 
        help = 'the type of predictions to run, either "rattaca" for colony 
                predictions or "ad_hoc" for predictions on requested RFIDs'),
    argument_def('--sample_type', type = 'character', required =  TRUE, 
        help = 'the SNP sampling strategy: "ldprune" or "random"bpar'),
    argument_def('--trait', type = 'character', required = TRUE, 
        help = 'trait name, as used in trait list and column headers'),
    argument_def('--trait_dict', type = 'character', required = TRUE, 
        help = 'path to trait dictionary csv'),
    argument_def('--rfid_list', type = 'character', required = FALSE,
        help = 'path to file listing all test RFIDs (for ad-hoc predictions)'),
    argument_def('--generation', type = 'character', required = TRUE,
        help = 'generation name/number'),
    argument_def('--train_geno_input', type = 'character', required = TRUE,
        help = 'path to genotype training data (w/o plink file extension)'),
    argument_def('--test_geno_input', type = 'character', required = TRUE,
        help = 'path to genotype test data (w/o plink file extension)'),
    argument_def('--train_pheno_file', type = 'character', required = TRUE,
        help = 'path to phenotype training data (csv)'),
    argument_def('--proj_dir', type = 'character', required = TRUE,
        help = 'top directory path for the rattaca generation or ad hoc request'),
    argument_def('--bpar_dir', type = 'character', required = TRUE,
        help = 'directory in which to save official bpar files for future predictions'),
    argument_def('--plink1', type = 'character', required = TRUE,
        help = 'path to plink v1 software'),
    argument_def('--plink2', type = 'character', required = TRUE,
        help = 'path to plink v2 software'), 
    argument_def('--maf_cutoff', type = 'double', default_val = 0.01,
        help = 'minor allele frequency cutoff: variants with MAF below this 
                threshold will be dropped from genotype data'),
    argument_def('--missing_cutoff', type = 'double', default_val = 0.1,
        help = 'genotype missingness cutoff: variants with missing data in rats 
                at a higher frequency than this cutoff will be dropped from 
                genotype data'),
    argument_def('--hwe_cutoff', type = 'double', default_val = 1e-6,
        help = 'Hardy-Weinberg p-value cutoff: variants whose p-value from a
                Hardy-Weinberg exact test are below this cutoff will be dropped 
                from genotype data'),
    argument_def('--ld_r2', type = 'double', default_val = 0.99,
        help = 'linkage disequilibrium pairwise r^2 threshold: if LD-pruning, 
                one variant per pair will be kept in the genotype dataset if the 
                pairwise correlation excedes this threshold'),
    argument_def('--ld_window_size', type = 'double', default_val = 1000,
        help = 'linkage disequilibrium window size: if LD-pruning, pairwise 
                correlations will be estimated between all variants in a genomic 
                window containing this number of SNPs'),
    argument_def('--ld_step_size', type = 'double', default_val = 100, 
        help = 'linkage disequlibrium step size: if LD-pruning, pairwise 
                correlations will be recalculated after sliding the genomic 
                window by this number of SNPs'),
    argument_def('--n_samples', type = 'double', default_val = 5, 
        help = 'number of random SNP samples to execute'),
    argument_def('--n_snps', type = 'double', default_val = 50000, 
        help = 'number of random SNPs to retain per random sample')
)

# parse the command-line arguments
args <- parser(commandArgs())
list2env(args, envir = .GlobalEnv)

# set up directories
results_dir <- file.path(proj_dir, 'results')
dir.create(results_dir, showWarnings = FALSE)
trait_dir <- file.path(results_dir, trait)
dir.create(trait_dir, showWarnings = FALSE)
train_data_dir <- file.path(trait_dir, 'train')
dir.create(train_data_dir, showWarnings = FALSE)
test_data_dir <- file.path(trait_dir, 'test')
dir.create(test_data_dir, showWarnings = FALSE)
pca_dir <- file.path(trait_dir, 'pca')
dir.create(pca_dir, showWarnings = FALSE)
snp_dir <- file.path(trait_dir, 'snps')
dir.create(snp_dir, showWarnings = FALSE)

# output directories for model parameter (.bpar) files
if (exists('trait_dict')) {
    dict <- read.csv(trait_dict)
    proj_name <- dict[dict$trait==trait, 'project_name']
    variable_used <- dict[dict$trait==trait, 'variable_used']
    trait_source <- dict[dict$trait==trait, 'source_file']
} else {
    proj_name <- format(Sys.time(), '%Y%m%d-%H:%M:%S')
    variable_used <- NULL
    trait_source <- NULL
}
if (is.null(variable_used)) {
    trait_bpar_dir <- file.path(bpar_dir, paste0(proj_name, '-', trait))
    bpar_file <- file.path(trait_bpar_dir, paste0(trait,'.bpar'))
} else {
    trait_bpar_dir <- file.path(bpar_dir, paste0(proj_name, '-', variable_used))
    bpar_file <- file.path(trait_bpar_dir, paste0(variable_used,'.bpar'))
}
dir.create(trait_bpar_dir, showWarnings = FALSE)

cat('\n')
printout(paste('Model fits for', trait))
cat('\n')
cat('train_pheno_file:', train_pheno_file, '\n')
cat('train_geno_input:', train_geno_input, '\n')
cat('proj_dir:', proj_dir, '\n')
cat('plink1:', plink1, '\n')
cat('plink2:', plink2, '\n')
cat('generation:', generation, '\n')


#################################################
## ---------------  DATA PREP  --------------- ##
#################################################


# # get all RFIDs to use in the test set
# test_rfids <- read.csv(rfid_list)$rfid
# test_rfids_file <- format_ids_file(test_rfids, test_data_dir, paste0(trait, '_test_rfids'))

# write phenotyped rfids to a file to enable subsetting using Plink
phtyped_rfids_file <- get_phenotyped_ids(
    phenotype_data = train_pheno_file, 
    id_column = 'rfid', 
    trait_column = trait,
    output_dir = train_data_dir)$ids_file

# make the phenotype training dataset
pheno_train <- prep_trait_data(train_pheno_file, 'rfid', trait, 'numeric')


get_common_snpset2 <- function(train_genotypes, test_genotypes, out_prefix=NULL) {
    
    if (length(train_genotypes) == 1 && is.character(train_genotypes) && file.exists(paste0(train_genotypes, '.bim'))) {
        # read in snp if provided from a plink bim file
        all_train_snps <- genio::read_bim(paste0(train_genotypes, '.bim'))
        all_train_snps <- all_train_snps$id
    } else if (length(train_genotypes) == 1 && is.character(train_genotypes) && 
              grepl('.bpar$', train_genotypes) && file.exists(train_genotypes)) {
        # read in snps if provided from a bpar file
        all_train_snps <- get_snps_from_bpar(train_genotypes)
    } else if (!is.list(train_genotypes) && length(train_genotypes) > 1) {
        # keep snps as-is if provided as a vector of variant IDs
        all_train_snps <- train_genotypes
    } else if (is.list(train_genotypes) && 'geno_file' %in% names(train_genotypes)) {
        # extract snps from genotype data if provided as a rattaca dataset
        all_train_snps <- colnames(train_genotypes$geno)
    } else if (is.matrix(train_genotypes)) {
        # extract snps from column names if provided as a genotype matrix
        all_train_snps <- colnames(train_genotypes)
    } else {
        cat('Check train_genotypes format: \n')
        print(str(train_genotypes))
        stop("Unrecognized format for train_genotypes")
    }   

    if (length(test_genotypes) == 1 && is.character(test_genotypes) && file.exists(paste0(test_genotypes, '.bim'))) {
        # read in snp if provided from a plink bim file
        all_test_snps <- genio::read_bim(paste0(test_genotypes, '.bim'))
        all_test_snps <- all_test_snps$id
    } else if (length(test_genotypes) == 1 && is.character(test_genotypes) && 
              grepl('.bpar$', test_genotypes) && file.exists(test_genotypes)) {
        # read in snps if provided from a bpar file
        all_test_snps <- get_snps_from_bpar(test_genotypes)
    } else if (!is.list(test_genotypes) && length(test_genotypes) > 1) {
        # keep snps as-is if provided as a vector of variant IDs
        all_test_snps <- test_genotypes
    } else if (is.list(test_genotypes) && 'geno_file' %in% names(test_genotypes)) {
        # extract snps from genotype data if provided as a rattaca dataset
        all_test_snps <- colnames(test_genotypes$geno)
    } else if (is.matrix(test_genotypes)) {
        # extract snps from column names if provided as a genotype matrix
        all_test_snps <- colnames(test_genotypes)
    } else {
        cat('Check test_genotypes format: \n')
        print(str(test_genotypes))
        stop("Unrecognized format for test_genotypes")
    }   

    # save all snps common to train/test sets from which to sample a common set of snps
    all_snps <- intersect(all_test_snps, all_train_snps)
    
    out <- list(snps = all_snps, file = NULL)

    # save a plink-formatted variant list
    if (!is.null(out_prefix)) {
        train <- basename(train_genotypes)
        test <- basename(test_genotypes)
        outfile <- paste0(out_prefix, '_common_snps')
        writeLines(all_snps, outfile)
        out$file <- outfile
    }

    return(out)
}

# get the common snp set between training and test genotype data
all_snps <- get_common_snpset2(train_geno_input, test_geno_input, out_prefix=file.path(snp_dir, trait))


# make the base training dataset
printout('Producing base training genotype dataset')

base_train_geno <- make_plink_dataset(
    input_genotypes = train_geno_input, 
    output_dir = train_data_dir,
    outfile_prefix = paste0('base_', trait),
    samples_to_keep = phtyped_rfids_file,
    snps_to_keep = all_snps$file,
    maf_cutoff = maf_cutoff,         
    missing_cutoff = missing_cutoff,
    hwe_cutoff = hwe_cutoff,        
    return_data = F)
base_train_file <- base_train_geno$geno_file

# produce the actual training genotype dataset
if (sample_type == 'ldprune') {
    
    printout('LD-pruning training genotypes')
    geno_train <- make_plink_dataset(
        input_genotypes = base_train_file, 
        output_dir = train_data_dir,
        outfile_prefix = trait,
        ld_prune = T,
        ld_r2 = ld_r2,
        ld_window_size = ld_window_size,
        ld_step_size = ld_step_size,
        snp_directory = snp_dir,
        return_data = T
        )
    train_geno_file <- geno_train$geno_file
    # geno_train <- list(geno_train)

} else if (sample_type == 'random') {

    # produce smaller snp samples & save files listing snps
    printout(paste('Producing', n_samples, 'random SNP samples'))
    use_snps <- sample_snps(
        input_snps = all_snps$snps, 
        keep = n_snps, 
        save_file = T, 
        output_dir = snp_dir, 
        iterations = n_samples
        )

    printout('Randomly sampling SNPs from training genotypes')
    geno_train <- sample_snps_from_plink_files(
        input_genotypes = train_geno_file,
        snp_directory = snp_dir,
        output_dir = train_data_dir,
        n_samples = 5
        )
}


#################################################
## --------------  MODEL FITS  --------------- ##
#################################################


# model fits for LD-pruned training data
if (sample_type == 'ldprune') {

    train_dat <- align_data(geno_train, pheno_train, trait)

    # perform k-fold cross-validation
    printout('Performing k-fold cross-validation')
    cv_out <- kfold_cv(train_dat, 5, trait_dir)
    save_cv_results(cv_out, trait_dir)

    # plot cross-validation results for the best performing model
    printout('Plotting cross-validation results')
    plot_kfold(cv_out, trait_dir)

    # fit a model to the entire training dataset
    train_fit <- fit(train_dat$pheno, train_dat$geno)

# model fits for randomly-sampled training data
} else if (sample_type == 'random') {

    train_dat <- list()
    cv_list <- list()

    # loop over training snp samples
    for (i in 1:length(geno_train)) {

        # align the training data for prediction
        train_dat[[i]] <- align_data(geno_train[[i]], pheno_train, trait)

        # perform k-fold cross-validation
        printout('Performing k-fold cross-validation')
        printout(paste('CV on snp set #', i))
        cv_list[[i]] <- kfold_cv(train_dat[[i]], 5)
    }

    # save the best performing model
    best_fold <- best_kfold_mod(cv_list)
    cv_out <- cv_list[[best_fold]]
    train_dat <- train_dat[[best_fold]]

    # fit a model to the entire training dataset using the best-performing SNP sample
    train_fit <- fit(train_dat$pheno, train_dat$geno)

}

write_pars2 <- function(
    trait, trait_file, genotype_prefix, pars, filename,
    genotype_source = NULL, # the training geno input file (eg. round 10.5)
    trait_source = NULL, # the pheno datafile used to produce the trait_file
    trait_variable = NULL # the official trait variable name from trait_source
) {
    META_PREFIX <- "##"
    META_KEY_VAL_DELIM <- "="
    META_VALID_KEY_REGEX <- ".+"
    META_VALID_VAL_REGEX <- ".+"
    RECORD_DELIM <- "\n"
    FIELD_DELIM <- ","
    ENCODING <- "UTF-8"
    HEADER_PREFIX <- "#"

    sess_info <- utils::sessionInfo()

    # append .bpar file extension to output file name if needed
    if (!grepl('.bpar$', filename)) {
        filename <- paste0(filename, '.bpar')
    }

    fconn <- file(description=filename,
                  open="wt",
                  encoding=ENCODING)

    writeLines(c(paste(paste0(META_PREFIX,"date"),
                       Sys.time(),
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "user"),
                       Sys.getenv("USER"),
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "trait_name"),
                       trait,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX,"trait_file"),
                       trait_file,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "trait_source"),
                       trait_source,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "trait_source_variable"),
                       trait_variable,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "genotypes_file"),
                       genotype_prefix,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "genotypes_source"),
                       genotype_source,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX,"R_version"),
                       paste(sess_info$R.version$major,
                             sess_info$R.version$minor,
                             sep="."),
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "platform"),
                       sess_info$R.version$platform,
                       sep=META_KEY_VAL_DELIM)),
               con=fconn)

    pkg_info <- vector(mode="character",
                       length(utils::sessionInfo()$otherPkgs))
    i <- 1
    for (tmp_pkg_info in sess_info$otherPkgs)
    {
        tmp <- paste(tmp_pkg_info$Package,
                     tmp_pkg_info$Version,
                     sep=":")

        pkg_info[i] <- paste(paste0(META_PREFIX, "package"),
                             tmp,
                             sep=META_KEY_VAL_DELIM)
        i <- i + 1
    }

    writeLines(pkg_info, con=fconn)
    
    writeLines(c(paste(paste0(META_PREFIX, "n_snps"),
                       length(pars$u),
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "variance_u"),
                       pars$Vu,
                       sep=META_KEY_VAL_DELIM),
                paste(paste0(META_PREFIX, "variance_e"),
                       pars$Ve,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "log_likelihood"),
                       pars$LL,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "intercept"),
                       pars$beta,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "intercept_se"),
                       pars$beta.SE,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "r_sq"),
                       pars$r_sq,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "pearson_corr"),
                       pars$pearson_corr,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "spearman_corr"),
                       pars$spearman_corr,
                       sep=META_KEY_VAL_DELIM)),
               con=fconn)

    writeLines(paste(paste0(HEADER_PREFIX, "var_id"),
                     "random_effect",
                     "random_effect_se",
                     sep=FIELD_DELIM),
               con=fconn)

    for (uname in names(pars$u))
    {
        writeLines(paste(uname,
                         pars$u[uname],
                         pars$u.SE[uname],
                         sep=FIELD_DELIM),
                   con=fconn)
    }
    close(fconn)
    cat('Model parameters written to', filename, '\n')

}

# save model parameters
printout('Writing model parameters to file')
trait_bpar_file <- file.path(trait_dir, paste0(trait, '.bpar')) # to save model parameters with trait results
paramfiles <- c(trait_bpar_file, bpar_file)


for (parfile in paramfiles) {

    write_pars2(
        trait = trait, 
        trait_file = train_pheno_file, 
        genotype_prefix = train_dat$geno_file, 
        pars = train_fit, 
        genotype_source = train_geno_input,
        trait_source = trait_source,
        trait_variable = variable_used,
        filename = parfile
    )
}

# save cross-validation results 
for (cv_outdir in c(trait_dir, trait_bpar_dir)) {
    save_cv_results(cv_out, cv_outdir)
    plot_kfold(cv_out, cv_outdir)
}


save_training_data <- function(bpar_file) {

    pars <- read_pars(bpar_file)
    
    # get variant IDs used in the model
    dat <- pars$data
    snp_ids <- rownames(dat)

    # get RFIDs and trait values used to train the model
    trait <- pars$meta$trait_name
    trait_var <- pars$meta$trait_source_variable
    trait_dat <- pars$meta$trait_file
    trait_dat <- read.csv(trait_dat)
    trait_dat <- trait_dat[,c('rfid',trait)]
    trait_dat <- trait_dat[complete.cases(trait_dat),]

    # write variants and pheno data to files
    outdir <- dirname(bpar_file)    
    snp_file <- file.path(outdir, paste0(trait_var, '_train_snps'))
    trait_file <- file.path(outdir, paste0(trait_var, '_train_pheno.csv'))
    
    writeLines(snp_ids, snp_file)
    write.table(trait_dat, trait_file, sep=',', row.names=F, col.names=F, quote=F, na='')

}

# save training variants and pheno data with the .bpar file
save_training_data(bpar_file = bpar_file)

printout(paste(trait, 'model fitting complete'))

