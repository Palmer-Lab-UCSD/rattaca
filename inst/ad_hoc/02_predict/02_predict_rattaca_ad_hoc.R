library(rattaca)
library(argparse)
library(viridis)
library(scales)

################################################
# ---------------  ARGUMENTS  --------------- ##
################################################

# create a parser closure for input arguments
parser <- argument_parser(
    argument_def('rattaca', type = 'character', default_val = 'rattaca',
        help = 'arbitrary positional argument included to enable argument parsing'),
    argument_def('--trait', type = 'character', required = TRUE, 
        help = 'trait name, as used in trait list and column headers'),
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
    argument_def('--bpar_file', type = 'character', required = TRUE,
        help = 'path to the bpar file w/ fitted BLUP estimates to use for prediction'),
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
                window by this number of SNPs')
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

# read bpar file: unless specified, use the file inside the trait directory
if (!exists('bpar_file')) {
    bpar_file <- file.path(trait_dir, paste0(trait, '.bpar'))
}
train_geno_file <- read_pars(bpar_file)$meta$genotypes_file
train_pheno_file <- read_pars(bpar_file)$meta$trait_file

# use_snps_file <- list.files(snp_dir, pattern=paste0('^', trait, '.*', 'prune.in$'), full.names=T)

cat('\n')
cat('Predictions for', trait)
cat('\n')
cat('trait:', trait, '\n')
cat('proj_dir:', proj_dir, '\n')
cat('generation:', generation, '\n')
cat('test_geno_input:', test_geno_input, '\n')
cat('train_geno_input:', train_geno_input, '\n')
cat('train_geno_file:', train_geno_file, '\n')
cat('train_pheno_file:', train_pheno_file, '\n')
# cat('use_snps_file:', use_snps_file, '\n')
cat('\n')


#################################################
## ---------------  DATA PREP  --------------- ##
#################################################

# get all RFIDs to use in the test set
test_rfids <- read.csv(rfid_list)$rfid
test_rfids_file <- format_ids_file(test_rfids, test_data_dir, paste0(trait, '_test_rfids'))

printout('Loading fitted model and training genotypes')

# read in the genotype training data
geno_train <- load_existing_plink_dataset(train_geno_file)

# read in the previously-fitted model
train_fit <- convert_bpar(bpar_file)

# ensure that train/test data use a common set of snps
get_snps_from_bpar <- function(bpar_file, outdir=NULL) {
    pars <- read_pars(bpar_file)
    dat <- pars$data
    snp_ids <- rownames(dat)
    out <- list(bpar_snps = snp_ids, bpar_snps_file = NULL)

    if (!is.null(outdir)) {
        file_prefix <- basename(pars$meta$plink_genotypes_prefix)
        filename <- paste0(file_prefix, '_from_bpar_snps')
        outfile <- file.path(outdir, filename)
        out$bpar_snps_file <- outfile
        writeLines(snp_ids, outfile)
    }
    return(out)
}

bpar_snps <- get_snps_from_bpar(bpar_file = bpar_file, outdir = snp_dir)

# make the test genotype dataset
printout('Producing the test genotype dataset')
geno_test <- make_plink_dataset(
    input_genotypes = test_geno_input, 
    output_dir = test_data_dir,
    outfile_prefix = trait,
    samples_to_keep = test_rfids_file,
    snps_to_keep = bpar_snps$bpar_snps_file,
    maf_cutoff = NULL,          
    missing_cutoff = NULL,  
    hwe_cutoff = NULL,      
    ld_prune = F,
    return_data = T
    )

get_common_snpset2 <- function(train_genotypes, test_genotypes) {
    
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
    
    return(all_snps)
}

# subset BLUPs to only those SNPs shared between train & test data
train_snps <- bpar_snps$bpar_snps
test_snps <- colnames(geno_test$geno)

if (!setequal(train_snps, test_snps)) {
    use_snps <- get_common_snpset2(train_snps, test_snps)
    cat('Using', length(use_snps), 'of',length(train_snps), 'training SNPs for prediction \n')
    geno_train$geno <- geno_train$geno[,use_snps]
    geno_test$geno <- geno_test$geno[,use_snps]
    use_snps_idx <- which(names(train_fit$u) %in% use_snps)
    train_fit$u <- train_fit$u[use_snps_idx]
    train_fit$u.SE <- train_fit$u.SE[use_snps_idx]
}


# #################################################
# ## --------------  PREDICTIONS  -------------- ##
# #################################################


# predict traits in the test set using marker effects from the fitted model
printout('Predicting traits and saving to file')
trait_preds <- predict_lmm(
    geno_test$geno, 
    train_fit$u, 
    train_fit$beta)

# save predictions and ranks to file
zscores <- get_ranks_zscores(
    predictions = trait_preds, 
    trait = trait, 
    output_dir = trait_dir)    

# plot PCA results for genotypes used in predictions
printout('Running and plotting genotype PCA')
trait_pca <- pca_plink_genotypes(
    geno_train$geno_file, 
    geno_test$geno_file, 
    trait)

plot_pca(
    trait_pca, 
    geno_train, 
    geno_test, 
    trait_dir)

# write parameters used for prediction to file
printout(paste('Saving', trait, paste0('gen', generation), 'BLUPs'))

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

if (setequal(train_snps, test_snps)) {
    n_snps <- length(train_snps)
    parfile <- paste0(trait, '_blups_used_gen', generation, '_n', n_snps, '.bpar')
} else {
    all_snps <- length(train_snps)
    n_snps <- length(use_snps)
    parfile <- paste0(trait, '_blups_used_gen', generation, '_n', n_snps, '_of_', all_snps, '.bpar')
}
parfile <- file.path(trait_dir, parfile)

bpar_md <- read_pars(bpar_file)$meta
train_geno_file <- bpar_md$genotypes_file
train_geno_source <- bpar_md$genotypes_source
train_pheno_file <- bpar_md$trait_file
train_pheno_source <- bpar_md$trait_source
variable_used <- bpar_md$trait_source_variable

write_pars2(
    trait = trait, 
    trait_file = train_pheno_file, 
    trait_source = train_pheno_source,
    trait_variable = variable_used,
    genotype_prefix = train_geno_file, 
    genotype_source = train_geno_source,
    pars = train_fit, 
    filename = parfile
)

printout(paste(trait, 'predictions complete!'))

# #################################################
# ## ------------  POWER ANALYSIS  ------------- ##
# #################################################

# cat('Running power analysis:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')

# power <- test_power_multi(
#     trait = trait,
#     genotypes = geno_test$geno,
#     predictions = trait_preds,
#     fitted_mod = use_fit,
#     reps = 10,
#     tests_per_rep = 1000,
#     outdir = trait_dir)
