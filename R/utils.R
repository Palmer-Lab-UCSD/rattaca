# Utility functions for RATTACA
#
# Author: Robert Vogel
# Date: 2023-08-20
#   
# Contributors:
#   
#

#' Test whether a matrix consists of genotype data
#'
#' @description
#' Valid genotype matrices should contain values in the
#' set \eqn{\{1,0,-1\}} if no imputation performed.  Otherwise,
#' genotype values \eqn{g_{ij} \in [-1, 1]}.
#'
#' @export
#'
#' @param x ((n, m) array | vector | matrix) of genotype
#'      values
#' @param imputed_vals (boolean) whether input contains
#'      imputed values
#'
#' @return (boolean)
#
is_genotype <- function(x, imputed_vals=TRUE)
{

    if (imputed_vals)
        return(all(x <= 1) && all(x >= -1))        

    return(length(setdiff(x, c(1, 0, -1))) == 0)
}


#' Load plink genotypes and perform mean imputation
#'
#' @export
#' 
#' @param prefix_name (character)
#'      prefix of plink bim, bed, and fam files that will be loaded
#'
#' @return ((n sample, q markers) matrix)
#'      of genotype values
load_and_prepare_plink_data <- function(prefix_name)
{
    # the genotype matrix with q marker rows and
    # n sample columns
    genotypes <- genio::read_plink(prefix_name)$X - 1

    return(impute(genotypes))
}


#' load trait data from csv and prepare
#'
#' @description
#'      Data are loaded assuming a header exists.  Using this 
#'      header, the row names are set to the sample id's and
#'      and samples with trait value equal to NA are removed.
#'
#' @export
#'
#' @param filename (character)
#'      path to trait csv table
#' @param id_column (character)
#'      name of column with animal IDs whose values will be the data index names
#' @param trait_name (character)
#'      name of column whose values consists of noramlized
#'      trait measurements
#' @param format (string)
#'      The data type in which to output results 
#'      (either 'data.frame' or 'numeric')
#' 
#' @return A dataframe with n_sample rows x 1 column, rows named with sample IDs, 
#'      or a numeric vector named with sample IDs.
#
prep_trait_data <- function(
    filename,
    id_column,
    trait_name,
    format = c('numeric','data.frame'))
{
    data <- utils::read.csv(filename)[,c(id_column, trait_name)]
    data <- data[!is.na(data[,trait_name]), ]
    rownames(data) <- data[, id_column]
    data[,id_column] <- NULL

    if (format == 'numeric') {
        ids <- rownames(data)
        data <- as.numeric(data[,trait_name])
        names(data) <- ids
    }

    return(data)
}


# TO DO: when processing multiple traits, exclude any rows 
# that have NAs for all traits

#' Save a Plink-formatted file listing all samples phenotyped
#' for one or all traits in a dataset.
#'
#' @export
#'
#' @param phenotype_data (string)
#'      File path to a csv file containing one named column of 
#'      sample IDs and one or more named columns with phenotype
#'      data.
#' 
#' @param id_column (string)
#'      Text string with the header for the column containing
#'      sample IDs.
#' 
#' @param trait_column (string)
#'      (default NULL) Text string with the header for the desired
#'      phenotype, if only processing one. If NULL, all non-ID
#'      columns in the file will be processed. 
#' @param output_dir (string)
#'      The directory in which the Plink-formatted IDs file will be 
#'      saved.
#'
#' @return A list containing (1) the file path to the Plink-formatted 
#'      IDs file and (2) a vector of phenotyped IDs.
#
get_phenotyped_ids <- function(phenotype_data, id_column, output_dir, trait_column = NULL) {

    if (class(phenotype_data) == 'data.frame') {
        pheno_dat <- phenotype_data
    } else {
        pheno_dat <- read.csv(phenotype_data)
    }

    if (is.null(trait_column)) {
    
        phtyped_ids_file <- file.path(output_dir, 'phtyped_ids_all')
    
        write.table(data.frame(fam = 0, id = pheno_dat[[id_column]]), 
                phtyped_ids_file, sep='\t', row.names=F, col.names=F, quote=F)

        return(list(ids_file = phtyped_ids_file, ids = as.character(pheno_dat[[id_column]])))

    } else {
    
        trait <- trait_column
        phtyped_ids_file <- file.path(output_dir, paste0('phtyped_ids_', trait))
    
        trait_dat <- pheno_dat[[trait]]
        names(trait_dat) <- pheno_dat[[id_column]]
        trait_dat <- trait_dat[!is.na(trait_dat)]
        
        write.table(data.frame(fam = 0, id = names(trait_dat)), 
                phtyped_ids_file, sep='\t', row.names=F, col.names=F, quote=F)

        return(list(ids_file = phtyped_ids_file, ids = as.character(names(trait_dat))))

        
    }
}

#' Save a Plink-formatted file listing all samples input in a vector.
#'
#' @export
#' 
#' @param ids (character)
#'      A vector of sample IDs.
#' @param output_dir (string)
#'      The directory in which the Plink-formatted IDs file will be 
#'      saved.
#' @param filename (string)
#'      (default 'plink_ids') The desired file name
#'
#' @return A list containing (1) the file path to the Plink-formatted 
#'      IDs file and (2) a vector of IDs.
#
format_ids_file <- function(ids, output_dir, filename='plink_ids') {
    ids_file <- file.path(output_dir, filename)
    write.table(data.frame(fam = 0, id = ids), ids_file,
        sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    return(ids_file)
}

#' Modify a Plink genotype dataset to produce a new Plink dataset
#' of system files, with option to read data into R
#'
#' @export
#'
#' @param input_genotypes (string)
#'      prefix of plink bim, bed, and fam files that will be
#'      loaded and modified
#' @param output_dir (string)
#'      Output directory in which to save the new Plink dataset
#' @param outfile_prefix (string) 
#'      (default '') The base file name for new Plink files.
#'      This will be appended with modifications to the original 
#'      dataset.
#' @param samples_to_keep (string) 
#'      (default NULL) If subsetting samples from the original
#'      dataset, the path to the file listing all samples desired
#'      for the new dataset. Must be in Plink .fam format, with 
#'      one row per sample and two columns: one column of family IDs 
#'      (or zeros), and one column of individual IDs, as produced
#'      by get_phenotyped_ids(). Using this option calls the Plink
#'      option '--keep'.
#' @param snps_to_keep (string) 
#'      (default NULL) If subsetting SNP variants from the original
#'      dataset, the path to the file listing all SNP variants desired
#'      for the new dataset. Must be formatted with one column of 
#'      variant IDs, one row per variant, as produced by either 
#'      sample_snps() or sample_snps_from_plink_files(). Using this 
#'      option calls the Plink option '--extract'.
#' @param maf_cutoff (double) 
#'      (default NULL) The desired minor allele frequency cutoff. All 
#'      variants with allele frequency in the original dataset below 
#'      the provided threshold will be excluded from the new dataset.
#'      Using this option calls the Plink option '--maf'.
#' @param missing_cutoff (double) 
#'      (default NULL) The desired missing call rate cutoff. All variants 
#'      with missing call rates exceeding the provided threshold in the 
#'      original dataset will be excluded from the new dataset. Using 
#'      this option calls the Plink option '--geno'.
#' @param hwe_cutoff (double) 
#'      (default NULL) The desired cutoff for Hardy-Weinberg exact test 
#'      P-values. All variants in the original dataset whose P-values 
#'      following a Hardy-Weinberg equilibrium exact test fall below the 
#'      provided threshold will be excluded from the new dataset. Using
#'      this option calls the Plink option '--hwe'.
#' @param ld_prune (bool) 
#'      (default FALSE) Whether or not to prune variants from the original
#'      dataset based on linkage disequilibrium with other variants. Using
#'      this option calls the Plink option '--indep-pairwise'.
#' @param ld_r2 (double) 
#'      (default NULL) If LD-pruning, The desired pairwise r^2 threshold.
#'      For all variant pairs within a desired window size, if the r^2
#'      of the pairwise correlation between variant allele counts exceeds
#'      this threshold, only one variant will be kept in the new dataset.
#' @param ld_window_size (integer) 
#'      (default NULL) If LD-pruning, the desired window size (in variant 
#'      count) in which to estimate pairwise linkage disequilibrium.
#' @param ld_step_size (integer) 
#'      (default NULL) If LD-pruning, the desired step size (in variant 
#'      count) by which the pruning window will slide before estimating 
#'      pairwise LD again.
#' @param snp_directory (string)
#'      (default NULL) If LD-pruning, the desired directory path in which to
#'      save pruned.in and pruned.out files
#' @param return_data (bool) 
#'      (default TRUE) Whether or not to return the new dataset in the R 
#'      environment. If TRUE, returns a list containing (1) the file 
#'      path to the new dataset and (2) the new genotype matrix. If 
#'      FALSE, returns only the file directory. This option is preferable 
#'      when producing large datasets that become computationally intractable 
#'      to process within R.
#'
#' @return A list containing (1) the file path to the new dataset and 
#'      (2) the new genotype matrix.
#
make_plink_dataset <- function(input_genotypes,       
                               output_dir,            
                               outfile_prefix = '',   
                               samples_to_keep = NULL,
                               snps_to_keep = NULL,   
                               maf_cutoff = NULL,     
                               missing_cutoff = NULL, 
                               hwe_cutoff = NULL,     
                               ld_prune = FALSE,      
                               ld_r2 = NULL,          
                               ld_window_size = NULL, 
                               ld_step_size = NULL,
                               snp_directory = NULL,   
                               return_data=TRUE)      
{
    
    # common arguments
    args <- c('-bfile', input_genotypes, '-make-bed')
    
    # set extra arguments based on user input
    
    # samples to keep
    if (!is.null(samples_to_keep)) {
        n_samples <- length(readLines(samples_to_keep))
        args <- c(args, '--keep', samples_to_keep)
        outfile_prefix <- paste0(outfile_prefix, '_n', n_samples)
    }
    # minor allele frequency cutoff
    if (!is.null(maf_cutoff)) {
        args <- c(args, '--maf', maf_cutoff)
        outfile_prefix <- paste0(outfile_prefix, '_maf', maf_cutoff)
    }
    
    # genotype missingness cutoff
    if (!is.null(missing_cutoff)) {
        args <- c(args, '--geno', missing_cutoff)
        outfile_prefix <- paste0(outfile_prefix, '_missing', missing_cutoff)
    }
    
    # HWE cutoff
    if (!is.null(hwe_cutoff)) {
        args <- c(args, '--hwe', hwe_cutoff)
        hwe_pval <- -log10(hwe_cutoff)
        outfile_prefix <- paste0(outfile_prefix, '_hwe', hwe_pval)
    }

    # SNP sampling
    if (!is.null(snps_to_keep)) {
        args <- c(args, '--extract', snps_to_keep)
        snp_n <- length(readLines(snps_to_keep)) / 1000
        outfile_prefix <- paste0(outfile_prefix, '_', snp_n, 'k_snps')
    }

    # LD pruning
    if (ld_prune) {
        
        if( is.null(ld_r2) | is.null(ld_window_size) | is.null(ld_step_size) | is.null(snp_directory)) {
            stop('Please provide LD-pruning parameters: r^2, window size, step size, and snp directory')
        }
        
        args <- c(args, '--indep-pairwise', ld_window_size, ld_step_size, ld_r2)
        outfile_prefix <- paste(outfile_prefix, 'ldprune', ld_r2, ld_window_size, ld_step_size, sep='_')
    }
    
    # set plink arguments
    plink_file_name <- file.path(output_dir, outfile_prefix)
    args <- c(args, '--out', plink_file_name)

    # print the plink call to the user
    cat('Plink call:', plink2, paste(args, collapse=' '), '\n')

    # run plink
    system2(plink2, args)

    # if LD-pruning, create the final dataset using only SNPs in linkage equilibrium
    if (ld_prune) {

        # reset arguments to read in the file that was just produced
        args <- c('-bfile', plink_file_name, '-make-bed')

        snps_to_keep <- paste0(plink_file_name,'.prune.in')
        snp_n <- length(readLines(snps_to_keep))
        ld_outfile_prefix <- paste0(outfile_prefix, '_', snp_n, '_LDpruned_snps')
        ld_plink_file_name <- file.path(output_dir, ld_outfile_prefix)

        args <- c(args, '--extract', snps_to_keep, '--out', ld_plink_file_name)

        # print the plink call to the user
        cat('Plink call:', plink2, paste(args, collapse=' '), '\n')

        # run plink
        system2(plink2, args)

        # remove the intermediate plink files
        system2('rm', args=c(paste0(plink_file_name,'.bed'), 
            paste0(plink_file_name,'.bim'), paste0(plink_file_name,'.fam')))
        
        # move LD-pruned SNP files to SNP directory
        system2('cp', args=c(snps_to_keep, snp_directory))

        # reset the plink file name to output to the R console
        plink_file_name <- ld_plink_file_name
    }
    
    if (return_data){
    
        # read in plink genotypes and replace NA values
        plink_dat <- load_and_prepare_plink_data(plink_file_name)
        out <- list(geno_file = plink_file_name, geno = plink_dat)
        
        if (ld_prune) {
            out$ldpruned_snps_file <- snps_to_keep
        }

    } else {
        out <- list(geno_file = plink_file_name, geno = NULL)
    }
    
    return(out)
}


#' Produce a rattaca geno data object from an existing Plink genotype dataset.
#'
#' @export
#'
#' @param file_stem (string)
#'      The file path/prefix for the Plink dataset.
#' 
#' @return A rattaca genotype dataset: a list of (1) the file stem for the 
#'      dataset and (2) the genotype matrix.
#
load_existing_plink_dataset <- function(file_stem) {

    plink_dat <- load_and_prepare_plink_data(file_stem)
    out <- list(geno_file = file_stem, geno =  plink_dat)
}


#' Identify the set of SNP variants common to both the training
#' and test datasets
#'
#' @export
#'
#' @param train_genotypes (string, list, or matrix)
#'      Alternatively, the file path/prefix for the training Plink dataset,
#'      a training genotype data object as produced by make_plink_dataset(), 
#'      or a named training genotype matrix with RFID rows and variant columns
#' 
#' @param test_genotypes (string, list, or matrix)
#'      Alternatively, the file path/prefix for the test Plink dataset,
#'      a test genotype data object as produced by make_plink_dataset(), 
#'      or a named test genotype matrix with RFID rows and variant columns
#'
#' @return A vector of all SNP variants found in both datasets.
#
get_common_snpset <- function(train_genotypes, test_genotypes) {

    # read in snp sets
    if (is.character(train_genotypes)) {
        all_train_snps <- genio::read_bim(paste0(train_genotypes, '.bim'))
        all_train_snps <- all_train_snps$id
    } else if ('geno_file' %in% names(train_genotypes)) {
        all_train_snps <- colnames(train_genotypes$geno)
    } else if (is.matrix(train_genotypes)) {
        all_train_snps <- colnames(train_genotypes)
    } else {
        cat('Check train_genotypes format: \n')
        print(str(train_genotypes))
    }   

    if (is.character(test_genotypes)) {
        all_test_snps <- genio::read_bim(paste0(test_genotypes, '.bim'))
        all_test_snps <- all_test_snps$id
    } else if ('geno_file' %in% names(test_genotypes)) {
        all_test_snps <- colnames(test_genotypes$geno)
    } else if (is.matrix(test_genotypes)) {
        all_test_snps <- colnames(test_genotypes)
    } else {
        cat('Check test_genotypes format: \n')
        print(str(test_genotypes))
    }


#' Produce one or multiple sets of SNPs randomly sampled
#' from an input SNP set
#'
#' @export
#'
#' @param input_snps (character)
#'      A vector of SNP variant names from which to sample. Names must be
#'      in format 'chromosome:position', e.g. 1:2456462
#' 
#' @param keep (int)
#'      (default 50000) The number of SNPs to keep
#'
#' @param iterations (int)
#'      (default 1) The number of random samples to produce      
#'
#' @param save_file (bool)
#'      (default FALSE) Whether or not to save the SNP sample(s) to file(s)
#'
#' @param output_dir (character)
#'      The output directory in which to save the SNP sample(s)
#'
#' @return A list of length 'iterations' with multiple random SNP samples
#
sample_snps <- function(input_snps,
                        keep = 50000,
                        iterations = 1,
                        save_file = FALSE,
                        output_dir = NULL){
    
    # ensure inputs/outputs are provided if files are desired
    if (save_file) {
        if (is.null(output_dir)){
            stop('Include an output directory in which to save your file(s)')
        }
    }
    
    # list to store SNP samples
    snp_list <- list()
    
    n_k <- paste0(keep/1000,'k')
    max_i_figs <- nchar(iterations)

    for (i in 1:iterations){
                
        # sample SNPs, save SNPs to a dataframe
        use_snps <- sample(input_snps, keep)
        use_snps <- strsplit(use_snps,':')
        chr <- sapply(use_snps, function(x) x[1])
        pos <- sapply(use_snps, function(x) x[2])
        snp_df <- data.frame(cbind(chr,pos))
        snp_df$pos <- as.numeric(snp_df$pos)
        extra_chrs <- c('MT', 'X', 'Y')
        snp_num <- snp_df[!snp_df$chr %in% extra_chrs,]
        snp_char <- snp_df[snp_df$chr %in% extra_chrs,]
        snp_num$chr <- as.numeric(snp_num$chr)
        snp_char$chr <- as.character(snp_char$chr)
        snp_num <- snp_num[order(snp_num$chr, snp_num$pos),]
        snp_char <- snp_char[order(snp_char$chr, snp_char$pos),]
        snp_df <- rbind(snp_num, snp_char)
        use_snps <- paste0(snp_df$chr, ':', snp_df$pos)

        snp_list[[i]] <- use_snps
    
        if (save_file){
            
            i_figs <- nchar(i)
            rep_0s <- paste0(rep(0, max_i_figs - i_figs),collapse='')
            write.table(use_snps, file.path(output_dir, paste0('snpset_', n_k, '_random_', rep_0s, i)), 
                        quote=F, row.names=F, col.names=F)    
        }

    }

    return(snp_list)          
}


#' Store model parameters from a bpar file.
#' 
#' @description
#' Produce one sample of n snps randomly sampled from a plink bim file.
#'
#' @export
#'
#' @param n_snps (int)
#'      The desired number of variants to sample.
#' 
#' @param in_prefix (string)
#'      The path and base filename (without extension) of the Plink dataset 
#'      from which to sample snps.
#' 
#' @param outdir (string)
#'      The directory path in which to save a Plink-formatted file of sampled 
#'      variants.
#' 
#' @return A list with (1) the path to the sampled variants file and (2) the 
#'      vector of sampled variants.
#
sample_plink_snps <- function(
    n_snps,
    in_prefix,
    outdir) 
{
    bim_dat <- read.table(paste0(in_prefix,'.bim'))
    all_snps <- bim_dat$V2
    
    use_snps <- sort(sample(all_snps, n_snps))
    use_snps <- strsplit(use_snps,':')
    chr <- sapply(use_snps, function(x) x[1])
    pos <- sapply(use_snps, function(x) x[2])
    snp_df <- data.frame(cbind(chr,pos))
    snp_df$pos <- as.numeric(snp_df$pos)
    extra_chrs <- c('MT', 'X', 'Y')
    snp_num <- snp_df[!snp_df$chr %in% extra_chrs,]
    snp_char <- snp_df[snp_df$chr %in% extra_chrs,]
    snp_num$chr <- as.numeric(snp_num$chr)
    snp_char$chr <- as.character(snp_char$chr)
    snp_num <- snp_num[order(snp_num$chr, snp_num$pos),]
    snp_char <- snp_char[order(snp_char$chr, snp_char$pos),]
    snp_df <- rbind(snp_num, snp_char)
    use_snps <- paste0(snp_df$chr, ':', snp_df$pos)

    n_k <- paste0(n_snps/1000,'k')
    filename <- paste0('snpset_', n_k, '_random')
    outfile <- file.path(outdir, filename)

    writeLines(use_snps, outfile)

    return(list(snp_file = outfile, snps = use_snps))
}


#' Produce one or multiple sets of SNPs randomly sampled from an input
#' SNP set read in from Plink files, then use the new sample(s) to produce
#' a new Plink dataset(s)
#'
#' @export
#'
#' @param input_genotypes (character)
#'      The path/prefix of Plink files from which to sample SNPs
#' 
#' @param snp_directory (character)
#'      The directory storing 1 or more text files listing sampled SNPS,
#'      as produced by sample_snps()
#'
#' @param output_dir (character)
#'      The output directory in which to save the new Plink dataset(s)  
#'
#' @param n_samples (int)
#'      (default 1) The number of random samples and Plink datasets 
#'      to produce
#'
#' @param keep_files (bool)
#'      (default TRUE) Whether or not to keep Plink genotype files.
#'      Use FALSE if only SNP files are desired
#'
#' @return A list of length n_samples of genotype datasets as produced
#'          by make_plink_dataset()
#
sample_snps_from_plink_files <- function(input_genotypes,   # genotype data in plink format: base filename
                                         snp_directory,     # directory of 1 or more txt files listing sampled SNPs
                                         output_dir,        # directory to hold the new datasets
                                         n_samples = 1,     # number of SNP samples to extract & datasets to produce
                                         keep_files = TRUE) # T: keep plink files; F: delete plink files
{
    input_filename <- basename(input_genotypes)
    all_samples <- list.files(snp_directory, full.names=T)
    geno_datasets <- list()
    
    for (i in 1:n_samples) {
        
        snp_sample <- all_samples[i]
        
        geno_datasets[[i]] <- make_plink_dataset(
            input_genotypes = input_genotypes,
            output_dir = output_dir,             # directory for the new dataset
            outfile_prefix = paste0(input_filename, '_', i),    # base filename for the new dataset
            snps_to_keep = snp_sample,    # text file of snps to keep from the input dataset
            return_data=TRUE)       # T: return the dataset in R; F: just produce the data files

        # delete unwanted files
        if (!keep_files) {
            
            system_call <- paste0('rm', ' ', geno_datasets[[i]]$geno_file, '.*')
            system(system_call)
            cat('System call:', system_call, '\n')
        }
    }
    
    return(geno_datasets)
}


#' Align genotype and phenotype datasets for prediction. 
#' 
#' @description 
#' Keeps only samples that are shared between datasets, maintaining the same 
#' order of sample IDs in both
#'
#' @export
#'
#' @param genotypes (list)
#'      A genotype dataset as produced by make_plink_dataset() or
#'      sample_snps_from_plink_files(): a list with elements $geno_file
#'      (the path/prefix to the Plink dataset) and $geno (a genotype matrix
#'      of n samples x q variants)
#' 
#' @param phenotypes (numeric)
#'      A named numeric vector of phenotype measurements
#' 
#' @param trait (character)
#'      The name of the trait being analyzed
#'
#' @return A list of (1) the trait name, (2) all sample IDs shared between
#'      datasets, (3) the path/prefix for the Plink dataset used in
#'      alignment, (4) the aligned genotype matrix, and (5) the aligned
#'      phenotype data
#
align_data <- function(genotypes,
                       phenotypes,
                       trait)
{
    
    geno_file <- genotypes$geno_file
    geno <- genotypes$geno
        
    # align genotype and phenotype data for prediction
    rat_ids <- intersect(rownames(geno), names(phenotypes))
    geno <- geno[rat_ids,]
    trait_dat <- phenotypes[rat_ids]
    names(trait_dat) <- rat_ids

    out <- list(trait = trait, ids = rat_ids, geno_file = geno_file, geno = geno, pheno = trait_dat)
    return(out)
} 


#' Format trait observations and predictions into a dataframe
#'
#' @export
#'
#' @param lst (list)
#'      A list with elements $obs (named trait observations) and $pred (trait
#'      predictions in the same order as $obs). Each element may itself be a 
#'      list of length k sets of observations or predictions, such as output
#'      by kfold_cv()
#' 
#' @return A dataframe of all IDs, observations, and predictions, concatenated
#'      across all k folds
#
format_obs_pred <- function(lst) # list with trait observations and predictions
{

    k <- numeric()
    rfids <- character()
    obs <- numeric()
    pred <- numeric()
    
    for (i in seq_along(lst)) {
    
        k <- c(k, rep(i, length(lst[[i]]$obs)))
        rfids <- c(rfids, names(lst[[i]]$obs))
        obs <- c(obs, lst[[i]]$obs)
        pred <- c(pred, lst[[i]]$pred)
    
    }
    
    df <- data.frame(kfold = k, rfid = rfids, obs = obs, pred = pred)
    return(df)
        
}


#' Given a set of trait predictions, get the rank order and the z-score
#' of each
#'
#' @export
#'
#' @param predictions (numeric)
#'      A named vector of trait predictions
#' 
#' @param trait (character)
#'      The name of the trait to be analyzed
#' 
#' #' @param output_dir (character)
#'      (default NULL) The directory in which to save a dataframe of 
#'      trait predictions, prediction ranks, and z-scores, if desired
#' 
#' @return A list containing (1) the trait analyzed, (2) all prediction ranks,
#'      and (3) all prediction z-scores
#
get_ranks_zscores <- function(predictions, # named vector of trait predictions
                      trait,       # character string
                      output_dir=NULL) # character directory path
{

    pred_rank <- rank(predictions)
    pred_zscore <- zscore(predictions)
    
    if (!is.null(output_dir)){
        out_df <- data.frame(rfid = names(predictions), pred = predictions, 
                             rank = pred_rank, zscore = pred_zscore)
        names(out_df) <- c('rfid', trait, paste0(trait, '_rank'), paste0(trait, '_zscore'))
        write.csv(out_df, paste0(output_dir, '/', trait, '_predictions.csv'),
                             row.names=F, quote=F)

    }
                       
    return(list(trait = trait, rank = pred_rank, z_score = pred_zscore))    

}


#' Plot the results of a Plink PCA to png files
#'
#' @export
#'
#' @param plink_pca (list)
#'      A list containing elements $trait (the trait to analyze) and $pca_file
#'      (the path to Plink .eigenvec PCA results), as output by 
#'      pca_plink_genotypes()
#' 
#' @param train_genotypes (character)
#'      A genotype matrix from the training set
#' 
#' #' @param test_genotypes (character)
#'      A genotype matrix from the test set
#' 
#' #' @param output_dir (character)
#'      The directory in which to save png plots 
#' 
#' @return None. Plots are saved to files, without output into the R 
#'      environment
#
plot_pca <- function(plink_pca, 
                     train_genotypes, 
                     test_genotypes,
                     output_dir)
{

    # read PCA results
    pca_df <- read.table(paste0(plink_pca$pca_file, '.eigenvec'), header=F)
    eigenval <- read.table(paste0(plink_pca$pca_file, '.eigenval'), header=F)$V1
    names(eigenval) <- paste0('pc', seq(1,10))
    names(pca_df) <- c('fid', 'rfid', paste0('pc', seq(1,10)))
    train_rfids <- rownames(train_genotypes$geno)
    test_rfids <- rownames(test_genotypes$geno)
    pca_df$dataset <- sapply(pca_df$rfid, function(x) 
        if  (x %in% train_rfids) {
            return('train')
        } else if (x %in% test_rfids) {
            return('test')
        } else {
            return(NA)
        }
        )

    file_name <- file.path(output_dir, paste0(plink_pca$trait, '_gtypes_pca.png'))
    png(file_name, width=10, height=8, units='in', res=300)
    par(mfrow=c(2, 3), oma=c(0, 0.3, 2.2, 0))

    plot(pca_df$pc1,pca_df$pc2, col=0,
         xlab=paste0('PC1  (', round(eigenval['pc1']/sum(eigenval)*100,2),' %)'),
         ylab=paste0('PC2  (', round(eigenval['pc2']/sum(eigenval)*100,2),' %)'))
    title(line=0.9, main='PC1 vs. PC2')
    points(pca_df[pca_df$dataset=='train',]$pc1,pca_df[pca_df$dataset=='train',]$pc2, col=alpha('blue',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='train',]$pc1,pca_df[pca_df$dataset=='train',]$pc2, lwd=1, cex=1)
    points(pca_df[pca_df$dataset=='test',]$pc1,pca_df[pca_df$dataset=='test',]$pc2, col=alpha('red',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='test',]$pc1,pca_df[pca_df$dataset=='test',]$pc2, lwd=1, cex=1)
    legend('topleft',inset=c(0.02,0.02), title='Dataset',legend=c('Train','Test'),bg='white',
           pch=16,cex=1,col=c(alpha('blue',0.5),alpha('red',0.5)))
    legend('topleft',inset=c(0.02,0.02), title='',legend=c('Train','Test'),bty='n',pch=1,pt.lwd=1,cex=1)

    plot(pca_df$pc1,pca_df$pc3, col=0,
         xlab=paste0('PC1  (', round(eigenval['pc1']/sum(eigenval)*100,2),' %)'),
         ylab=paste0('PC3  (', round(eigenval['pc3']/sum(eigenval)*100,2),' %)'))
    title(line=0.9, main='PC1 vs. PC3')
    points(pca_df[pca_df$dataset=='train',]$pc1,pca_df[pca_df$dataset=='train',]$pc3, col=alpha('blue',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='train',]$pc1,pca_df[pca_df$dataset=='train',]$pc3, lwd=1, cex=1)
    points(pca_df[pca_df$dataset=='test',]$pc1,pca_df[pca_df$dataset=='test',]$pc3, col=alpha('red',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='test',]$pc1,pca_df[pca_df$dataset=='test',]$pc3, lwd=1, cex=1)

    plot(pca_df$pc1,pca_df$pc4, col=0,
         xlab=paste0('PC1  (', round(eigenval['pc1']/sum(eigenval)*100,2),' %)'),
         ylab=paste0('PC4  (', round(eigenval['pc4']/sum(eigenval)*100,2),' %)'))
    title(line=0.9, main='PC1 vs. PC4')
    points(pca_df[pca_df$dataset=='train',]$pc1,pca_df[pca_df$dataset=='train',]$pc4, col=alpha('blue',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='train',]$pc1,pca_df[pca_df$dataset=='train',]$pc4, lwd=1, cex=1)
    points(pca_df[pca_df$dataset=='test',]$pc1,pca_df[pca_df$dataset=='test',]$pc4, col=alpha('red',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='test',]$pc1,pca_df[pca_df$dataset=='test',]$pc4, lwd=1, cex=1)

    plot(pca_df$pc2,pca_df$pc3, col=0,
         xlab=paste0('PC2  (', round(eigenval['pc2']/sum(eigenval)*100,2),' %)'),
         ylab=paste0('PC3  (', round(eigenval['pc3']/sum(eigenval)*100,2),' %)'))
    title(line=0.9, main='PC2 vs. PC3')
    points(pca_df[pca_df$dataset=='train',]$pc2,pca_df[pca_df$dataset=='train',]$pc3, col=alpha('blue',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='train',]$pc2,pca_df[pca_df$dataset=='train',]$pc3, lwd=1, cex=1)
    points(pca_df[pca_df$dataset=='test',]$pc2,pca_df[pca_df$dataset=='test',]$pc3, col=alpha('red',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='test',]$pc2,pca_df[pca_df$dataset=='test',]$pc3, lwd=1, cex=1)

    plot(pca_df$pc2,pca_df$pc4, col=0,
         xlab=paste0('PC2  (', round(eigenval['pc2']/sum(eigenval)*100,2),' %)'),
         ylab=paste0('PC4  (', round(eigenval['pc4']/sum(eigenval)*100,2),' %)'))
    title(line=0.9, main='PC2 vs. PC4')
    points(pca_df[pca_df$dataset=='train',]$pc2,pca_df[pca_df$dataset=='train',]$pc4, col=alpha('blue',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='train',]$pc2,pca_df[pca_df$dataset=='train',]$pc4, lwd=1, cex=1)
    points(pca_df[pca_df$dataset=='test',]$pc2,pca_df[pca_df$dataset=='test',]$pc4, col=alpha('red',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='test',]$pc2,pca_df[pca_df$dataset=='test',]$pc4, lwd=1, cex=1)

    plot(pca_df$pc3,pca_df$pc4, col=0,
         xlab=paste0('PC3  (', round(eigenval['pc3']/sum(eigenval)*100,2),' %)'),
         ylab=paste0('PC4  (', round(eigenval['pc4']/sum(eigenval)*100,2),' %)'))
    title(line=0.9, main='PC3 vs. PC4')
    points(pca_df[pca_df$dataset=='train',]$pc3,pca_df[pca_df$dataset=='train',]$pc4, col=alpha('blue',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='train',]$pc3,pca_df[pca_df$dataset=='train',]$pc4, lwd=1, cex=1)
    points(pca_df[pca_df$dataset=='test',]$pc3,pca_df[pca_df$dataset=='test',]$pc4, col=alpha('red',0.5), cex=1, pch=16)
    points(pca_df[pca_df$dataset=='test',]$pc3,pca_df[pca_df$dataset=='test',]$pc4, lwd=1, cex=1)

    mtext(paste(plink_pca$trait, 'genotypes PCA'), outer=TRUE, font=2, cex=1.15, line=0)

    dev.off() 

    cat('PCA plots saved to file', file_name, '\n') 
}


#' Given results from one k-fold cross-validation, identify the
#' fold with the best model performance.
#'
#' @export
#'
#' @param cv_results (list)
#'      A list with results from one k-fold cross-validation, as output by 
#'      kfold_cv(). 
#' 
#' @param metric (character)
#'      (Default 'pearson') The desired statistic to use to identify the best 
#'      model fit. The 'best' fold is defined as that whose performance 
#'      maximizes this metric.
#' 
#' @return The integer value identifying the 'test' list element (i.e., the 
#'      fitted fold) with the best model performance, per the desired 
#'      performance metric.
#
best_fold <- function(cv_results,
                     metric='pearson') {

    test_results <- cv_results$test
    perf <- c()
    
    for (k in 1:length(test_results)) {
        if (metric == 'pearson') {
            perf <- c(perf, test_results[[k]]$pearson_corr)    
        } else if (metric == 'spearman') {
            perf <- c(perf, test_results[[k]]$spearman_corr)    
        } else if (metric == 'r_sq') {
            perf <- c(perf, test_results[[k]]$r_sq)                
        }
    }
    
    return(which.max(perf))
}


#' Given results from multiple k-fold cross-validations, identify the
#' model with the best mean cross-validation performance.
#'
#' @export
#'
#' @param crossval_list (list)
#'      A list with results from multiple k-fold cross-validations. Each
#'      element must be a list as produced by kfold_cv() tested on
#'      a different genomic dataset 
#' 
#' @return The integer value identifying the list element with the best
#'      mean cross-validation performance
#
best_kfold_mod <- function(crossval_list)
{
    if (length(crossval_list) == 1)
        stop('Must include results from >1 cross-validation run')
    
    mean_cv_perf <- numeric()
    
    # get the mean pearson correlation coef for each CV fold
    for (i in 1:length(crossval_list)){
    
        cv <- crossval_list[[i]]
        cv_test <- cv$test
        fold_r <- numeric()

        for (j in 1:length(cv_test)){
            
            fold_r <- c(fold_r, cv_test[[j]]$pearson_corr)
        }
        
        mean_r <- mean(fold_r)
        mean_cv_perf <- c(mean_cv_perf, mean_r)
    }
    
    # identify the fold with the best performance
    top_mod <- which.max(mean_cv_perf)
    
    return(top_mod)
}


#' Save the results from one k-fold cross-validation to a csv.
#'
#' @export
#'
#' @param cv_results (list)
#'      A list with results from one k-fold cross-validations, as output by 
#'      kfold_cv(). 
#' #' @param output_dir (character)
#'      The directory in which to save results
#' 
#' @return None. Results are saved to file
#
save_cv_results <- function(cv_results, output_dir) {
    
    test_results <- cv_results$test
    n_folds <- length(test_results)
    
    obs <- c()
    pred <- c()
    r_sq <- c()
    r <- c()
    rho <- c()
    fold <- c()

    for (k in 1:length(test_results)) {
        out <- cv_results$test[[k]]
        obs <- c(obs, out$obs)
        pred <- c(pred, out$pred)
        fold <- c(fold, rep(k, length(out$obs)))
        r_sq <- c(r_sq, rep(out$r_sq, length(out$obs)))
        r <- c(r, rep(out$pearson_corr, length(out$obs)))
        rho <- c(rho, rep(out$spearman_corr, length(out$obs)))
    }

    df1 <- data.frame(
        trait = rep(trait, length(obs)),
        fold = fold,
        obs = obs,
        pred = pred,
        r_sq = r_sq,
        r = r,
        rho = rho)

    df2 <- data.frame(
        trait = unique(df1$trait),
        fold = unique(df1$fold),
        r_sq = unique(df1$r_sq),
        r = unique(df1$r),
        rho = unique(df1$rho))

    write.csv(df1, file.path(output_dir, paste0(trait, '_', n_folds, 'fold_cv_results.csv')),
        row.names=F, quote=F, na='')
    write.csv(df2, file.path(output_dir, paste0(trait, '_', n_folds, 'fold_cv_summary.csv')),
        row.names=F, quote=F, na='')
}

#' Plot the results of a k-fold cross validation to files
#'
#' @export
#'
#' @param kfold_results (list)
#'      A list with results from a k-fold cross-validation, as produced
#'      by kfold_cv(). Must contain elements $trait (the trait name) and
#'      $test (a list of k sets oftrait observations, predictions, and 
#'      performance metrics from a k-fold cross-validation)
#' 
#' @param output_dir (string)
#'      The directory in which the png files will be saved
#' 
#' @return None. Plots are saved to files, without output into the R 
#'      environment
#
plot_kfold <- function(kfold_results, output_dir) {

    test_dat <- kfold_results$test
    num_folds <- length(test_dat)
    trait <- kfold_results$trait
    test_df <- format_obs_pred(test_dat)
    lm_slope <- numeric()
    lms <- list()
    
    # plot 1: zoomed to fit
    png(file.path(output_dir, paste0(trait, '_', num_folds, 'fold_crossval_zoomed.png')), width=7, height=7, units='in', res=300)

    # empty plot
    with(test_df, plot(obs, pred, col=0,
        xlab='Observed', ylab='Predicted'))
    title(main = paste0(trait, ' ', num_folds, '-fold cross-validation'),
         line = 2.2)
    plot_cols <- viridis(num_folds, 1, 0, 0.6, 1)
    
    # plot points
    for (i in 1:num_folds){
        
        plot_dat <- test_dat[[i]]
        plot_df <- test_df[test_df$kfold==i,]
        plot_col <- 'black' #viridis(0.1,0.1,0.1) # make this relate to i somehow
        plot_col <- plot_cols[i]
        lm <- lm(pred ~ obs, plot_df)
        lms[[i]] <- lm
        lm_slope <- c(lm_slope, lm$coef[2])
        with(plot_df, points(obs, pred, col = plot_col, pch = 16))
        with(plot_df, points(obs, pred))
        abline(lm, lwd = 2, col = plot_col)
    }
    
    # plot lines
    for ( i in 1:length(lms)){
        abline(lms[[i]], lty=2, col=plot_cols[i])
    }
    
    # legend
    legend('bottomright', title='Fold', legend=seq.int(1,num_folds), 
           bg='white', cex=0.8, lty=1, lwd=2, col=plot_cols)

    # margin text: mean performance metrics
    r_sq <- mean(sapply(test_dat, function(x) x$r_sq))
    r <- mean(sapply(test_dat, function(x) x$pearson_corr))
    rho <- mean(sapply(test_dat, function(x) x$spearman_corr))                 
    r_sq <- paste0('r_sq: ', round(r_sq, 3))
    r <- paste0('r: ', round(r,3))
    rho <- paste0("rho: ", round(rho,3))
    m <- paste0("m: ", round(mean(lm_slope),3))
    str <- paste('mean', paste(r_sq, r, rho, m, sep = "  |  "))
    # abline(lm, lwd =2)
    mtext(str,side=3,adj=0.05,line=0.2,cex=1.1)

    dev.off()
                       

    # plot 2: axes scaled 1:1
    trait_min <- min(c(test_df$obs, test_df$pred))
    trait_max <- max(c(test_df$obs, test_df$pred))

    png(file.path(output_dir, paste0(trait, '_', num_folds, 'fold_crossval_scaled.png')), width=7, height=7, units='in', res=300)

    with(test_df, plot(obs, pred, col=0, xlim=c(trait_min, trait_max), ylim=c(trait_min, trait_max),
        xlab='Observed', ylab='Predicted'))
    title(main = paste0(trait, ' ', num_folds, '-fold cross-validation'),
         line = 2.2)
    plot_cols <- viridis(num_folds, 1, 0, 0.6, 1)
    
    # plot points
    for (i in 1:num_folds){
        
        plot_dat <- test_dat[[i]]
        plot_df <- test_df[test_df$kfold==i,]
        plot_col <- 'black' #viridis(0.1,0.1,0.1) # make this relate to i somehow
        plot_col <- plot_cols[i]
        lm <- lm(pred ~ obs, plot_df)
        lms[[i]] <- lm
        lm_slope <- c(lm_slope, lm$coef[2])
        with(plot_df, points(obs, pred, col = plot_col, pch = 16))
        with(plot_df, points(obs, pred))
        abline(lm, lwd = 2, col = plot_col)
    }
    
    # plot lines
    for ( i in 1:length(lms)){
        abline(lms[[i]], lty=2, col=plot_cols[i])
    }
    
    # legend
    legend('bottomright', title='Fold', legend=seq.int(1,num_folds), 
           lty=1, lwd=2, col=plot_cols)

    # margin text: mean performance metrics
    r_sq <- mean(sapply(test_dat, function(x) x$r_sq))
    r <- mean(sapply(test_dat, function(x) x$pearson_corr))
    rho <- mean(sapply(test_dat, function(x) x$spearman_corr))                 
    r_sq <- paste0('r_sq: ', round(r_sq, 3))
    r <- paste0('r: ', round(r,3))
    rho <- paste0("rho: ", round(rho,3))
    m <- paste0("m: ", round(mean(lm_slope),3))
    str <- paste('mean', paste(r_sq, r, rho, m, sep = "  |  "))
    # abline(lm, lwd =2)
    mtext(str,side=3,adj=0.05,line=0.2,cex=1.1)

    dev.off()

}


#' Given a directory containing multiple csv files of predictions for 
#' different traits, merge all files into a single csv of results for
#' all traits
#'
#' @export
#'
#' @param directory (character)
#'      The directory path housing prediction results for multiple traits.
#'      Within this directory, each individual trait's results must be 
#'      inside a separate sub-directory identified by the trait name, with
#'      a csv file named as <trait_name>_predictions.csv (as output by 
#'      get_ranks_zscores()), e.g.: 
#'      <directory>/<trait_name>/<trait_name>_predicitons.csv
#' 
#' @param traits (list or character)
#'      Either a list or vector of trait names whose respective trait
#'      predictions will be merged
#' 
#' @param composite_traits (character)
#'      (default NULL) A vector of composite trait names whose respective 
#'      composite trait scores will be merged
#' 
#' @param output_dir (character)
#'      (default NULL) The directory in which the csv of merged trait
#'      predictions will be saved, if desired.
#' 
#' @param basename (character)
#'      (default NULL) The stem of the filename used for the merged csv
#'      file, if desired. Final results will be saved to
#'      <output_dir>/<basename>_merged_predictions.csv
#' 
#' @return The dataframe of all merged trait predictions
#
# function to merge all trait predictions into one file
merge_preds <- function(directory, # directory containing trait-named directories
                        traits,    # list or vector of trait names,
                        composite_traits=NULL, # list or vector of composite trait names
                        output_dir=NULL, # directory to save the file, if desired
                        basename=NULL)
{

    list_of_dfs <- list()
    
    # read in predictions for each trait
    for (trait in traits){
        trait_df <- read.csv(file.path(directory, trait, paste0(trait, '_predictions.csv')))
        list_of_dfs[[trait]] <- trait_df
    }

    # read in scores for each composite trait
    for (trait in composite_traits){
        trait_df <- read.csv(file.path(directory, trait, paste0(trait, '_composite_scores.csv')))
        list_of_dfs[[trait]] <- trait_df
    }
    
    # merge all dataframes
    all_preds <- Reduce(function(x,y) merge(x, y, all=T), list_of_dfs)
    
    if (!is.null(output_dir)){
        if (is.null(composite_traits)) {
            outfile <- file.path(output_dir, paste0(basename, '_merged_predictions.csv'))
        } else {
            outfile <- file.path(output_dir, paste0(basename, '_merged_preds_composite_scores.csv'))
        }
        write.csv(all_preds, outfile, row.names=F, quote=F, na='')
    }

    return(all_preds)
                        
}

#' Get assignment group designations (e.g., high vs low) for a set of 
#'      predictions.
#'
#' @export
#'
#' @param preds (list or dataframe)
#'      A list of trait predictions or dataframe with at least one column of 
#'      trait predictions.
#' 
#' @param col (character)
#'      (default NULL) If preds is a dataframe, the name of the dataframe column 
#'      on which to designate group assignments.
#' 
#' @param n_groups (int)
#'      (default 2) The number of desired groups to split assignments into. The
#'      default 2 splits into 'high' vs. 'low' groups, 3 splits into 
#'      'low'/'mid'/'high', and other numbers split predictions into n_groups 
#'      numeric quantiles.
#' 
#' @return A vector of group assignments in the same order as col.
#
trait_groups <- function(preds, col=NULL, n_groups=2) {

    # ensure n_groups is numeric
    n_groups = as.numeric(n_groups)

    # extract predictions
    if (is.data.frame(preds) && !is.null(col)) {
        preds <- preds[[col]]
    }
    idx <- !is.na(preds)

    trait_quantiles <- quantile(preds[idx], probs = seq(0, 1, by = 1/n_groups))
    if (n_groups==2) {
        groups <- sapply(preds, function(x) ifelse(x > median(preds), 'high', 'low'))
    } else if (n_groups==3) {
        groups <- cut(preds, breaks = trait_quantiles, labels=c('low','mid','high'), include.lowest=T)
    } else if (n_groups > 3) {
        groups <- cut(preds, breaks = trait_quantiles, labels=F, include.lowest=T)
    }
    return(groups)
}


#' Split a phenotypic dataset into train and test samples, and produce files
#' of RFIDs for either sample. 
#'
#' @export
#'
#' @param phtypes_file (character)
#'      The path to a csv file with phenotype data. The file must have an 'rfid'
#'      column.
#' 
#' @param trait (character)
#'      The name of the dataframe column (the trait name) on which to conduct
#'      the train/test split
#' 
#' @param train_split (float)
#'      (default 0.7) The desired frequency of samples to randomly include in 
#'      the training set. The remainder (1 - train_split) will be included in 
#'      the test set.
#' 
#' @param train_dir (character)
#'      The directory path in which to store splitted training data.
#' 
#' @param test_dir (character)
#'      The directory path in which to store splitted test data.
#' 
#' @return A list with file paths to train and test data, file paths to text 
#'      files listing RFIDs present in either sample, and dataframes of 
#'      training and test data.
#
train_test_split <- function(
    phtypes_file, 
    trait,
    train_split=0.7,
    train_dir,
    test_dir) 
{
    test_split = 1 - train_split
    pheno_dat <- read.csv(phtypes_file)
    trait_dat <- pheno_dat[,c('rfid',trait)]
    trait_dat <- trait_dat[complete.cases(trait_dat),]
    rownames(trait_dat) <- trait_dat$rfid
    trait_dat <- prep_trait_data(phtypes_file, 'rfid', trait, 'data.frame')

    # set up training data
    train_rfids <- sample(rownames(trait_dat), train_split * nrow(trait_dat))
    train_df <- trait_dat[train_rfids,, drop=FALSE]
    train_df$rfid <- rownames(train_df)
    train_df <- train_df[,c(2,1)]

    # set up test data
    test_rfids <- setdiff(rownames(trait_dat), train_rfids)
    test_df <- trait_dat[test_rfids,, drop=FALSE]
    test_df$rfid <- rownames(test_df)
    test_df <- test_df[,c(2,1)]

    # save data
    train_path <- file.path(train_dir, paste0(trait,'_train_', train_split, '.csv'))
    test_path <- file.path(test_dir, paste0(trait,'_test_', test_split, '.csv'))
    write.csv(train_df, train_path, row.names=F, quote=F, na='')
    write.csv(test_df, test_path, row.names=F, quote=F, na='')

    # save ID lists in Plink format to allow downstream genotype data filtering
    train_ids_df <- data.frame(fam = 0, id = train_rfids)
    test_ids_df <- data.frame(fam = 0, id = test_rfids)
    train_ids_file <- file.path(train_dir, paste0('phtyped_ids_', trait, '_train_', train_split))
    test_ids_file <- file.path(test_dir, paste0('phtyped_ids_', trait, '_test_', test_split))     
    write.table(train_ids_df, train_ids_file , sep='\t', row.names=F, col.names=F, quote=F)
    write.table(test_ids_df, test_ids_file , sep='\t', row.names=F, col.names=F, quote=F)

    return(list(train_file = train_path, test_file = test_path,
                train_ids_file = train_ids_file, test_ids_file = test_ids_file,
                train_ids = train_rfids, test_ids = test_rfids,
                train = train_df, test = test_df))
}


#' Store model parameters from a bpar file.
#' 
#' @description
#' Reads a bpar file and stores model parameters as a rattaca model object, 
#' formatted as output by fit().
#'
#' @export
#'
#' @param bpar_file (character)
#'      The path to the bpar file with desired parameter values.
#' 
#' @return A list with model parameters, BLUPs, and performance statistics, 
#'      as output by fit().
#
convert_bpar <- function(bpar_file) {
    
    pars <- read_pars(bpar_file)
    md <- pars$meta
    dat <- pars$data
    out <- list()

    out$Vu <- md$variance_u
    out$Ve <- md$variance_e
    out$beta <- md$intercept
    out$beta.SE <- md$intercept_se
    out$u <- dat[,1]
    out$u.SE <- dat[,2]
    out$LL <- md$log_likelihood
    out$r_sq <- md$r_sq
    out$pearson_corr <- md$pearson_corr
    out$spearman_corr <- md$spearman_corr

    return(out)
}


#' Produce a file summarizing predictions model performance.
#' 
#' @description
#' Reads prediction results and a csv of trait metadata from multiple trait-
#' specific directories to merge metadata with model performance statistics.
#'
#' @export
#'
#' @param traits (character)
#'      A vector of trait names(s) (corresponding to directory names with 
#'      prediction results) to incorporate into the summary
#' 
#' @param results_dir (string)
#'      The directory path to the results folder housing output directories 
#'      named per trait. The function will search for files within each of 
#'      the directories named in 'traits'.
#' 
#' @param basename (string)
#'      The base name for the output summary file. The final output will be 
#'      named <basename>_summary_<datestamp>.csv.
#' 
#' @param pheno_dict (string)
#'      (Default NULL) The path to a phenotype data dictionary csv. The file 
#'      must have columns 'trait', 'description', 'covariates', 'heritability', 
#'      'variable_used', and 'source_file'.
#' 
#' @return A list with model parameters, BLUPs, and performance statistics, 
#'      as output by fit().
#
summarize_preds <- function(
    traits,            # vector of trait names, as used for predictions
    results_dir,       # directory containing trait-named results folders 
    basename = 'all_traits',
    pheno_dict = NULL) # path to csv with trait, variable_used, description columns
{
    summary <- data.frame(
        trait = traits,
        n_train = NA,
        n_test = NA,
        mean_r_sq = NA,
        mean_r = NA,
        mean_rho = NA)

    if (!is.null(pheno_dict)) {
        pheno_dict <- read.csv(pheno_dict, quote = '"')
        pheno_dict$description <- gsub(',', ';', pheno_dict$description)
        pheno_dict$covariates <- gsub(',', '|', pheno_dict$covariates)
        summary$heritability <- NA
        summary$variable_used <- NA
        summary$description <- NA
        summary$covariates <- NA
        summary$source_file <- NA
    }

    for (i in 1:nrow(summary)) {
        trait <- summary$trait[i]

        train_fam_file <- list.files(file.path(results_dir, trait, 'train'), pattern = '.fam$', full.names = T)
        test_fam_file <- list.files(file.path(results_dir, trait, 'test'), pattern = '.fam$', full.names = T)
        cv_sum <- read.csv(list.files(file.path(results_dir, trait), pattern = 'cv_summary.csv$', full.names = T))
        
        summary$n_train[i] <- system(paste('cat', train_fam_file, '| wc -l'), intern = T)
        summary$n_test[i] <- system(paste('cat', test_fam_file, '| wc -l'), intern = T)
        summary$mean_r_sq[i] <- mean(cv_sum$r_sq)
        summary$mean_r[i] <- mean(cv_sum$r)
        summary$mean_rho[i] <- mean(cv_sum$rho)

        if (!is.null(pheno_dict)) {
            dict <- pheno_dict[pheno_dict$trait == trait,]
            summary$heritability[i] <- dict$heritability
            summary$variable_used[i] <- dict$variable_used
            summary$description[i] <- dict$description
            summary$covariates[i] <- dict$covariates
            summary$source_file[i] <- dict$source_file
        }
    }
    if (!is.null(pheno_dict)) {
        col_order <- c('trait','heritability','n_train','n_test','mean_r_sq','mean_r','mean_rho',
                       'description','covariates','variable_used','source_file')
        summary <- summary[,col_order]
    }
    datestamp <- format(Sys.time(), '%Y%m%d')
    write.csv(summary, file.path(results_dir, paste0(basename, '_summary_', datestamp, '.csv')), row.names=F, quote=F, na='')
    return(summary)
}


#' Produce a "composite" trait from the predictions of multiple individual 
#' traits, to be used in downstream assignment.
#' 
#' @description
#' Reads trait predictions and modifies them to negate or invert their ranks 
#' or z-scores, then averages traits per individual to produce the composite. 
#' This allows prioritizing multiple traits at once for assignment, and trait 
#' modifications allow prioritization of rats that are, say, intermediate for 
#' one trait or whose set of traits are most consistently extreme, whether 
#' positively or negatively correlated as desired.
#'
#' @export
#'
#' @param preds (dataframe or character)
#'      A dataframe or csv path to a file containing merged RATTACA predictions 
#'      for multiple traits.
#' 
#' @param traits (character)
#'      Vector of trait names to be included in calculating the composite trait.
#' 
#' @param negate (logical)
#'      A logical vector denoting, in the order of 'traits', whether to negate 
#'      the predicted values for each trait. Negating multiplies each prediction 
#'      by -1, effectively turning low ranks (or z-scores) high and high ranks
#'      low.
#' 
#' @param invert (logical)
#'      A logical vector denoting, in the order of 'traits', whether to invert 
#'      the predicted values for each trait. Inverting divides 1/prediction, 
#'      effectively turning intermediate z-scores extreme and extreme z-scores 
#'      intermediate, or large ranks small and small ranks large.
#' 
#' @param metric (string)
#'      Which prediction metric to use (either 'rank' or 'zscore') in 
#'      calculating the composite trait. 
#' 
#' @param stat (string)
#'      Which summary statistic to use (either 'mean' or 'median') in 
#'      calculating the composite trait. That is, the composite will be 
#'      calculated as either the mean or median of all individual traits 
#'      (per rat) after negation and/or inversion (or neither).
#' 
#' @param new_trait (string)
#'      The name of the new composite trait
#' 
#' @param plot_mar (numeric)
#'      (Defualt c(0,0,8,0)) A vector of plotting margins used to plot trait 
#'      correlations This may need tweaking for aesthetic plotting depending on 
#'      trait and variable names.
#' 
#' @param plot_oma (numeric)
#'      (Default c(0,0,2,2)) A vector of outer plotting margins used to plot 
#'      trait correlations This may need tweaking for aesthetic plotting 
#'      depending on trait and variable names.
#' @param title_line (numeric)
#'      (Default 7) The line above the correlation plot on which to print the
#'      plot title. This may need tweaking for aesthetic plotting.
#' @param output_dir (string)
#'      (Default NULL) The directory path in which to save output files
#' 
#' @return Returns nothing to the R console. Saves csv files for trait 
#'      correlations and composite trait values, plus png files plotting trait
#'      correlations.
#
make_composite_trait <- function(
    preds, # df or path to predictions csv
    traits, # vector of trait names desired to calculate a composite trait
    negate, # logical vector of length(traits): whether to multiply each trait by -1 to align desired high/low values (turns high to low, low to high)
    invert, # logical vector of length(traits): whether to perform 1/trait (turns intermediate values high/low and high/low values intermediate)
    metric = c('rank', 'zscore'), # metric to use in calculating the composite trait
    stat = c('mean', 'median'), # the statistic to use in calculating the composite trait
    new_trait, # new variable name for the composite trait
    plot_mar=c(0, 0, 8, 0), 
    plot_oma=c(0, 0, 2, 2),
    title_line=7,
    output_dir = NULL)
{
    library(corrplot) # for plotting composite trait correlations

    metric <- match.arg(metric)
    stat <- match.arg(stat)

    if (is.character(preds)) {
        preds <- read.csv(preds)
    }

    traits_df <- preds[,traits]

    # matrix to store data for relevant traits
    tmp_mat <- matrix(nrow = nrow(preds), ncol = length(traits))
    colnames(tmp_mat) <- traits
    rownames(tmp_mat) <- as.character(preds$rfid)

    # negate data as needed, save to matrix
    trait_dat <- list()
    new_names <- c()
    for (i in 1:length(traits)) {
        
        trait <- traits[i]
        trait_name <- trait
        needs_negation <- negate[i]
        needs_inversion <- invert[i]
        rank_col <- paste0(trait, '_rank')
        zscore_col <- paste0(trait, '_zscore')
        use_col <- ifelse(metric == 'rank', rank_col, zscore_col)

        trait_vals <- as.vector(preds[[use_col]])

        # negation: switches orientation of predictions
        # high ranks become low ranks, low ranks become high
        # negative zscores become positive, positive zscores become negative while maintaining magnitude
        if (needs_negation) {
            trait_name <- paste0(trait_name, '_negated')
            if (metric == 'rank') {
                max_pred <- max(trait_vals) 
                trait_vals <- max_pred + 1 - trait_vals
            } else if (metric == 'zscore') {
                trait_vals <- trait_vals * -1
            }
        }

        # inversion: switches priority from extreme to intermediate predictions
        # high and low ranks become intermediate, intermediate become high or low while maintaining orientation
        # small magnitude zscores become large, large magnitude zscores become small while maintaining high/low orientation
        if (needs_inversion) {
            trait_name <- paste0(trait_name, '_inverted')
            if (metric == 'rank') {
                median_pred <- median(trait_vals)
                trait_vals <- 1/(trait_vals - median_pred)
            } else if (metric == 'zscore') {
                trait_vals[which(trait_vals==0)] <- 1/1e9
                trait_vals <- 1/trait_vals           
            }
        }
    
        tmp_mat[,i] <- trait_vals
        new_names[i] <- trait_name

    } # end of traits loop

    # finalize modified trait names
    new_names <- paste0(new_names, '_', metric)
    colnames(tmp_mat) <- new_names
    colnames(traits_df) <- paste0(traits, '_', metric)

    # calculate the composite trait
    new_trait_vals <- apply(tmp_mat, 1, stat)

    new_trait_metrics <- get_ranks_zscores(
        predictions = new_trait_vals,
        trait = new_trait)

    # save new trait values, ranks, zscores to df
    new_vals_df <- data.frame(
        rfid = names(new_trait_vals),
        new_trait = new_trait_vals)
    names(new_vals_df)[2] <- new_trait

    new_metrics_df <- data.frame(
        rfid = names(new_trait_metrics$rank),
        rank = new_trait_metrics$rank,
        zscore = new_trait_metrics$z_score)
    names(new_metrics_df) <- c('rfid', paste0(new_trait, '_rank'), paste0(new_trait, '_zscore'))

    out_df <- merge(new_vals_df, new_metrics_df, by = 'rfid')

    # add the new trait to the traits dataframe/matrix
    # to plot trait correlations
    traits_df[[new_trait]] <- new_trait_vals
    traits_mat <- cbind(tmp_mat, new_trait_vals); colnames(traits_mat)[i+1] <- new_trait
    corr_df <- cor(traits_df, method='spearman')
    corr_mat <- cor(traits_mat, method='spearman')
    
    if (!is.null(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE)
        cat('Saving composite trait data to', paste0(output_dir,'/'), '\n')
        write.csv(out_df, file.path(output_dir, paste0(new_trait, '_composite_scores.csv')),
                  row.names=F, quote=F, na='')
        write.csv(corr_df, file.path(output_dir, paste0(new_trait, '_corr_raw_traits.csv')),
                  row.names=T, quote=F)
        write.csv(corr_mat, file.path(output_dir, paste0(new_trait, '_corr_altered_traits.csv')),
                  row.names=T, quote=F)
        png(file.path(output_dir, paste0(new_trait, '_corrplots.png')), width=11, height=7, units='in', res=300)
        par(mfrow=c(1,2), mar=plot_mar, oma=plot_oma) 
        corrplot(corr_df, type = 'upper', order = 'original', tl.col = 'black')
        mtext('raw sub-traits', side=3, line=title_line, cex=1.5, font=2)
        mtext(new_trait, side=3, line=title_line + 1.6, cex=1.5, font=2)
        corrplot(corr_mat, type = 'upper', order = 'original', tl.col = 'black')
        mtext('modified sub-traits', side=3, line=title_line, cex=1.5, font=2)
        mtext(new_trait, side=3, line=title_line + 1.6, cex=1.5, font=2)
        dev.off()
    }
    return(out_df)
}


#' Produce an assignment 'S' plot: trait rank vs. trait prediction.
#' 
#' @description
#' Reads a dataframe of predictions, ranks, and group assignments to plot 
#' prediction rank against predicted values, colored by assigned high/low group.
#'
#' @export
#'
#' @param df (dataframe)
#'      A dataframe of trait predictions for a trait of interest (as produced 
#'      by get_ranks_zscores()), plus a column of group assignments for that 
#'      trait (as produced by trait_groups()).
#' 
#' @param trait (string)
#'      The name of the trait of interest. Must be the column header for the 
#'      trait predictions to be plotted.
#' 
#' @param assignment_col (string)
#'      The column header for the column of group assignments to be plotted.
#' 
#' @param jitter (numeric)
#'      (Default c(80,40)) A vector of jitter magnitudes for non-assigned and 
#'      assigned points, respectively. This may need tweaking for aesthetic
#'      plotting.
#' 
#' @param gen (int)
#'      The RATTACA generation being plotted.
#' 
#' @param random_seed (int)
#'      (Default 1) An integer value with which to set the seed for random
#'      jittering. This allows consistent placement of jittered points and may 
#'      need tweaking for aesthetic plotting.
#' 
#' @param trait_name (string)
#'      (Default NULL) The desired trait name to use in the figure title and 
#'      file name. Use this option to simplify naming when a trait has a long 
#'      or unintuitive variable name. 
#' 
#' @param outdir (string)
#'      (Defualt NULL) The directory path in which to save the figure. If NULL, 
#'      the figure will plot to the R console.
#' 
#' @return Plots to the R console if outdir=NULL. Saves a png file if outdir is
#'          not NULL.
#
plot_assignments <- function(df, 
                            trait, 
                            assignment_col, 
                            jitter=c(80,40), # jitter for (background, assigned) points; set NULL to turn off jitter
                            gen, 
                            random_seed=1,
                            trait_name=NULL,
                            outdir=NULL){

    if (is.null(trait_name)) {
        trait_name <- trait
    }

    if (!is.null(outdir)) {
        if (!is.null(jitter)) {
            out_stem <- file.path(outdir, paste0('rattaca_gen', gen, '_', trait_name, '_assignment_jittered.png'))
            png(out_stem, width=7, height=5, units='in', res=300)
        } else {
            out_stem <- file.path(outdir, paste0('rattaca_gen', gen, '_', trait_name, '_assignment_.png'))
            png(out_stem, width=7, height=5, units='in', res=300)
        }
        png(out_stem, width=7, height=5, units='in', res=300)
    }
    
    trait_df <- df[df[[assignment_col]]==1 | df[[assignment_col]]==T | df[[assignment_col]]=='True',]
    trait_rank <- paste0(trait,'_rank')
    group_col <- paste0(trait,'_group')
    
    set.seed(random_seed)
    
    if (!is.null(jitter)) {
        x_all <- jitter(df[[trait_rank]], jitter)
        y_all <- jitter(df[[trait]], jitter[1])
        x_trait_high <- jitter(trait_df[trait_df[[group_col]]=='high', trait_rank], jitter[2])
        y_trait_high <- jitter(trait_df[trait_df[[group_col]]=='high', trait], jitter[2])
        x_trait_low <- jitter(trait_df[trait_df[[group_col]]=='low', trait_rank], jitter[2])
        y_trait_low <- jitter(trait_df[trait_df[[group_col]]=='low', trait], jitter[2])        
    } else {
        x_all <- df[[trait_rank]]
        y_all <- df[[trait]]
        x_trait_high <- trait_df[trait_df[[group_col]]=='high', trait_rank]
        y_trait_high <- trait_df[trait_df[[group_col]]=='high', trait]
        x_trait_low <- trait_df[trait_df[[group_col]]=='low', trait_rank]
        y_trait_low <- trait_df[trait_df[[group_col]]=='low', trait]
    }

    par(mar=c(6,6,4,2))
    plot(x_all, y_all, pch=16, cex=1.3, col=alpha(1,0.15), xlab='', ylab='')
    points(x_trait_high, y_trait_high, cex=1.6, pch=16, col=inferno(1,1,0.7))
    points(x_trait_high, y_trait_high, cex=1.6, lwd=1.4)
    points(x_trait_low, y_trait_low, cex=1.6, pch=16, col=inferno(1,1,0.25))
    points(x_trait_low, y_trait_low, cex=1.6, lwd=1.4)
    title(line=2.6, xlab=bquote(bold(.(trait)))); title(line=4, xlab=expression(bold('prediction rank')))
    title(line=4, ylab=bquote(bold(.(trait)))); title(line=2.6, ,ylab=expression(bold('prediction')))
    title(line=2.3, main=paste('RATTACA gen', gen))
    title(line=0.9, main=paste(trait, 'assignments'))
    title(line=2.6, xlab=bquote(bold(.(trait)))); title(line=4, xlab=expression(bold('prediction rank')))
    title(line=4, ylab=bquote(bold(.(trait)))); title(line=2.6, ,ylab=expression(bold('prediction')))
    title(line=2.3, main=paste('RATTACA gen', gen))
    title(line=0.9, main=paste(trait_name, 'assignments'))

    if (!is.null(outdir)) {
        dev.off()
    }
} 



#' Produce an assignment density plot: trait z-score density.
#' 
#' @description
#' Reads a dataframe of predictions, z-scores, and group assignments to plot the
#' density of prediction z-scores, with a rug of individual z-scores colored by 
#' assigned high/low group.
#'
#' @export
#'
#' @param df (dataframe)
#'      A dataframe of trait predictions for a trait of interest (as produced 
#'      by get_ranks_zscores()), plus a column of group assignments for that 
#'      trait (as produced by trait_groups()).
#' 
#' @param trait (string)
#'      The name of the trait of interest. Must be the column header for the 
#'      trait predictions to be plotted.
#' 
#' @param assignment_col (string)
#'      The column header for the column of group assignments to be plotted.
#' 
#' @param gen (int)
#'      The RATTACA generation being plotted.
#' 
#' @param trait_name (string)
#'      (Default NULL) The desired trait name to use in the figure title and 
#'      file name. Use this option to simplify naming when a trait has a long 
#'      or unintuitive variable name. 
#' 
#' @param outdir (string)
#'      (Defualt NULL) The directory path in which to save the figure. If NULL, 
#'      the figure will plot to the R console.
#' 
#' @return Plots to the R console if outdir=NULL. Saves a png file if outdir is
#'          not NULL.
#
density_assignments <- function(
    df, 
    trait, 
    assignment_col, 
    gen = NULL, 
    trait_name=NULL,
    outdir=NULL){

    if (is.null(trait_name)) {
        trait_name <- trait
    }
    if (!is.null(outdir)) {
        out_stem <- file.path(outdir, paste0('rattaca_gen', gen, '_', trait_name, '_assignment_density.png'))
        png(out_stem, width=7, height=5, units='in', res=300)
    }

    trait_rank <- paste0(trait,'_rank')
    group_col <- paste0(trait,'_group')
    zscore_col <- paste0(trait,'_zscore')

    trait_df <- df[df[[assignment_col]]==1 | df[[assignment_col]]==T | df[[assignment_col]]=='True',]
    high_df <- trait_df[trait_df[[group_col]]=='high',]
    low_df <- trait_df[trait_df[[group_col]]=='low',]

    d <- density(df[[zscore_col]])
    y_at_zero <- approx(d$x, d$y, xout=0)$y  # interpolate y at x=0

    plot(d, xlab='', ylab='', lwd=2, main='')
    polygon(d, col=alpha(1,0.15)) 
    segments(x0=0, y0=0, x1=0, y1=y_at_zero)    
    rug(df[[zscore_col]], col=alpha(1,0.2))
    rug(trait_df[[zscore_col]], lwd=1.6)
    points(low_df[[zscore_col]], rep(0, nrow(low_df)), col=inferno(1,1,0.25),pch=16, cex=1.4)
    points(high_df[[zscore_col]], rep(0, nrow(high_df)), col=inferno(1,1,0.7),pch=16, cex=1.4)
    points(trait_df[[zscore_col]], rep(0, nrow(trait_df)), cex=1.4, lwd=1.6)
    title(line=2.6, xlab=bquote(bold(.(trait_name)))); title(line=4, xlab=bquote(bold('prediction Z-score')))
    title(line=2.6, ,ylab=bquote(bold('Density')))
    if (!is.null(gen)) {
        title(line=2.3, main=paste('RATTACA gen', gen))
        title(line=0.9, main=paste(trait_name, 'assignments'))
    } else {
        title(line=1.2, main=paste(trait_name, 'assignments'))
    }

    # close the png device if opened
    if (!is.null(outdir)) {
        dev.off()
    }

}

#' Print time-stamped messages to stdout.
#' 
#' @export
#'
#' @param str (character)
#'      A text string to print.

printout <- function(str) {
    cat(paste0('\n[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ']'), 
        str, '\n\n')
}
