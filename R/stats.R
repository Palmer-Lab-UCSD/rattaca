# Statistical helper functions for RATTACA
#
# By: Robert Vogel
# Date: 2023-08-20
#
# Contributors
#

# Create a rattaca power analysis object
new_power_class <- function(pval, group_n, n_clones, trait_sim) {
    structure(pval, class = 'power_analysis', 
        'n_clones' = n_clones, 
        'group_n' = group_n, 
        'trait_sim' = trait_sim)
}

# function to check if a vector of IDs includes clones
# assumes clones are identified with an underscore and integer, e.g. id_1, id_2
count_clones <- function(ids) {
    
    id_counts <- table(sapply(strsplit(ids, '_'), `[[`, 1))

    if (any(duplicated(ids))) {
        cat('Error: some IDs are duplicated: \n')
        print(ids[duplicated(ids)])
        return()
    }
    
    if (length(unique(id_counts)) == 1) {
        n_clones <- id_counts[1]; names(n_clones) <- NULL
        if (n_clones > 1) {        
            cat('Cloned dataset:',n_clones,'clones per ID \n')
        }
        return(n_clones)
    } else {
        cat('Error: some IDs have different numbers of clones: \n')
        print(id_counts)
        return()
    }
    
}


#' Power analysis given genotypes and model simulator
#'
#' @export
#'
#' @param geno_low ((n_low samples, q markers) array | matrix) 
#'      genotypes of samples whose mean predicted value is
#'      less than geno_high
#' @param geno_high ((n_high sample, q marker) array | matrix)
#'      genotypes of samples whose mean predicted value is 
#'      larger than geno_low
#' @param sim (function) a function from simulator closure.
#' @param significance_level (double)
#' @param m_power_reps (int)
#'
#' @return (power_analysis) An object of class power_analysis: an estimate of
#' statistical power (class double), with attributes for sample sizes and 
#' simulated data
#
power_analysis <- function(geno_low, geno_high, sim,
                  significance_level=0.05,
                  m_power_reps=100,
                  trait=NULL)
{

    # get samples sizes of the two groups
    n_low <- nrow(geno_low)
    n_high <- nrow(geno_high)
    # group_n <- c(n_high, n_low); names(group_n) <- c('high','low')
    
    # check if power is being conducted on a cloned dataset
    n_clones_low <- count_clones(rownames(geno_low))
    n_clones_high <- count_clones(rownames(geno_high))

    if (n_clones_low != n_clones_high) {
        cat('Error: Different numbers of clones between high and low groups')
        return()
    } 
    n_clones <- n_clones_low

    
    # list to store simulated trait values
    sim_pheno <- list()

    # counter for experimental successes (differences between high vs low)
    num_reject_null <- 0
    
    for (i in seq(m_power_reps)) {
                
        tmp_high <- sim(geno_high)
        tmp_low <- sim(geno_low)

        # sim_pheno <- list('high' = tmp_high, 'low' = tmp_low)
        sim_pheno_df <- data.frame(
            trait = ifelse(is.null(trait),NA,trait),
            rep = i,
            group = c(rep('high',length(tmp_high)), rep('low',length(tmp_low))),
            rfid = c(names(tmp_high),names(tmp_low)),
            sim_pred = c(tmp_high, tmp_low)
        )
        # identify clones in the simulation df
        sim_pheno_df$clone <- ifelse(grepl("_", sim_pheno_df$rfid), 
                                     sub(".*_", "", sim_pheno_df$rfid), 1)
        sim_pheno_df$rfid  <- sub("_.*", "", sim_pheno_df$rfid)
        col_order <- c('trait','rep','group','rfid','clone','sim_pred')
        sim_pheno_df <- sim_pheno_df[,col_order]
        rownames(sim_pheno_df) <- NULL

        sim_pheno[[i]] <- sim_pheno_df
        
        result <- stats::t.test(tmp_high, tmp_low,
                             alternative="greater")

        if (result$p.value < significance_level)
            num_reject_null <- num_reject_null + 1
    
    } # end of m_power_reps loop

    # store simulated data and summaries
    sim_pheno_df <- do.call(rbind, sim_pheno); row.names(sim_pheno_df) <- NULL
    group_means <- aggregate(sim_pred ~ rep + group, data = sim_pheno_df, FUN = mean)
    colnames(group_means)[3] <- 'mean'

    power <- num_reject_null / m_power_reps

    out <- new_power_class(power, 
        trait = trait, n_high = n_high, n_low = n_low, n_clones = n_clones, trait_sim = sim_pheno_df)
   
    return(out)
}


#' Compute the coefficient of determination (\eqn{R^2})
#'
#' @export
#'
#' @param true_vals ((n,) array or vector)
#' @param predictions ((n,) array or vector)
#'
#' @return (double)
#
compute_r_sq <- function(true_vals, predictions)
{
    if (length(true_vals) != length(predictions))
        stop("Number of true_vals and predictions must be equal")

    ssr <- sum((true_vals - predictions)^2)
    stot <- sum((true_vals - mean(true_vals))^2)

    return(1 - ssr /stot)
}


#' Mean imputation of NA values in a vector
#'
#' @description
#' Consider an \eqn{n} length vector \eqn{v}. If not all elements
#' of \eqn{v} are not NA, but a subset are.  Replace each NA
#' element by the mean of the elements of \eqn{v} that are not NA.
#'
#' @param v ((n,) vector)
#'
#' @return ((n,) vector) in which NA values are replaced by
#'      imputed values.
mean_imputation <- function(v)
{
    num_na_vals <- sum(is.na(v))

    if (num_na_vals == 0)
        return(v)
    else if (num_na_vals < length(v)) {

        v[is.na(v)] <- mean(v, na.rm=TRUE)
        return(v)
    }        

    # all values are NA, stop evaluation
    stop("All values are NA")
}


#' Impute mean genotypes of NA values in genotype matrix
#'
#' @description
#' Consider the genotype \eqn{g_{ij}} of sample \eqn{j} and
#' genetic marker (e.g. SNP) \eqn{i}.  If \eqn{g_{ij} = }NA, then 
#' replace \eqn{g_{ij}} by the mean of \eqn{\bar{g}_i} of marker
#' \eqn{i} over all samples \eqn{k} such that \eqn{g{ik}\neq}NA.
#'
#' @export
#'
#' @param genotypes ((q markers, n samples) array | matrix)
#'
#' @return ((n samples, q markers) array | matrix) data with
#'      imputed values.
#
impute <- function(genotypes)
{
    if (length(dim(genotypes)) != 2)
        stop("Genotypes must by an (n sample, q marker) matrix")

    return(apply(genotypes, 1, mean_imputation))
}


#' Perform a k-fold cross-validation on rrBLUP predictions for one trait
#'
#' @export
#'
#' @param data (list)
#'      An aligned genotype/phenotype dataset, as produced by rattaca::align:
#'      A list with elements $trait (character string naming the trait to
#'      analyze), $pheno (a named vector of phenotype observations), and $geno
#'      (a genotype matrix aligned with $pheno)
#' 
#' @param num_folds (int)
#'      The desired number of folds (k) used for k-fold cross-validation
#' 
#' @param out_dir (character)
#'      The path to the directory in which to save output files summarizing
#'      cross-validation performance
#'
#' @return A list of (1) the trait name, (2) the lists of RFIDs including in the
#'      train and test sets of each fold, (3) the list of k model fits onto each
#'      training sample (as produced by fit()), and (3) the list of k model
#'      validations of each test sample (as produced by validate_test_preds())
#
kfold_cv <- function(data, num_folds, out_dir)
{
            
    # set up k-fold cross validation: train/test split on RATS
    # by randomly assigning individuals to 1 of k groups
    trait_df <- data.frame(trait_val = data$pheno, kfold = NA)    
    colnames(trait_df)[1] <- data$trait
    obs_per_fold <- ceiling(length(data$pheno) / num_folds)
    remaining_indices <- 1:length(data$pheno)

    for (fold in 1:num_folds) {

        if (fold < num_folds) {

            fold_size <- obs_per_fold

      } else {

            fold_size <- length(remaining_indices)
      }

        fold_indices <- sample(remaining_indices, size = fold_size, replace = FALSE)
        trait_df$kfold[fold_indices] <- fold
        remaining_indices <- setdiff(remaining_indices, fold_indices)
    }

    # empty lists to store model parameters and performance
    train_fits <- list()
    test_fits <- list()
    splits <- list()

    # k-fold cross validation
    for (k in 1:num_folds){
        
        cat(data$trait, 'CV fold', paste0(k,':'), format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')

        # set up phenotype train/test sets
        test_rfids <- rownames(trait_df[trait_df$kfold == k,])
        train_rfids <- rownames(trait_df[trait_df$kfold != k,])

        splits[[k]] <- list(train = train_rfids, test = test_rfids)

        # split phenotype data
        pheno_train <- trait_df[train_rfids, data$trait]
        names(pheno_train) <- train_rfids
        pheno_test <- trait_df[test_rfids, data$trait]
        names(pheno_test) <- test_rfids

        # split genotype data
        geno_train <- data$geno[train_rfids,]    
        geno_test <- data$geno[test_rfids,]

        # run rrBLUP on the training set
        train_fits[[k]] <- fit(pheno_train, geno_train)

        # re-fit the trained model to the test data
        test_fits[[k]] <- validate_test_preds(pheno_test, geno_test, train_fits[[k]])
                
    } # end of k-fold loop
    
    model_results <- list(trait = data$trait, splits = splits, train = train_fits, test = test_fits)

    # write cross-validation results to files
    save_cv_results(model_results, out_dir)

    # test_results <- model_results$test
    
    # obs <- c()
    # pred <- c()
    # r_sq <- c()
    # r <- c()
    # rho <- c()
    # fold <- c()

    # for (k in 1:length(test_results)) {
    #     out <- model_results$test[[k]]
    #     obs <- c(obs, out$obs)
    #     pred <- c(pred, out$pred)
    #     fold <- c(fold, rep(k, length(out$obs)))
    #     r_sq <- c(r_sq, rep(out$r_sq, length(out$obs)))
    #     r <- c(r, rep(out$pearson_corr, length(out$obs)))
    #     rho <- c(rho, rep(out$spearman_corr, length(out$obs)))
    # }

    # cv_df <- data.frame(
    #     rfid = names(obs),
    #     trait = rep(trait, length(obs)),
    #     fold = fold,
    #     obs = obs,
    #     pred = pred,
    #     r_sq = r_sq,
    #     r = r,
    #     rho = rho)
    # outfile <- paste0(data$trait,'_',num_folds,'fold_cv.csv')
    # outfile <- file.path(out_dir, outfile)
    # write.csv(cv_df, outfile, row.names=F, quote=F, na='')
    # cat('Cross-validation dataset written to', outfile, '\n')

    # summary_df <- data.frame(
    #     trait = trait,
    #     fold = unique(fold),
    #     r_sq = unique(r_sq),
    #     r = unique(r),
    #     rho = unique(rho))
    # outfile <- paste0(data$trait,'_',num_folds,'fold_cv_summary.csv')
    # outfile <- file.path(out_dir, outfile)
    # write.csv(summary_df, outfile, row.names=F, quote=F, na='')
    # cat('Cross-validation summary written to', outfile, '\n\n')

    return(model_results)
} 


#' Conduct a PCA on train- and test-set genotypes using Plink
#'
#' @export
#'
#' @param train_genotypes (character)
#'      The file path/prefix for the training Plink dataset.
#' 
#' @param test_genotypes (character)
#'      The file path/prefix for the test Plink dataset.
#'
#' @param trait (character)
#'      The name of the trait to be analyzed
#' 
#' @return A list of (1) the trait name, (2) the path to the Plink file
#'      of PCA results
#
pca_plink_genotypes <- function(train_genotypes, test_genotypes, trait)
{

    train_file <- train_genotypes
    test_file <- test_genotypes

    # merge train/test genotype data
    merge_args <- c('-bfile', train_file, 
                    '--bmerge', test_file,
                    '--make-bed',
                    '--out', file.path(pca_dir, paste0(trait, '_merged'))
                    )
    
    # print the plink call to the user
    cat('Plink merge call:', plink2, paste(merge_args, collapse=' '), '\n')

    # run plink
    system2(plink1, merge_args)

    # conduct PCA
    pca_args <- c('-bfile', file.path(pca_dir, paste0(trait, '_merged')),
                  '--pca',
                  '--out', file.path(pca_dir, paste0(trait, '_pca'))
                  )
    
    # run plink
    system2(plink2, pca_args)

    # print the plink call to the user
    cat('Plink PCA call:', plink2, paste(pca_args, collapse=' '), '\n')
    
    return(list(trait = trait, pca_file = file.path(pca_dir, paste0(trait, '_pca'))))

}


#' Standard-scale a dataset into Z-scores
#'
#' @export
#'
#' @param data (numeric)
#'      Any numeric vector
#' 
#' @return A vector of Z-scores 
#
zscore <- function(data){

    z <- (data - mean(data)) / sd(data)
    return(z)

}

#' Estimate per-group means for trait simulations produced during power analysis
#'
#' @description
#' Estimate per-group means for trait simulations produced during power analysis
#'
#' @export
#'
#' @param pwr_list (list)
#'      A list of power_analysis objects 
#'
#' @return A dataframe of simulated group means for each group x replicate 
#'      provided in the input
#

get_group_means <- function(pwr_list) {
    
    pwr_means <- lapply(seq_along(pwr_list), function(i) { 
        x <- pwr_list[[i]]                                 
        trait_sim <- attr(x, 'trait_sim')
        n_clones <- attr(x, 'n_clones')
        n_high <- attr(x, 'n_high')
        n_low <- attr(x, 'n_low')
        means <- aggregate(sim_pred ~ group, data = trait_sim, FUN = mean)
        se    <- aggregate(sim_pred ~ group, data = trait_sim,
                           FUN = function(d) sd(d) / sqrt(length(d)))
        var   <- aggregate(sim_pred ~ group, data = trait_sim, FUN = var)
        sd    <- aggregate(sim_pred ~ group, data = trait_sim, FUN = sd)
        group_means <- data.frame(
            n_clones = n_clones,
            group = means$group,
            mean  = means$sim_pred,
            se    = se$sim_pred,
            var = var$sim_pred,
            sd = sd$sim_pred
        )
        group_means$group_n <- sapply(group_means$group, function(x) {
            ifelse(x=='high', n_high, n_low)
        })
        group_means$rep <- i                        
        col_order <- c('rep','group','group_n','n_clones','mean','se','var','sd')
        group_means <- group_means[, col_order]
        return(group_means)
    })
    pwr_means <- do.call(rbind, pwr_means); rownames(pwr_means) <- NULL
    return(pwr_means)
}


#' Conduct a RATTACA power analysis
#' 
#' @description
#' This is a wrapper for the basic power_analysis() function. Can conduct 
#' multiple experimental replicates of a power analysis (each with a given 
#' number of bootstrap replicates) and enables writing results, including 
#' simulated data, to file.
#'
#' @export

#' @param geno_low ((n_low samples, q markers) array | matrix) 
#'      Genotypes of samples whose mean predicted value is
#'      less than geno_high
#' 
#' @param geno_high ((n_high sample, q marker) array | matrix)
#'      Genotypes of samples whose mean predicted value is 
#'      larger than geno_low
#' 
#' @param sim (function) 
#'      A function from simulator closure, as output by gen_trait_sim_closure().
#' 
#' @param significance_level (double)
#'      (default 0.05) The desired alpha level at which to test differences
#'      between groups.
#' 
#' @param tests_per_rep (int)
#'      The number of bootstrap replicates to be conducted by each
#'      experimental replicate power analysis.
#' 
#' @param reps (int)
#'      The number of experimental replicates to conduct
#' 
#' @param trait (character)
#'      Name of the trait being analyzed
#' 
#' @param outdir (character)
#'      (default NULL) The directory path in which to save results, if provided.
#' 
#' @return (dataframe) A dataframe of metadata and statistical power for all
#'      experimental replicates of the analysis.
#
rattaca_power <- function(
    geno_low, 
    geno_high, 
    sim, 
    alpha=0.05, 
    tests_per_rep, 
    reps, 
    trait, 
    outdir=NULL) 
{

    pwr_list <- replicate(reps,
                          
        power_analysis2(
            geno_low = geno_low, 
            geno_high = geno_high, 
            sim = sim, 
            sig = alpha, 
            m_power_reps = tests_per_rep,
            trait = trait
        ), 
        simplify=FALSE
    )

    summary <- list()

    pwr_df <- lapply(pwr_list, function(pwr) {
        data.frame(
            trait = trait,
            n_low = attr(pwr, 'n_low'),
            n_high = attr(pwr, 'n_high'),
            n_clones = attr(pwr, 'n_clones'),
            alpha = alpha,
            tests_per_rep = tests_per_rep,
            power = pwr
        )
    })
    summary <- c(summary, pwr_df)
    summary <- do.call(rbind, summary)
    summary$rep <- 1:nrow(summary)
    col_order <- c('trait','rep','tests_per_rep','alpha','n_clones','n_low','n_high','power')
    summary <- summary[,col_order]

    if(!is.null(outdir) && dir.exists(outdir)) {
        
        # save power summary to file
        outfile <- paste0(trait,'_power_n',tests_per_rep,'_tests_n',reps,'_reps_summary.csv')
        outfile <- file.path(outdir, outfile)
        write.csv(summary, outfile, row.names=F, quote=F, na='')
        cat(trait, 'power summary written to', outfile, '\n')

        # save group means to file
        group_means <- get_group_means(pwr_list)
        outfile <- paste0(trait,'_power_n',tests_per_rep,'_tests_n',reps,'_reps_group_means.csv')
        outfile <- file.path(outdir, outfile)
        write.csv(group_means, outfile, row.names=F, quote=F, na='')
        cat(trait, 'power simulated trait means written to', outfile, '\n')

        # save trait simulations to file
        lapply(seq_along(pwr_list), function(i) { 

            rep_str <- paste0('rep_',paste0(rep(0,nchar(reps)-nchar(i)),collapse=''), i)
            x <- pwr_list[[i]]         
            trait_sim <- attr(x, 'trait_sim')
            sim_dir <- file.path(outdir, 'simulations')
            dir.create(sim_dir, recursive=TRUE, showWarnings=FALSE)
            outfile <- paste0(trait,'_power_n',tests_per_rep,'_tests_',rep_str,'.csv')
            outfile <- file.path(sim_dir, outfile)
            write.csv(trait_sim, outfile, row.names=F, quote=F, na='')
            cat(trait, 'power simulations written to', outfile, '\n')

        })
    }
    
    return(summary)

}
