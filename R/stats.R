# Statistical helper functions for RATTACA
#
# By: Robert Vogel
# Date: 2023-08-20
#
# Contributors
#

# TODO
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
#' @return (double) statistical power
#
power_analysis <- function(geno_low, geno_high, sim,
                  significance_level=0.05,
                  m_power_reps=100) #stat_test="ttest")
{
    num_reject_null <- 0

    #if (test == "ttest")
    #    compute_stat <- stats::t.test
    # else if (test == "lmm")
    #     
    # else
    #     stop(paste("Provided stat_test is not one of"
    #                "the supported - ttest, lmm - tests"))

    for (i in seq(m_power_reps))
    {
        tmp_high <- sim(geno_high)
        tmp_low <- sim(geno_low)

        out <- stats::t.test(tmp_high, tmp_low,
                             alternative="greater")

        if (out$p.value < significance_level)
            num_reject_null <- num_reject_null + 1
    }

    return(num_reject_null / m_power_reps)
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
#' @return A list of (1) the trait name, (2) the list of k model fits onto each
#'      training sample (as produced by fit()), and (3) the list of k model
#'      validations of each test sample (as produced by validate_test_preds())
#
kfold_cv <- function(data, num_folds)
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
    test_preds <- list()
    test_out <- list()

    # k-fold cross validation
    for (k in 1:num_folds){
    
        # set up phenotype train/test sets
        test_rfids <- rownames(trait_df[trait_df$kfold == k,])
        train_rfids <- rownames(trait_df[trait_df$kfold != k,])

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
        test_fit <- validate_test_preds(pheno_test, geno_test, train_fits[[k]]$u, train_fits[[k]]$beta)
        
        test_out[[k]] <- list(obs = test_fit$obs, pred = test_fit$pred, r_sq = test_fit$r_sq, 
                              pearson_corr = test_fit$pearson_corr, spearman_corr = test_fit$spearman_corr)
        
    } # end of k-fold loop
    
    model_results <- list(trait = data$trait, train = train_fits, test = test_out)
    return(model_results)
    
} 
