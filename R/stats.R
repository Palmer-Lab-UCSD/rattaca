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


#' Perform power analyses on a set of hypothetical assignment groups.
#' 
#' @description
#' Genotypes are sampled into hypothetical assignment (high/low) groups 
#' according to desired parameters, then used to simulate trait predictions 
#' under the given model for power analyses under a combination of parameters.
#'
#' @export
#'
#' @param trait (character)
#'      The trait name
#' 
#' @param genotypes (matrix)
#'      A named genotype matrix for all predicted animals
#' 
#' @param predictions (numeric)
#'      The named vector of trait predictions that will be used to separate
#'      animals into respective high or low assignment groups.
#' 
#' @param fitted_mod (list)
#'      A fitted model object, as output by rattaca::fit()
#' 
#' @param group_size (numeric)
#'      (default c(0.05, 0.1, 0.25, 0.33)) A numeric vector reflecting the 
#'      sample size(s) desired for each assigned group to compare in a power 
#'      analysis. Can be input either as proportions of the total predicted 
#'      sample or sample counts. For example, 0.1 will compare the highest 10% 
#'      and lowest 10% of predictions, and 20 will compare the highest 20 and 
#'      lowest 20 predictions.
#' 
#' @param low_ids (character)
#'      (default NULL) A vector of RFIDs that have already been assigned to a
#'      'low' group based on trait predictions.
#' 
#' @param high_ids character)
#'      (default NULL) A vector of RFIDs that have already been assigned to a
#'      'high' group based on trait predictions.
#' 
#' @param alpha (numeric)
#'      (default 0.05) The desired alpha value(s) (significance cutoff(s)) to 
#'      use when determining statistical power to compare assigned groups.
#' 
#' @param tests (numeric)
#'      (default 100) The desired number of replicate simulations to conduct
#'      for each power analysis
#' 
#' @param outdir (character)
#'      (default NULL) The desired output directory in which to save results
#' 
#' @return A dataframe with one row per combination of desired group sample 
#'      size and alpha, and corresponding statistical power.
#
# function to do a power analysis on different theoretical assigned samples
test_power <- function(
    trait,
    genotypes,   # geno matrix of all predicted samples
    predictions, # named vector of trait predictions
    fitted_mod,  # a fitted model (as output by fit()) 
    group_size = c(0.01, 0.025, 0.05, 0.1), # fraction or sample count for assignment group size
    low_ids = NULL,  # vector of IDs already assigned to the low sample
    high_ids = NULL, # vector of IDs already assigned to the high sample
    alpha = 0.05,    # vector of desired alpha cutoffs to test \
    tests = 100,      # number of replicate tests per power analysis (the number of comparisons to make btwn samples)
    outdir = NULL)
{
    preds <- sort(predictions)
    total_analyses <- length(group_size) * length(alpha) * tests
    milestone_analysis <- ceiling(total_analyses/10)
    progress <- 0 
    n_analyses <- 0 

    # empty vectors to store results
    sample_pct <- c()
    sample_n <- c()
    sample_type <- c()
    power <- c() 

    # create an LMM simulation closure for the fitted model
    sim <- gen_trait_sim_closure(
        intercept = fitted_mod$beta,
        u = fitted_mod$u,
        intercept_se = fitted_mod$beta.SE,
        u_se = fitted_mod$u.SE,
        sd_error = sqrt(fitted_mod$Ve))

    # read in assigned IDs
    if (!is.null(low_ids) & !is.null(high_ids)) {
        
        n_low <- length(low_ids)
        n_high <- length(high_ids)
        n_assigned <- max(c(n_low, n_high))
        group_size <- c(group_size, n_assigned)
        names(group_size)[1:length(group_size)] <- 'hypothetical'
        names(group_size)[length(group_size)] <- 'assigned'
    }

    # process each desired group size
    for (i in 1:length(group_size)) {

        group_sample <- group_size[i]
        group_type <- names(group_size)[i]
        
        # if groups are input as a fraction, calculate sample sizes
        if (group_sample < 1) {
            pct <- group_sample
            group_n <- ceiling(length(preds)*group_sample)
        } else { 
            pct <- round(group_sample/length(preds),3)
            group_n <- group_sample
        }

        # produce high & low samples: the n most extreme high or low IDs
        low_rfids <- names(preds)[1:group_n]
        high_rfids <- names(preds)[(length(preds)-group_n):length(preds)]

        # if testing an assigned group, use the actual assigned IDs
        if (group_type == 'assigned') {
            low_rfids <- low_ids
            high_rfids <- high_ids
        }

        # subset genotypes to the assigned sample
        geno_low <- genotypes[low_rfids,]
        geno_high <- genotypes[high_rfids,]

        # power analyses on each combination of group size x alpha
        for (j in 1:length(alpha)) {

            sig <- alpha[j]
            n_analyses = n_analyses + 1
            progress <- progress + 10

            if (n_analyses %% milestone_analysis == 0) {
                cat(paste0('\t[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ']'),
                    'Power analysis', paste0(n_analyses,':'), 
                    'group size =', paste0(group_sample, ','), 
                    'sigma =', sig, 
                    paste0('(',progress,'% complete)'), '\n')
            }
            
            pwr <- power_analysis(
                geno_low = geno_low, 
                geno_high = geno_high, 
                sim = sim, 
                sig = sig, 
                m_power_reps = tests) 
            power <- c(power, pwr)
        }

        sample_pct <- c(sample_pct, rep(pct, length(alpha)))
        sample_n <- c(sample_n, rep(group_n, length(alpha)))
        sig_cutoff <- rep(alpha, length(group_size))
    
    } # end of group_size for loop

    if (is.null(names(group_size))) {
        out_df <- data.frame(
            pct_cutoff = sample_pct,
            group_n = sample_n,
            alpha = sig_cutoff,
            tests = tests,
            power = power)
        out_df <- out_df[with(out_df, order(group_n, alpha)), ]
    } else {
        out_df <- data.frame(
            sample_type = rep(names(group_size), each = length(alpha)),
            pct_cutoff = sample_pct,
            group_n = sample_n,
            alpha = sig_cutoff,
            tests = tests,
            power = power)
        out_df <- out_df[with(out_df, order(group_n, sample_type, alpha)), ]
    }


    if (!is.null(outdir)) {
        write.csv(out_df, file.path(outdir, paste0(trait, '_test_power_n', 
            tests,'_tests.csv')),
              row.names=F, quote=F, na='')
    }

    return(out_df)

}


#' Perform replicated power analyses on a set of hypothetical assignment groups.
#' 
#' @description
#' Conducts replicate runs of test_power() to calculate the mean and standard 
#' error of power analysis summary statistics. 
#'
#' @export
#'
#' @param trait (character)
#'      The trait name
#' 
#' @param genotypes (matrix)
#'      A named genotype matrix for all predicted animals
#' 
#' @param predictions (numeric)
#'      The named vector of trait predictions that will be used to separate
#'      animals into respective high or low assignment groups.
#' 
#' @param fitted_mod (list)
#'      A fitted model object, as output by rattaca::fit()
#' 
#' @param group_size (numeric)
#'      (default c(0.05, 0.1, 0.25, 0.33)) A numeric vector reflecting the 
#'      sample size(s) desired for each assigned group to compare in a power 
#'      analysis. Can be input either as proportions of the total predicted 
#'      sample or sample counts. For example, 0.1 will compare the highest 10% 
#'      and lowest 10% of predictions, and 20 will compare the highest 20 and 
#'      lowest 20 predictions.
#' 
#' @param low_ids (character)
#'      (default NULL) A vector of RFIDs that have already been assigned to a
#'      'low' group based on trait predictions.
#' 
#' @param high_ids character)
#'      (default NULL) A vector of RFIDs that have already been assigned to a
#'      'high' group based on trait predictions.
#' 
#' @param alpha (numeric)
#'      (default 0.05) The desired alpha value(s) (significance cutoff(s)) to 
#'      use when determining statistical power to compare assigned groups.
#' 
#' @param reps (numeric)
#'      (default 10) The desired number of replicated power analyses to conduct 
#' 
#' @param tests_per_rep (numeric)
#'      (default 100) The desired number of replicate simulations to conduct
#'      for each power analysis. The number of comparisons to make between 
#'      samples (per replicate).
#' 
#' @param outdir (character)
#'      (default NULL) The desired output directory in which to save results
#' 
#' @return A dataframe with one row per combination of desired group sample 
#'      size and alpha, and corresponding statistical power.
#
# function to do multiple reps of test_power
test_power_multi <- function(
    trait,
    genotypes,       # geno matrix of all predicted samples
    predictions,     # named vector of trait predictions
    fitted_mod,      # a fitted model (as output by fit()) 
    group_size = c(0.01, 0.025, 0.05, 0.1), # fraction or sample count for assignment group size
    low_ids = NULL,  # vector of IDs already assigned to the low sample
    high_ids = NULL, # vector of IDs already assigned to the high sample
    alpha = 0.05,    # vector of desired alpha cutoffs to test
    reps = 10,       # number of outer reps: number of repeated runs of test_power()
    tests_per_rep = 100,     # number of internal replicate tests for test_power() (the number of comparisons to make btwn samples)
    outdir = NULL)
{
    # conduct multiple runs of test_power()
    results <- replicate(reps, test_power(
        trait = trait,
        genotypes = genotypes,
        predictions = predictions,
        fitted_mod = fitted_mod,
        group_size = group_size,
        low_ids = low_ids,
        high_ids = high_ids,
        alpha = alpha,
        tests = tests_per_rep,
        outdir = NULL
    ), simplify = FALSE)

    # combine results into a dataframe
    results <- do.call(rbind, results)

    # calculate mean and SE of power from each analysis
    if ('sample_type' %in% colnames(results)) {
        mean_power <- aggregate(power ~ sample_type + pct_cutoff + group_n + alpha, data = results, FUN = mean)
        sd_power <- aggregate(power ~ sample_type + pct_cutoff + group_n + alpha, data = results, FUN = sd)
        n_power <- aggregate(power ~ sample_type + pct_cutoff + group_n + alpha, data = results, FUN = length)
        se_power <- sd_power$power / sqrt(n_power$power)
    } else {
        mean_power <- aggregate(power ~ pct_cutoff + group_n + alpha, data = results, FUN = mean)
        sd_power <- aggregate(power ~ pct_cutoff + group_n + alpha, data = results, FUN = sd)
        n_power <- aggregate(power ~ pct_cutoff + group_n + alpha, data = results, FUN = length)
        se_power <- sd_power$power / sqrt(n_power$power)
    }

    summary <- data.frame(
        pct_cutoff = mean_power$pct_cutoff,
        group_n = mean_power$group_n,
        alpha = mean_power$alpha,
        reps = reps,
        tests_per_rep = tests_per_rep,
        power_mean = round(mean_power$power,3),
        power_se = round(se_power,3))

    if ('sample_type' %in% colnames(results)) {
        n_groups <- nrow(summary)
        sample_types <- results$sample_type[1:n_groups]
        sumcols <- colnames(summary)
        summary$sample_type <- sample_types
        summary <- summary[,c('sample_type',sumcols)]
    }

    if (!is.null(outdir)) {
        
        results_file <- file.path(outdir, paste0(trait, '_test_power_multi_n',
            tests_per_rep, '_tests_n', reps, '_reps.csv'))
        summary_file <- file.path(outdir, paste0(trait, '_test_power_multi_n',
            tests_per_rep, '_tests_n', reps, '_reps_summary.csv'))
        write.csv(results, results_file, row.names=F, quote=F, na='')
        write.csv(summary, summary_file, row.names=F, quote=F, na='')
    }

    cat('Replicate results written to', results_file, '\n')
    cat('Summarized results written to', summary_file, '\n')

    return(list(results = results, summary = summary))
}
