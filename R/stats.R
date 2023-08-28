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
                  m_power_reps=100)
{
    num_reject_null <- 0

    for (i in seq(m_power_reps))
    {
        tmp_high <- sim(geno_high)
        tmp_low <- sim(geno_low)

        out <- stats::wilcox.test(tmp_high, y=tmp_low,
                                  alternative=c("greater"),
                                  paired=FALSE, mu=0)

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

    return(1 - stats::var(true_vals - predictions) / stats::var(true_vals))
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
