% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{impute}
\alias{impute}
\title{Impute mean genotypes of NA values in genotype matrix}
\usage{
impute(genotypes)
}
\arguments{
\item{genotypes}{((q markers, n samples) array | matrix)}
}
\value{
((n samples, q markers) array | matrix) data with
     imputed values.
}
\description{
Consider the genotype \eqn{g_{ij}} of sample \eqn{j} and
genetic marker (e.g. SNP) \eqn{i}.  If \eqn{g_{ij} = }NA, then 
replace \eqn{g_{ij}} by the mean of \eqn{\bar{g}_i} of marker
\eqn{i} over all samples \eqn{k} such that \eqn{g{ik}\neq}NA.
}
