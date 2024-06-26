% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{gen_trait_sim_closure}
\alias{gen_trait_sim_closure}
\title{Generate an LMM simulation closure for a trait}
\usage{
gen_trait_sim_closure(intercept, u, sd_error, intercept_se = NULL, u_se = NULL)
}
\arguments{
\item{intercept}{(double)}

\item{u}{((q markers,) vector) 
vector of BLUPs, that is the expected value of the
random marker effect sizes}

\item{sd_error}{(double)
standard deviation of the model error.}

\item{u_se}{(NULL | (q_markers,) vector)
of standard error of the BLUPs}

\item{int_se}{(NULL | double)
standard deviation of the LMM random error term.}
}
\value{
(function)
     for generating samples from genotypes under the 
     rrBLUP LMM
}
\description{
Given the parameters of the LMM return a function to
generate independent samples from a given genotype matrix.
}
