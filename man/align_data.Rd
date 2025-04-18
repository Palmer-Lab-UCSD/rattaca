% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{align_data}
\alias{align_data}
\title{Align genotype and phenotype datasets for prediction.}
\usage{
align_data(genotypes, phenotypes, trait)
}
\arguments{
\item{genotypes}{(list)
A genotype dataset as produced by make_plink_dataset() or
sample_snps_from_plink_files(): a list with elements $geno_file
(the path/prefix to the Plink dataset) and $geno (a genotype matrix
of n samples x q variants)}

\item{phenotypes}{(numeric)
A named numeric vector of phenotype measurements}

\item{trait}{(character)
The name of the trait being analyzed}
}
\value{
A list of (1) the trait name, (2) all sample IDs shared between
     datasets, (3) the path/prefix for the Plink dataset used in
     alignment, (4) the aligned genotype matrix, and (5) the aligned
     phenotype data
}
\description{
Keeps only samples that are shared between datasets, maintaining the same 
order of sample IDs in both
}
