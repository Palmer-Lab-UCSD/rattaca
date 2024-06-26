% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{argument_parser}
\alias{argument_parser}
\title{Parse command line arguments}
\usage{
argument_parser(..., description = NULL)
}
\arguments{
\item{...}{((character) key = (list) value) pairs where
keys are command line options (including -- prefix) in
quotes, and value is a list generated by rattaca::arg}

\item{description}{(character)}
}
\value{
(list) ((character) key = (object)  value) from
     option defaults of specified values from the command
     line
}
\description{
Parse command line arguments
}
