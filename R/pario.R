# RATTACA parameter IO
#
# Author: Robert Vogel
# Date: 2023-08-21
#
# Contributors:
#

META_PREFIX <- "##"
META_KEY_VAL_DELIM <- "="
META_VALID_KEY_REGEX <- ".+"
META_VALID_VAL_REGEX <- ".+"
RECORD_DELIM <- "\n"
FIELD_DELIM <- ","
ENCODING <- "UTF-8"
HEADER_PREFIX <- "#"


#' Read the parameters written to file
#'
#' @export
#' 
#' @param filename (character)
#'
#' @return list with elements:
#' \itemize{
#'      \item{$meta}{list of key value pairs of meta data}
#'      \item{$pars}{labeled matrix of parameters}
#' }
read_pars <- function(filename)
{

    fconn <- file(description=filename,
                  open="rt")

    # load meta data

    # user R's implementation of regex capture groups for
    # extracting key value pairs.
    proto <- data.frame(key=character(), val=character())
    meta_str_pattern <- paste0("^",
                               META_PREFIX,
                               "(",
                               META_VALID_KEY_REGEX,
                               ")",
                               META_KEY_VAL_DELIM,
                               "(",
                               META_VALID_VAL_REGEX,
                               ")",
                               "$")

    meta <- list()

    line <- utils::strcapture(meta_str_pattern,
                              readLines(con=fconn, n=1),
                              proto, perl=TRUE)
    while(!any(is.na(line)))
    {
        meta[[line$key]] <- type.convert(line$val, as.is=TRUE)

        raw_line <- readLines(con=fconn, n=1)

        line <- utils::strcapture(meta_str_pattern,
                                  raw_line,
                                  proto, perl=TRUE)
    }
    

    # check whether header exists and load
    proto <- data.frame(index=character(),
                        u=character(),
                        u_se=character())
    header = utils::strcapture(paste0("^",
                                   HEADER_PREFIX,
                                   "(.+),(.+),(.+)$"),
                            raw_line,
                            proto)

    # parameter values
    data <- utils::read.csv(fconn, sep=",", header=FALSE)
    close(fconn)

    rownames(data) <- data$V1
    data$V1 <- NULL
    colnames(data) <- c(header$u, header$u_se)
    data <- as.matrix(data) 

    return(list("meta"=meta,
                "header"=header,
                "data"=data))
}


# TODO: how to close file connections in R upon an error
# TODO: add misc. key val pair inputs using ellipses to input args

#' Write inferred parameters to file
#'
#' @export
#'
#' @param trait (character) string specifying trait fit
#' @param trait_file (character) path to trait file used in fit
#' @param genotype_prefix (character) path to plink
#'      genotype files used for fitting
#' @param pars (list) output from rattaca::fit
#' @param filename (character) path of file to write
#'
#
write_pars <- function(trait, trait_file, genotype_prefix,
                      pars, filename)
{
    sess_info <- utils::sessionInfo()

    fconn <- file(description=filename,
                  open="wt",
                  encoding=ENCODING)

    writeLines(c(paste(paste0(META_PREFIX,"date"),
                       Sys.time(),
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "user"),
                       Sys.getenv("USER"),
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "trait"),
                       trait,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX,"trait_file"),
                       trait_file,
                       sep=META_KEY_VAL_DELIM),
                 paste(paste0(META_PREFIX, "plink_genotypes_prefix"),
                       genotype_prefix,
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
    
    writeLines(c(paste(paste0(META_PREFIX, "variance_u"),
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

    for (rname in rownames(pars$u))
    {
        writeLines(paste(rname,
                         pars$u[rname],
                         pars$u.SE[rname],
                         sep=FIELD_DELIM),
                   con=fconn)
    }
    close(fconn)

}
