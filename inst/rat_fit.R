# Run rrBLUP on specified data and save parameters to file
#
# By: Robert Vogel
# Date: 2023-08-10
#

# TODO: Existing files will be overwritten!
# TODO: Need to remove the addition of library path

library(rattaca)


VERSION = "1.0.2"
LOG_DELIMITER = " : "


# Predict trait values under LMM model and rrBLUP pkg.

main <- function(args)
{
    
    if (!dir.exists(args[["--out_dir"]]))
    {
        dir.create(args[["--out_dir"]], 
                   mode = "0750")
    }

    output_filename <- paste0(file.path(args[["--out_dir"]],
                                        args[["--trait"]]),
                             ".bpar")
    # load data
    print(paste(Sys.time(),
                "Loading data",
                sep=LOG_DELIMITER))

    genotypes <- rattaca::load_and_prepare_plink_data(
                         args[["--plink_genotypes_prefix"]])

    trait <- rattaca::load_and_prepare_trait_data(
                               args[["--trait_file"]],
                               args[["--trait_rat_id_colname"]],
                               args[["--trait"]])

    # align data tables
    # TODO: what is the exact data structure for genotypes and traits
    #       an array, data frame, matrix?
    rat_ids <- intersect(rownames(trait), rownames(genotypes))
    genotypes <- genotypes[rat_ids,]
    trait <- trait[rat_ids,]
    

    print(paste(Sys.time(),
                "Running rrBLUP to fit LMM model.",
                sep=LOG_DELIMITER))

    out_pars <- rattaca::fit(trait[, args[["--trait"]]],
                             genotypes)

    print(paste(Sys.time(), 
                "Done fitting, writing results to disk",
                sep=LOG_DELIMITER))


    out <- rattaca::write_pars(args[["--trait"]],
                            args[["--trait_file"]],
                            args[["--plink_genotypes_prefix"]],
                            out_pars,
                            output_filename)
}


# Parse input arguments using the options in the function
#
# Args:
#   args (vector) of character strings
#
# Returns:
#   (list)
parse_args <- function(args)
{

    parser <- rattaca::argument_parser(
            "--plink_genotypes_prefix" = rattaca::argument(
                help=paste("path and filename prefix of plink",
                            "bed, bim, and fam files.",
                            sep=" ")),
            "--trait_file" = rattaca::argument(
                help=paste("path and filename of csv table",
                           "with measured traits per animal.",
                           sep=" ")),
            "--trait" = rattaca::argument(
                help="name of trait used in trait table."),
            "--trait_rat_id_colname" = rattaca::argument(
                default_val="rfid",
                help="Column name of rat id's in trait file.",
                required=FALSE),
            "--genotype_rat_id_colname" = rattaca::argument(
                default_val="id",
                help="Column name of rat id's in fam file.",
                required=FALSE),
            "--out_dir" = rattaca::argument(
                help="Path to directory to write results"),
            description=paste(
                    "Train a Linear Mixed Model (LMM) using",
                    "the rrBLUP package. The model learns the",
                    "random effect sizes of each SNP, prediction",
                    "error, etc.",
                     sep="\n")
            )

    return(parser(args))
}



# run main program if R script, otherwise just load
# functions
if (!interactive())
{
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    
    if (is.null(args))
    {
        quit()
    }

    for (arg_item in args)
    {
        if (is.null(arg_item))
        {
            stop("Please include all required arguments, see --help.")
        }
    }

    out <- main(args)
}
