

test_that("is_genotype return correct boolean vals",
          {
              genotypes <- matrix(seq(9), nrow=3, ncol=3)
              testthat::expect_false(rattaca::is_genotype(genotypes))
              testthat::expect_false(rattaca::is_genotype(genotypes,
                                                          imputed_vals=FALSE))


              genotypes <- matrix(sample(c(1,0,-1), 8, replace=TRUE),
                                  nrow=2,
                                  ncol=4)
              testthat::expect_true(rattaca::is_genotype(genotypes))
              testthat::expect_true(rattaca::is_genotype(genotypes,
                                                         imputed_vals=FALSE))

              genotypes[2,3] <- 0.3
              genotypes[1,2] <- -0.7
              testthat::expect_true(rattaca::is_genotype(genotypes))
              testthat::expect_false(rattaca::is_genotype(genotypes,
                                                          imputed_vals=FALSE))
          }
)


test_that("argument: default inputs",
          {
              # test that a list is returned and that default values are NULL
              default <- list(default_val=NULL, help=NULL,
                              required=TRUE, type="character")

              arg_out <- rattaca::argument()

              testthat::expect_true(is.list(arg_out))

              for (key in names(arg_out))
              {
                testthat::expect_identical(default[[key]], arg_out[[key]])
              }
          }
)

test_that("argument: valid specified input",
          {
              # double value input
              arg_in <- list(val=42, help="help string", type="double", required=TRUE)

              arg_out <- rattaca::argument(default_val = arg_in$val,
                                           help = arg_in$help,
                                           type = arg_in$type)
              for (key in names(arg_out))
              {
                  if (key == "required")
                      testthat::expect_true(arg_out[[key]])
                  else
                    testthat::expect_identical(arg_out[[key]], arg_in[[key]])
              }


              # character value input, no help specified
              arg_in <- list(val="cat", type="character", required=TRUE)

              arg_out <- rattaca::argument(default_val = arg_in$val,
                                           type = arg_in$type)
              for (key in names(arg_out))
              {
                  if (key == "required")
                      testthat::expect_true(arg_out[[key]])
                  else
                    testthat::expect_identical(arg_out[[key]], arg_in[[key]])
              }


              # character value input, no help specified, not required
              arg_in <- list(val="cat", type="character", required=FALSE)

              arg_out <- rattaca::argument(default_val = arg_in$val,
                                           type = arg_in$type,
                                           required=FALSE)
              for (key in names(arg_out))
              {
                  if (key == "required")
                      testthat::expect_false(arg_out[[key]])
                  else
                    testthat::expect_identical(arg_out[[key]], arg_in[[key]])
              }

          }
)


test_that("argument: invalid input",
          {
              testthat::expect_error(rattaca::argument(default_val=42.3, type="integer"))
          }
)


# TODO: How do we test stanard out when --help is specified?
# test_that("argument_parser: test specified default values",
#           {
#               # note that command line inputs will be character
#               args_target <- list("--x"=42, "--y"="a_string", "--z"=NULL)
#               args_in <- c("--x", "42", "--y", "a_string")
# 
#               parser <- rattaca::argument_parser(
#                         "--x" = rattaca::argument(default_val=args_target[["--x"]],
#                                        help_string="This is help",
#                                        type="double"),
#                         "--y" = rattaca::argument(default_val=args_target[["--y"]],
#                                        help_string="test"),
#                         "--z" = rattaca::argument(help_string="z test",
#                                        required=FALSE)
#                        )
#               
#               testthat::expect_true(is.function(parser))
# 
#               args_out <- parser(args_in)
# 
#               for (key in names(args_target))
#               {
#                   testthat::expect_identical(args_out[[key]], arg_in[[key]])
#               }
#           }
# )
