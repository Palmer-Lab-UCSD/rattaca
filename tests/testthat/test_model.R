# test model.R functions
#
# Author: Robert Vogel
# Date: 2023-08-24
#
# Contributors:
#
# sim_data <- function(n_samples, q_markers)
# {
#     q_markers <- 100
#     n_samples <- 50
#     genotypes <- sample(c(-1, 0, 1), q_markers * n_samples,
#                         replace=TRUE)
#     genotypes <- matrix(genotypes,
#                         ncol=q_markers,
#                         nrow=n_samples)
#     
#     marker_random_effects <- stats::rnorm(q_markers, 0, 1)
# }

# test_that("predict_lmm: no imputed genotypes",
#           {
#               marker_random_effects 
#           }
# )
