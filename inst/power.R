# Run power analysis on fit data.
#
# Perfomrs analysis na prints figures to specified directory.
#
# Author: Robert Vogel
# Date: 2023-09-29
#
# Contributors:
#

library(rattaca)


VERSION <- "0.0.1"
M_POWER_REPS <- 100
R_REPS <- 25
N_GROUPS <- 7
GROUP_SIZES <- 2^seq(N_GROUPS)

sem <- function(x)
    return(sd(x) / sqrt(length(x)))

errorbar <- function(x, y, xerr=NULL, yerr=NULL,
                     new=FALSE,
                     xlim=NULL, ylim=NULL,
                     xlab=NULL, ylab=NULL,
                     plot_line=FALSE,
                     pch=19,
                     col="black",
                     frac=0.05,
                     main=NULL)
{
    if (is.null(xlim) && is.null(xerr))
    {
        delta <- (max(x) - min(x)) * frac
        xlim <- c(min(x) - delta, max(x) + delta)
    } else if (is.null(xlim) && !is.null(xerr)) {
        delta <- (max(x + xerr) - min(x - xerr)) * frac
        xlim <- c(min(x - xerr) - delta, max(x + xerr) + delta)
    }
    
    if (is.null(ylim) && is.null(yerr))
    {
        delta <- (max(y) - min(y)) * frac
        ylim <- c(min(y) - delta, max(y) + delta)
    } else if (is.null(ylim) && !is.null(yerr)) {
        delta <- (max(y + yerr) - min(y - yerr)) * frac
        ylim <- c(min(y - yerr) - delta, max(y + yerr) + delta)
    }
        
    
    # c(bottom, left, top, right)
    if (new)
    {
        par(mar=c(5,5,4,2))
        plot.new()
        plot.window(xlim=xlim, ylim=ylim)
        axis(1, cex.axis=1.5)
        axis(2, cex.axis=1.5)
    }

    if (!is.null(main))
        mtext(main, side=3, cex=2)
    
    if (!is.null(xlab))
        mtext(xlab, side=1, line=3, cex=2)
    
    if (!is.null(ylab))
        mtext(ylab, side=2, line=3, cex=2)
    
    points(x, y, col=col, pch=pch, cex=1.5)

    if (plot_line)
        lines(x, y, col=col, lwd=2.5)
    
    if (!is.null(xerr))
        arrows(x-xerr, y, x+xerr, y, angle=90, code=3,
               length=0.05, lwd=2.5, col=col)
    
    if(!is.null(yerr))
        arrows(x, y-yerr, x, y+yerr, angle=90, code=3,
               length=0.05, lwd=2.5, col=col)
}



main <- function(args)
{

    genotypes <- rattaca::load_and_prepare_plink_data(
                        args[["--plink_genotypes_prefix"]])
    pars <- rattaca::read_pars(args[["--par_file"]])

    sim <- rattaca::gen_trait_sim_closure(pars$meta$intercept,
                        pars$data[,1],
                        sqrt(pars$meta$variance_e),
                        intercept_se = pars$meta$intercept_se,
                        u_se = pars$data[,2])

    # get ascening order of animals by breeding value
    idx <- sort.list(rattaca::predict_lmm(genotypes,
                                          pars$data[,1],
                                          pars$meta$intercept))

    # This is general yet, the point is that the maximum
    # group size needs to be less than or equal to 1/2
    # the number of samples
    i <- N_GROUPS
    n_samples <- dim(genotypes)[1]
    max_group_size <- n_samples / 2
    while (i > 0 && GROUP_SIZES[i] > max_group_size)
        i <- i-1

    if (i == 0)
        stop("Insufficient num. samples to perform analysis.")

    group_sizes <- GROUP_SIZES[1:i]
    n_group <- length(group_sizes)

    power <- list("mean"=vector(mode="double", length=n_group),
                  "sem"=vector(mode="double", length=n_group))

    delta <- list("mean"=vector(mode="double", length=n_group),
                  "sem"=vector(mode="double", length=n_group))

    sig <- list("mean"=vector(mode="double", length=n_group),
                "sem"=vector(mode="double", length=n_group))

    tp <- vector(mode="double", length=R_REPS)
    td <- vector(mode="double", length=R_REPS)
    ts <- vector(mode="double", length=R_REPS)

    for (i in seq(n_group))
    {
        idx_low <- idx[1:group_sizes[i]]
        idx_high <- idx[(n_samples-group_sizes[i]+1):n_samples]

        for (r in seq(R_REPS))
        {
            tp[r] <- rattaca::power_analysis(genotypes[idx_low,],
                            genotypes[idx_high,],
                            sim,
                            m_power_reps=M_POWER_REPS)

            t_tmp <- stats::t.test(sim(genotypes[idx_high,]),
                                   sim(genotypes[idx_low,]),
                                   alternative="greater")
            names(t_tmp$statistic) <- NULL
            td[r] <- t_tmp$statistic * t_tmp$stderr
            ts[r] <- t_tmp$stderr

        }

        power[["mean"]][i] <- mean(tp)
        power[["sem"]][i] <- sem(tp)
        delta[["mean"]][i] <- mean(td)
        delta[["sem"]][i] <- sem(td)
        sig[["mean"]][i] <- mean(ts)
        sig[["sem"]][i] <- sem(ts)
    }
        


    # Generate plots
    # plot power results

    if (!dir.exists(args[["--out_dir"]]))
        dir.create(args[["--out_dir"]])        


    pdf(file.path(args[["--out_dir"]], "power.pdf"))
    errorbar(group_sizes, power[["mean"]],
             yerr=power[["sem"]],
             xlab="Group size",
             ylab="Power",
             plot_line=TRUE,
             new=TRUE,
             ylim= c(0,1),
             main=paste("Trait:", pars$meta$trait,
                        "Pop. size:",
                        n_samples,
                        "\nErrorbar Reps:",
                        R_REPS,
                        "Power Reps:", M_POWER_REPS))
    dev.off()


    coluer <- hcl.colors(n_group, palette="Earth")

    pdf(file.path(args[["--out_dir"]], "mean_sd.pdf"))
    errorbar(sig[["mean"]][1], delta[["mean"]][1],
             xerr=sig[["sem"]][1],
             yerr=delta[["sem"]][1],
             new=TRUE,
             col=coluer[1],
             xlab="s.d.",
             ylab="mean(g_high)-mean(g_low)",
             main=paste("Trait:", pars$meta$trait,
                        "Pop. size:",
                        n_samples,
                        "\nErrorbar Reps:",
                        R_REPS),
             xlim=c(min(c(0, sig[["mean"]]-sig[["sem"]])),
                    max(sig[["mean"]]+sig[["sem"]])),
             ylim=c(min(c(0, delta[["mean"]]-delta[["sem"]])),
                    max(delta[["mean"]]+delta[["sem"]])))

    for (i in seq(2, n_group))
        errorbar(sig[["mean"]][i], delta[["mean"]][i],
             xerr=sig[["sem"]][i],
             yerr=delta[["sem"]][i],
             col=coluer[i])

    lines(sig[["mean"]], delta[["mean"]])

    ymin <- min(delta[["mean"]])
    yint <- ymin - sig[["mean"]][delta[["mean"]] == ymin]
    abline(yint, 1)

    legend_names <- vector(mode="character",
                           length=n_group + 1)
    for (i in seq(n_group))
        legend_names[i] <- sprintf("Group Size : %d",
                                   group_sizes[i])
    legend_names[length(legend_names)] <- "Slope = 1"

    legend((max(sig[["mean"]]) - min(sig[["mean"]])) * 0.75
           + min(sig[["mean"]]),
           (max(delta[["mean"]]) - ymin) * 0.3 + ymin,
           legend_names, 
           col=c(coluer, "black"),
           pch=c(rep(19, n_group), -1),,
           lty = c(rep(0, n_group), 1))
    dev.off()


    Z <- rbind(genotypes[idx_low, ],
               genotypes[idx_high,])
    G <- Z %*% t(Z)
    set_coluer <- c(rep(coluer[1], length(idx_low)),
                    rep(coluer[length(coluer)], length(idx_high)))

    pdf(file.path(args[["--out_dir"]], "grm.pdf"))
    heatmap(G,
           Rowv=NA, Colv=NA,
           symm=TRUE, scale="none",
           RowSideColors = set_coluer,
           ColSideColors = set_coluer,
           labRow=NA, labCol=NA,
           main="Breeding value order",
           xlab="Sample", ylab="Sample")
    dev.off()

    low_idx <- as.vector(lower.tri(G))
    pdf(file.path(args[["--out_dir"]], "grm_hist.pdf"))
    par(mar=c(5,5,4,2))
    hist(as.vector(G)[low_idx],
         xlab="Lower triangle elements of G",
         ylab="Occurences",
         cex.lab=1.75,
         cex=1.5,
         main=NULL)
    dev.off()
}


parse_args <- function(args)
{
    parser <- rattaca::argument_parser(
            "--par_file" = rattaca::argument(
                required=TRUE,
                help="path and filename of parameter file"),
            "--plink_genotypes_prefix" = rattaca::argument(
                help=paste("path and filename prefix of plink",
                            "bed, bim, and fam files.",
                            sep=" ")),
            "--out_dir" = rattaca::argument(
                default_val="figs",
                required=FALSE,
                help="directory to print figures."),
            description=paste("Make default plots of",
                              "power analysis")
            )

    return(parser(args))
}



# run main program if R script, otherwise just load
# functions
if (!interactive())
{
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    
    if (is.null(args))
        quit()

    for (arg_item in args)
    {
        if (is.null(arg_item))
            stop(paste("Please include all required arguments,",
                       "see --help."))
    }

    out <- main(args)
}
