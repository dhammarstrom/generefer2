#' Analyze reference gene stability in a mixed model framework -- Alternative function
#'
#' The function is an implementation of the method proposed by Dai, H., et al. (2013) with modifications.
#' The function determines gene stability by iterative fitting of potential gene combinations and calculation of a bootstrap confidence interval
#' around the intraclass correlation. Inference about fixed effects are by default drawn from LRT-tests based where
#' the full model is compared to a reduced model not containing the effects of interest.
#'
#' @references Dai, H., et al. (2013). "Mixed modeling and sample size calculations for identifying housekeeping genes." Stat Med 32(18): 3115-3125.
#' @param data a data.frame() in tidy format, each row an observation and each column a variable
#' @param target character name of column containing gene/target identifier
#' @param response character Specifying the respons column, typically log-transformed expression values.
#' @param fixed.effects character vector specifying fixed effects included in the model
#' @param random.effect character vector specifying the random effect term, should be specified as used in lme4 e.g. (1|participant)
#' @param p.threshold numeric Specifying the p-value threshold for hypothesis testing in LRT, default is 0.05.
#' @param icc.interval numeric a single fraction specifying the bootstrap confidence interval, defaults to 0.95
#' @param icc.type character specifying which method to use for bootstrap confidence interval calculation
#' ("norm", "basic" or "perc"), default = "norm". See boot::boot.ic for details.
#' @param n.genes numeric A vector specifying the number of genes in possible combinations of genes to be evaluated.
#' Can be a vector of e.g. c(2,3), which will perform the algorithm using combinations of two and three genes/targets.
#' @param n.sims numeric Specifies how many bootstraps to perform for the calculation of CI. This process is time consuming when large numbers are used. Defaults to 500 simulations.
#' @param cores Number of cores used in parallel execution, if NULL, all avalable cores are used.
#' @param progress Logical default to TRUE gives a progressbar.
#' @return A data frame with with bootstrap confidence intervals for each possible combination of n.genes without significant fixed effects.
#' @import foreach
#' @export
#'
stab_mixed <- function(data,
                        target,
                        response,
                        fixed.effects,
                        random.effect,
                        p.threshold = 0.05,
                        icc.interval = 0.95,
                        icc.type="norm",
                        n.genes=2,
                        n.sims=500,
                        cores = NULL,
                       progress=TRUE){

        data <- data.frame(data)

        # Store results from each n.genes
        results_combined <- list()

        for(n in 1:length(n.genes)){



                # possible combinations of size n.genes
                combinations <- combn(unique(data[, which(colnames(data) == target)]), n.genes[n])


                results_combined[[n]] <-  stab_mixed_parallel(data = data,
                                                              combinations = combinations,
                                                              cores = cores,
                                                              progress = progress,

                                                              response = response,
                                                              target = target,
                                                              fixed = fixed.effects,
                                                              random = random.effect,
                                                              p.threshold = p.threshold,
                                                              icc.interval = icc.interval,
                                                              icc.type = icc.type,
                                                              n.sims = n.sims)



        }

        return(dplyr::bind_rows(results_combined))

}
