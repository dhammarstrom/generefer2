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

        results <- list()

        data <- data.frame(data)

        for(n in 1:length(n.genes)){

                # Store results from each n.genes
                results_combined <- list()

                # possible combinations of size n.genes
                combinations <- combn(unique(data[, which(colnames(data) == target)]), n.genes[n])

                # function for calculating ICC based on user defined random effects and bootstrap model
                icc.calc <- function(model){
                        icc <- as.numeric(data.frame(lme4::VarCorr(model))[1,4]) / sum(as.numeric(data.frame(lme4::VarCorr(model))[,4]))
                        return(icc)
                }



                # set up parallel processing using the foreach package
                if(is.null(cores)) {
                        cores_detected <- parallel::detectCores()
                        cl <- parallel::makeCluster(cores_detected[1])
                        doSNOW::registerDoSNOW(cl)
                }

                if(!is.null(cores)) {
                        if(!is.numeric(cores)) stop("Number of cores must be specified as an integer or NULL")
                        cl <- parallel::makeCluster(cores)
                        doSNOW::registerDoSNOW(cl)
                }



                # Progress bar (should work on linux also)
                if(progress==TRUE){
                        iterations <- ncol(combinations)
                        pb <- txtProgressBar(max = iterations, style = 3)
                        progress <- function(n) setTxtProgressBar(pb, n)
                        opts <- list(progress = progress)
                } else {

                        iterations <- ncol(combinations)
                        opts <- list()

                }



                results <- foreach::foreach(i = 1:ncol(combinations),
                                   .packages = c("lme4", "dplyr"),
                                   .options.snow = opts) %dopar% {



                         # A result data frame for icc
                         results.icc <- data.frame(gene.combination = NA,
                                                   icc    = NA ,
                                                   icc.l  = NA ,
                                                   icc.u  = NA )



                         # A data frame to store p-values from global test
                         results.lrt <- data.frame(p.val =  NA)





                        genes <- combinations[,i] # genes in the specified combination
                        results.icc[1,1] <- paste(genes[1:n.genes[n]], collapse=":") # gene combination
                        subset.data <- data[data[, which(colnames(data) == target)] %in% genes,] # subset of genes used in iteration



                                response.target <- paste0(response, "~", target)
                                full.f <- paste(response.target, fixed.effects, sep="*", collapse = "+")

                                # combine to formula
                                full.f <- as.formula(paste(full.f, random.effect, sep = "+"))
                                reduced.f <- as.formula(paste(response.target, random.effect, sep = "+"))

                                reduced.mod <- lme4::lmer(reduced.f, data = subset.data,
                                                          REML = FALSE,
                                                          control = lme4::lmerControl(calc.derivs = FALSE))
                                full.mod <- lme4::lmer(full.f, data=subset.data,
                                                       REML=FALSE,
                                                       control = lme4::lmerControl(calc.derivs = FALSE))
                                results.lrt[1, 1] <- data.frame(anova(full.mod, reduced.mod))[2,8]



                                if(any(results.lrt[1, 1] < p.threshold)) h.test <- TRUE else h.test <- FALSE



                        # If an hypothesis test is significant
                        # do not calculate bootstrap ICC
                        if(h.test == TRUE) {
                                results.icc[1,2]<-NA
                                results.icc[1,3]<-NA
                                results.icc[1,4]<-NA
                        } else {
                                # If no t-values over critical t
                                # calculate bootstrap ICC

                                # fit model for ICC calculations
                                icc.mod <- lme4::lmer(reduced.f,
                                                      REML = TRUE,
                                                      control = lme4::lmerControl(calc.derivs = FALSE),
                                                      data = subset.data)

                                b <- lme4::bootMer(icc.mod, icc.calc, use.u = FALSE, nsim = n.sims)

                                ci <- boot::boot.ci(b, conf = icc.interval, type = icc.type)

                                if(icc.type=="basic"){
                                        cis<-c(data.frame(ci[4])[1,4], data.frame(ci[4])[1,5])
                                }
                                if(icc.type=="norm"){
                                        cis<-c(data.frame(ci[4])[1,2], data.frame(ci[4])[1,3])
                                }
                                if(icc.type=="perc"){
                                        cis<-c(data.frame(ci[4])[1,4], data.frame(ci[4])[1,5])
                                }

                                results.icc[1, 2] <- b$t0
                                results.icc[1, 3] <- cis[1]
                                results.icc[1, 4] <- cis[2]

                        }

                                r <- cbind(results.icc, results.lrt)
                                r$n.genes <- n.genes[n]

                                r

                }


        results_combined[[i]] <- dplyr::bind_rows(results)

        }

        return(dplyr::bind_rows(results_combined))

}
