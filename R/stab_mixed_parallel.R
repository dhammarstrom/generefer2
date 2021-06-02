#' Mixed model stability parallel computation
#'
#' This function returns a data frame of ICC calculated from gene sets of size n. To be used with the
#' exported wrapper function.
#'
#' @param data, a tidy data set provided in the wrapper function
#' @param combinations, a set of combination as created by the wrapper function
#' @param cores, number of cores to be used in parallel execution
#' @param progress, should a progress bar be shown
#' @param response, the response variable (e.g. log(expression))
#' @param target, character vector for target column in data
#' @param fixed, formula for the fixed effects (not including target)
#' @param random, formula in lme4-style for the random effect (only one accepted)
#' @param p.threshold, threshold for effect of Beta (design variables)
#' @param icc.interval, specifying interval for the ICC CI
#' @param icc.type, type of ICC passed to boot function
#' @param n.sims, number of simulations in the ICC boot estimation
#'
#' @return A data frame with cominations of genes and their ICC confidence interval.


stab_mixed_parallel <- function(data,
                                combinations,
                                cores,
                                progress,
                                response,
                                target,
                                fixed,
                                random,
                                p.threshold,
                                icc.interval,
                                icc.type,
                                n.sims) {


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


        ### Create formulas #######################################################
        response.target <- paste0(response, "~", target)
        full.f <- paste(response.target, fixed, sep="*", collapse = "+")

        # combine to formula
        full.f <- as.formula(paste(full.f, random, sep = "+"))
        reduced.f <- as.formula(paste(response.target, random, sep = "+"))





        results <- foreach::foreach(i = 1:ncol(combinations),
                                    .errorhandling='pass',
                            .packages = c("lme4", "dplyr"),
                            .options.snow = opts) %dopar% {




                 tryCatch(
                        expr = {




                                    # A result data frame for icc
                                    results.icc <- data.frame(gene.combination = NA,
                                                              icc    = NA ,
                                                              icc.l  = NA ,
                                                              icc.u  = NA )



                                    # A data frame to store p-values from global test
                                    results.lrt <- data.frame(p.val =  NA)





                                    genes <- combinations[,i] # genes in the specified combination
                                    results.icc[1,1] <- paste0(combinations[,i], collapse = ":") # gene combination
                                    subset.data <- data[data[, which(colnames(data) == target)] %in% genes,] # subset of genes used in iteration



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

                                            b <- lme4::bootMer(icc.mod, icc_calc, use.u = FALSE, nsim = n.sims)

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
                                    r$n.genes <- nrow(combinations)

                                    r

                        },
                        error = function(e){
                                message('** ERR at ', Sys.time(), " **")
                                print(e)


                        })
                }

        close(pb)
        parallel::stopCluster(cl)

        # Test if all results are data frames, keep only combinations with results



        is_df <- sapply(results, is.data.frame)

        finaldf <- bind_rows(results[is_df])

        return(finaldf)
} # end function
