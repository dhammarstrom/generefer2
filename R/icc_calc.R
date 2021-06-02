#' ICC calculation.
#'
#' This function returns the sum of two values.
#' @param model A lme4 model.
#' @return Numeric.

# function for calculating ICC based on user defined random effects and bootstrap model
icc_calc <- function(model){
        icc <- as.numeric(data.frame(lme4::VarCorr(model))[1,4]) / sum(as.numeric(data.frame(lme4::VarCorr(model))[,4]))
        return(icc)
}
