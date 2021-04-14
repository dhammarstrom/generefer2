#' refset - A set of potential reference genes
#'
#' Raw abundance data of 11 potential reference genes and 2 target genes
#' from a study with data nested within participants.
#'
#'
#' @docType data
#'
#' @usage data(refset)
#'
#' @keywords datasets
#'
#' @references Hammarström et al. (2020) J Physiol, 598: 543-565. https://doi.org/10.1113/JP278455
#' @format A data frame with 3484 rows and 7 variables:
#' \describe{
#'    \item{subject}{participant id}
#'    \item{time}{time-point, w0, week 0 (baseline); w2pre, week 2 before acute exercise; w2post, week 2 after acute exercise; w12, week 12}
#'    \item{condition}{volume condition used in exercises single and multiple set exercises}
#'    \item{target}{gene identifier}
#'    \item{cq}{quantification cycle}
#'    \item{efficiency}{amplification efficiency, averaged over target}
#'    \item{type}{type of target, potential reference genes are "ref", genes of interest are "goi"}
#' }
#'
#' @examples
#'
"refset"
