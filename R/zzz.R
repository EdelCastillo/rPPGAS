#' rPPGAS
#' 
#'  MSI data peak picking
#' 
#' @docType package
#' @author Esteban del Castillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib rPPGAS
#' @name rPPGAS
NULL  

.onUnload <- function (libpath)
{
  library.dynam.unload("rPPGAS", libpath)
}