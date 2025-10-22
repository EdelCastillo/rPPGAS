#' rPPGAS
#' 
#'  MSI data peak picking
#' 
#' @name rPPGAS
#' @aliases rPPGAS-package 
#' @author Esteban del Castillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib rPPGAS
NULL  

.onUnload <- function (libpath)
{
  library.dynam.unload("rPPGAS", libpath)
}