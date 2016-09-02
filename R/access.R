#' @include sparseHiC-class.R
NULL

#' Grab a sample's specific resolution
#'
#' \code{getResolution} takes a \code{sparseHiCdatum} object
#' and returns a list of sparse Hi-C matrices. If a list is 
#' returned, then the names are the chromosomes. 
#'
#' @param obj A \code{sparseHiCdatum} object
#' @param res Resolution(s) desired
#' @param a.list = TRUE Return a list of matrices? If FALSE, return another
#' \code{sparseHiCdatum} object
#'
#' @return Returns object subsetted by resolution
#'
#' @examples
#' resolutions <- "1000000"

#' @export
setGeneric(name = "getResolution", def = function(obj, res, a.list = TRUE)
    standardGeneric("getResolution"))

#' @rdname getResolution
setMethod("getResolution", signature("sparseHiCdatum", "ANY", "ANY"),
          definition = function(obj, res, a.list = TRUE) {
              options(scipen=999)
              res <- as.character(res)
              stopifnot(res %in% names(obj@resolutionNamedList))
              if(a.list){
                  return(obj@resolutionNamedList[[res]])
              } else {
                    obj@resolutionNamedList <- obj@resolutionNamedList[res]
                    return(obj)
              }
})

#' Grab a specific chromosome
#'
#' \code{getChromosome} takes a \code{sparseHiCdatum} object
#' and returns a list of sparse Hi-C matrices. If a list is
#' returned, then the names are the resolutions. If the 
#' length of chr > 1 (multiple chromosomes), then
#' a list of lists will be returned if a.list is TRUE
#'
#' @param obj A \code{sparseHiCdatum} object
#' @param chr Chromosome(s) desired
#' @param a.list = TRUE Return a list of matrices? If FALSE,
#' return another \code{sparseHiCdatum} object
#'
#' @return Returns object subsetted by resolution
#'
#' @examples
#' chr <- "chr7"

#' @export
setGeneric(name = "getChromosome", def = function(obj, chr, a.list = TRUE)
    standardGeneric("getChromosome"))

#' @rdname getChromosome
setMethod("getChromosome", signature("sparseHiCdatum","character", "ANY"),
          definition = function(obj, chr, a.list = TRUE) {
              stopifnot(chr %in% names(obj@resolutionNamedList[[1]]))
              simp <- lapply(obj@resolutionNamedList, function(x){
                  x[names(x) %in% chr]
              })
              if(a.list & length(chr) > 1){
                return(simp)
              } else if (a.list & length(chr) == 1) {
                  ns <- names(simp)
                  out <- unlist(simp)
                  names(out) <- ns
                  return(out)
              } else {
                obj@resolutionNamedList <- simp
                return(obj)
              }
})

#' Grab a specific sparse Hi-C matrix
#'
#' \code{getHiCMatrix} takes a \code{sparseHiCdatum} object
#' and returns a single matrix. Thus, the length of chr
#' and res must be 1
#'
#' @param obj A \code{sparseHiCdatum} object
#' @param chr Chromosome desired
#' @param res Resolution desired
#'
#' @return Returns a Hi-C Matrix
#'
#' @examples
#' chr <- "chr7"
#' res <- "1000000"

#' @export
setGeneric(name = "getHiCMatrix", def = function(obj, chr, res)
    standardGeneric("getHiCMatrix"))

#' @rdname getHiCMatrix
setMethod("getHiCMatrix", signature("sparseHiCdatum", "ANY", "ANY"),
          definition = function(obj, chr, res) {
              stopifnot(length(chr) == 1)
              stopifnot(length(res) == 1)
              return(obj@resolutionNamedList[[res]][[chr]])
})