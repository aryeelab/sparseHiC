#' @include sparseHiC-class.R
NULL

#' hESC1 Sample
#'
#' A sparseHiCdatum object from the ESC replicate 1 
#'
#' @format A small loops object 
#' \describe{
#'   \item{sampleName}{hESC1}
#'   \item{resolutionNamedList}{1Mb with chr20-22 }
#'   \item{metaData}{Shows alignment to hg19}
#' }
#' @return A sparseHiCdatum
#' @source # hESC1
#' matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' e1 <- import.HiCPro(matrix.files, bed.files, resolutions, sampleName = "hESC1", genomeBuild)
#' hESCdatum1 <- getChromosome(e1, c("chr20", "chr21", "chr22"), FALSE)
#' saveRDS(hESCdatum1, "inst/rds/hESCdatum1.rds")
#' 
"hESCdatum1"

#' hESC2 Sample
#'
#' A sparseHiCdatum object from the ESC replicate 2
#'
#' @format A small loops object 
#' \describe{
#'   \item{sampleName}{hESC2}
#'   \item{resolutionNamedList}{1Mb with chr20-22 }
#'   \item{metaData}{Shows alignment to hg19}
#' }
#' @return A sparseHiCdatum
#' @source # hESC2
#' matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' "HiC-Pro/hESC_Rep2/hESC_Rep2_1000000_iced.matrix", sep = "/")
#' bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' "HiC-Pro/hESC_Rep2/hESC_Rep2_1000000_abs.bed", sep = "/")
#' e2 <- import.HiCPro(matrix.files, bed.files, resolutions, sampleName = "hESC1", genomeBuild)
#' hESCdatum2 <- getChromosome(e2, c("chr20", "chr21", "chr22"), FALSE)
#' saveRDS(hESCdatum2, "inst/rds/hESCdatum2.rds")
"hESCdatum2"

#' hESC Samples
#'
#' A sparseHiCdata object for both the ESC replicates
#'
#' @format A small loops object 
#' \describe{
#'   \item{HiCSamplesList}{hESC1 and hESC2 each with 1Mb with chr20-22 }
#'   \item{metaData}{Shows alignment to hg19}
#' }
#' @return A sparseHiCdata object for samples
#' @source # hESC
#' hESCdatum1 <- readRDS("inst/rds/hESCdatum1.rds")
#' hESCdatum2 <- readRDS("inst/rds/hESCdatum2.rds")
#' hESCdata <- hESCdatum1 + hESCdatum2
#' saveRDS(hESCdata, "inst/rds/hESCdata.rds")
"hESCdata"

