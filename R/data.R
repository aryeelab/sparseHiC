#' @include sparseHiC-class.R
NULL

#' hESC1 Sample
#'
#' A sparseHiCdatum object from the ESC replicate 1 
#' @name hESCdatum1
#' @docType data
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
NULL

#' hESC2 Sample
#'
#' A sparseHiCdatum object from the ESC replicate 2
#' 
#' @name hESCdatum2
#' @docType data
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
NULL

#' IMR90 Sample
#'
#' A sparseHiCdatum object from the IMR90 replicate1
#' @name IMR90datum1
#' @docType data
#' 
#' @format A small loops object 
#' \describe{
#'   \item{sampleName}{IMR90_1}
#'   \item{resolutionNamedList}{1Mb with chr20-22 }
#'   \item{metaData}{Shows alignment to hg19}
#' }
#' @return A sparseHiCdatum
#' @source # IMR90_1
#' matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' "HiC-Pro/hIMR90_Rep1/hIMR90_Rep1_1000000_iced.matrix", sep = "/")
#' bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' "HiC-Pro/hIMR90_Rep1/hIMR90_Rep1_1000000_abs.bed", sep = "/")
#' i1 <- import.HiCPro(matrix.files, bed.files, resolutions, sampleName = "IMR90_1", genomeBuild)
#' IMR90datum1 <- getChromosome(i1, c("chr20", "chr21", "chr22"), FALSE)
#' saveRDS(IMR90datum1, "inst/rds/IMR90datum1.rds")
NULL

#' hESC Samples
#'
#' A sparseHiCdata object for both the ESC replicates
#' @name hESCdata
#' @docType data
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
NULL

#' Human/mouse exon locations
#'
#' A dataframe used for plotting annotation for human
#' and mouse. Each loaded .rda has the same variable
#' called "geneinfo" (so don't co-load these), but 
#' the files differ by an m orh
#'
#' @name geneinfo
#' @docType data
#'
#' @format A GRanges object 
#' \describe{
#'   \item{chrom}{Chromosomes without "chr"}
#'   \item{start}{exon start location}
#'   \item{stop}{exon end location}
#'   \item{gene}{Gene Name}
#'   \item{score}{dummy column there for sushi}
#'   \item{strand}{+1 or -1 to indicate side of DNA}
#'   ...
#' }
#' @return A data.frame
#' @source biomaRt July 2015 stable build
NULL
