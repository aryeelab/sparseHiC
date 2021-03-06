#' @include hicpro.R
NULL

#' Export sparseHiCdatum to DNAlandscapeR files
#'
#' \code{exportDNAlandscapeR} takes a \code{sparseHiCdatum} object
#' and produces two files, a .rds and a .sparseHiC.meta file
#' that can be used in DNAlandscapeR (dnalandscaper.aryeelab.org). If
#' splitByChr is specified as TRUE, will create a folder with the sample 
#' name and an .rds of each chromosome for faster i/o. These split files
#' can only be linked via a bucket in DNAlandscapeR.
#'
#' @param dat A \code{sparseHiCdatum} object
#' @param newSampleName = "" New sample name for the object, if desired. Current
#' one in the object is used by default.
#' @param out.dir = "" The output directory to save the new tarball
#' @param splitByChr = FALSE  Split .rds files by chromosome for faster i/o in browser?
#'
#' @return Will produce two files in the specified out.dir
#'
#' @examples
#' # matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' #"HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' #bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' #"HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' resolutions <- "1000000"
#' #exportDNAlandscapeR(matrix.files, bed.files, resolutions,
#' #sampleName = "hESC1", genomeBuild, splitByChr = TRUE)


#' @export
setGeneric(name = "exportDNAlandscapeR", def = function(dat, newSampleName = "", out.dir = "", splitByChr = FALSE)
    standardGeneric("exportDNAlandscapeR"))

#' @rdname exportDNAlandscapeR
setMethod("exportDNAlandscapeR", c("sparseHiCdatum"), function(dat, newSampleName = "", out.dir = "", splitByChr = FALSE) {
    
    if(newSampleName == "") newSampleName <- dat@sampleName  
    
    # Set up the output for the metadata
    ss <- paste0(newSampleName, "-HiC")
    op <- paste0(out.dir, ss)
    dir.create(op)
    res <- names(dat@resolutionNamedList)
    warn <- "Do not edit this file; keep it in the same directory as the .rds file(s) for visualizing in DNAlandscapeR"
    sbc <- paste0("splitByChr=", as.character(splitByChr))
    outText <- c(warn, ss, sbc, paste(res, collapse = ","))
    
    if(splitByChr) {
        for(chr in names(dat@resolutionNamedList[[1]])){
            nd <- names(dat@resolutionNamedList)
            lo <- lapply(dat@resolutionNamedList, function(d){
                d[chr]
            })
            names(lo) <- nd
            rdsFile<-file(paste0(op, "/", ss, "-", chr, ".rds"))
            obj <- new("sparseHiCdatum", sampleName = newSampleName, resolutionNamedList = lo, metaData = dat@metaData)
            saveRDS(obj, file = rdsFile)
            close(rdsFile)
        }
    } else {
        rdsFile<-file(paste0(op, "/", ss, ".rds"))
        saveRDS(dat, file = rdsFile)
    }    
    
    # Output .txt file
    fileConn<-file(paste0(op, "/", ss, ".sparseHiC.meta"))
    writeLines(outText, fileConn)
    close(fileConn)
    
    tar(paste0(op,".tgz"),op,compression='gzip'); unlink(op, recursive = TRUE)
})


#' Export DNAlandscapeR files straight from HiC Pro Output
#'
#' \code{HiCPro2DNAlandscapeR} takes the input files from
#' HiC Pro output and creates a .rds and a .resolutions.txt
#' set of files that can be visualized in
#' DNAlandscapeR (dnalandscaper.aryeelab.org). If
#' splitByChr is specified as TRUE, will create a folder with the sample 
#' name and an .rds of each chromosome for faster i/o. These split files
#' can only be linked via a bucket in DNAlandscapeR.
#'
#' @param matrix.files Path to .matrix files from Hi-C Pro Output
#' @param bed.files Path to .bed files from Hi-C Pro Output
#' @param resolutions Character vector of the resolutions of the Hi-C output
#' @param sampleName Single string of Hi-C data sample that is imported. "-HiC" will 
#' automatically be appended.
#' @param genomeBuild = NA Can specify one of c("hg38", "hg19", "hg18", "mm10", "mm9", "mm8")
#' that are built-in options for Hi-C chromosomes and distances. If not of these options
#' are suitable, then use the \code{manual.chr} and \code{manual.dist} parameters
#' @param splitByChr = FALSE Split .rds files by chromosome for faster i/o in browser?
#' @param out.dir = "" The output directory to save the new tarball
#' @param drop.chrom = c("chrY", "chrM") Chromosomes dropped from compression
#' @param manual.chr = NA Specify a vector of chromosome names in the data
#' @param manual.dist = NA Specify a same length vector as manual.chr with the 
#' chromosomal distances corresponding to each element in the manual.chr
#' @param tempFile = TRUE Create a temporary file in Awk to make i/o faster
#' and more memory efficient. 
#' @param BPPARAM = bpparam() Parameters to pass to bplapply
#'
#' @return Will produce two files in the specified out.dir
#'
#' @examples
#' #' # matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' # bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' resolutions <- "1000000"
#' sampleName = "hESC1"
#' #HiCPro2DNAlandscapeR(matrix.files, bed.files, resolutions, sampleName = "hESC1", genomeBuild)
#' #HiCPro2DNAlandscapeR(matrix.files, bed.files, resolutions,
#' #sampleName = "hESC1", genomeBuild, splitByChr = TRUE)

#' @export
setGeneric(name = "HiCPro2DNAlandscapeR", 
           def = function(matrix.files, bed.files, resolutions, sampleName, genomeBuild = NA, splitByChr = FALSE,
                          out.dir = "", drop.chrom = c("chrY", "chrM"), manual.chr = NA, manual.dist = NA,
                          tempFile = TRUE, BPPARAM = BiocParallel::bpparam())
               standardGeneric("HiCPro2DNAlandscapeR"))

#' @rdname HiCPro2DNAlandscapeR
setMethod("HiCPro2DNAlandscapeR",
         def = function(matrix.files, bed.files, resolutions, sampleName, genomeBuild = NA, splitByChr = FALSE,
                          out.dir = "", drop.chrom = c("chrY", "chrM"), manual.chr = NA, manual.dist = NA,
                          tempFile = TRUE, BPPARAM = BiocParallel::bpparam()) {
    
    
    dat <- import.HiCPro(matrix.files, bed.files, resolutions, sampleName, genomeBuild,
                         n = 40, drop.chrom, manual.chr, manual.dist, tempFile, BPPARAM)
    # Set up the output for the metadata
    ss <- paste0(sampleName, "-HiC")
    op <- paste0(out.dir, ss)
    dir.create(op)
    warn <- "Do not edit this file; keep it in the same directory as the .rds file(s) for visualizing in DNAlandscapeR"
    sbc <- paste0("splitByChr=", as.character(splitByChr))
    outText <- c(warn, ss, sbc, paste(resolutions, collapse = ","))
    
    if(splitByChr) {
        
        for(chr in names(dat@resolutionNamedList[[1]])){
            nd <- names(dat@resolutionNamedList)
            lo <- lapply(dat@resolutionNamedList, function(d){
                d[chr]
            })
            names(lo) <- nd
            rdsFile<-file(paste0(op, "/", ss, "-", chr, ".rds"))
            obj <- new("sparseHiCdatum", sampleName = sampleName, resolutionNamedList = lo, metaData = dat@metaData)
            saveRDS(obj, file = rdsFile)
        }
    } else {
        rdsFile<-file(paste0(op, "/", ss, ".rds"))
        saveRDS(dat, file = rdsFile)
        close(rdsFile)
    }  
    
    # Output .txt file
    fileConn<-file(paste0(op, "/", ss, ".sparseHiC.meta"))
    writeLines(outText, fileConn)
    close(fileConn)
    
    tar(paste0(op,".tgz"),op,compression='gzip'); unlink(op, recursive = TRUE)
})