#' @include coreCompress.R
NULL

#' Compress Hi-C Pro output data
#'
#' \code{import.HiCPro} takes output data from HiC-Pro and compresses
#' it. Only intrachromosomal interactions will be preserved. 
#'
#' Supply three vectors of equivalent length; one specifies the location
#' of the .matrix files, one the location of the .bed files, and one
#' the resolution of the samples for each element in that order. 
#' 
#' This only works for one sample. Apply the function over multiple samples
#' and then combine with the \code{combine} function to make a single object
#' with multiple samples.
#' 
#' @param matrix.files Path to .matrix files from Hi-C Pro Output
#' @param bed.files Path to .bed files from Hi-C Pro Output
#' @param resolutions Character vector of the resolutions of the Hi-C output
#' @param sampleName Single string of Hi-C data sample that is imported
#' @param genomeBuild = NA Can specify one of c("hg38", "hg19", "hg18", "mm10", "mm9", "mm8")
#' that are built-in options for Hi-C chromosomes and distances. If not of these options
#' are suitable, then use the \code{manual.chr} and \code{manual.dist} parameters
#' @param n = 40 Number of off-diagonal rows to retain features to retain. 
#' If n = 0, retain the full data. 
#' @param drop.chrom = c("chrY", "chrM") Chromosomes dropped from compression
#' @param manual.chr = NA Specify a vector of chromosome names in the data
#' @param manual.dist = NA Specify a same length vector as manual.chr with the 
#' chromosomal distances corresponding to each element in the manual.chr
#' @param tempFile = TRUE Create a temporary file in awk to make i/o faster
#' and more memory efficient. 
#' @param BPPARAM = bpparam() Parameters to pass to bplapply
#'
#' @return sparseHiCdatum of all intrachromosomal interactions of n diagonals
#'
#' @examples
#' # matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' # bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' resolutions <- "1000000"
#' sampleName = "hESC1"
#' # x <- import.HiCPro(matrix.files, bed.files, resolutions, sampleName = "hESC1", genomeBuild)
#' 
#' @import GenomicRanges
#' @importFrom utils read.table combn tar
#' @import Matrix
#' @importFrom BiocParallel bplapply
#' @importFrom readr read_tsv
#' @importFrom reshape2 acast dcast
#' @importFrom stats aggregate

#' @export
setGeneric(name = "import.HiCPro",
          def = function(matrix.files, bed.files, resolutions, sampleName, genomeBuild = NA, 
                         n = 40, drop.chrom = c("chrY", "chrM"), manual.chr = NA, manual.dist = NA,
                         tempFile = TRUE, BPPARAM = BiocParallel::bpparam())
               standardGeneric("import.HiCPro"))

#
#matrix.files="/data/aryee/bernstein/hic/runs/hicpro_output/hic_results/matrix/HCT116-rao-1/raw/100000/HCT116-rao-1_100000.matrix"
#bed.files="/data/aryee/bernstein/hic/runs/hicpro_output/hic_results/matrix/HCT116-rao-1/raw/100000/HCT116-rao-1_100000_abs.bed"
#resolutions="100000"
#sampleName="prueba"
#genomeBuild="hg19"
#drop.chrom <- c("chrY", "chrM")
#manual.chr=NA
#manual.dist=NA
#tempFile=TRUE
#BPPARAM <- BiocParallel:::SerialParam()


#' @rdname import.HiCPro
setMethod(f = "import.HiCPro",
          def = function(matrix.files, bed.files, resolutions, sampleName, genomeBuild = NA,
                         n = 40, drop.chrom = c("chrY", "chrM"), manual.chr = NA, manual.dist = NA,
                         tempFile = TRUE, BPPARAM = BiocParallel::bpparam()){
              
    # Early fails; more in the chrDistBuild function that will serve globally
    stopifnot(is.character(resolutions))
    stopifnot( n >= 0 )
    stopifnot(length(matrix.files) == length(bed.files) & length(bed.files) == length(resolutions))
    
    # Configure some of the user parameters
    dist <- .chrDistBuild(genomeBuild, manual.chr, manual.dist)
    dist <- dist[!(names(dist) %in% drop.chrom)]
    # Import data for each supplied resolution
    collectedRes <- lapply(1:length(resolutions), function(i){
        bed.GRanges <- GRanges(data.frame(
            read_tsv(bed.files[i], col_names = c("chr", "start", "stop", "region"), col_types="ciii")))
        matrix.file <- matrix.files[i]
        if(tempFile & n != 0 ){
            temp <- paste0(sampleName, ".tempFile.sparseHiC.txt")
            cmd <- paste0('awk \'$1+', as.character(n),' >= $2 {print $0}\' ', 
                          matrix.file, ' > ', temp)
            system(cmd) 
            matrix.file <- temp
        }
        dat.long <-read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))
        if( tempFile & n != 0 ) file.remove(temp)
        list.dat <- bplapply(names(dist), function(chr) {
            Matrix(.matrixBuild.intra2(chr, bed.GRanges, dat.long, resolutions[i], dist, n))
        }, BPPARAM = BPPARAM)
        
        names(list.dat) <- names(dist)
        list.dat
    })
    names(collectedRes) <- resolutions
    
    if(is.na(genomeBuild)) genomeBuild <- "custom"
    md <- data.frame(genomeBuild)
    obj <- new("sparseHiCdatum", sampleName = sampleName, resolutionNamedList = collectedRes, metaData = md)
    return(obj)
})


#' Compress all Hi-C Pro output data, including interchromosomal interactions
#'
#' \code{import.HiCPro.full} takes output data from HiC-Pro and compresses
#' it. All interactions will be preserved including interchromosomal. No values 
#' can be zeroed out like in the \code{import.HiCPro} function. Not recommended
#' function to utilize for memory disk concerns but we do support importing the
#' full data into this framework. Additionally, a temporary file cannot be generated
#' to facilite memory constraints, so this command may fail for high resolution data. 
#' 
#' @param matrix.files Paths to .matrix files from Hi-C Pro Output
#' @param bed.files Paths to .bed files from Hi-C Pro Output
#' @param resolutions Character vector of the resolutions of the Hi-C output
#' @param sampleName Single string of Hi-C data sample that is imported
#' @param genomeBuild = NA Can specify one of c("hg38", "hg19", "hg18", "mm10", "mm9", "mm8")
#' that are built-in options for Hi-C chromosomes and distances. If not of these options
#' are suitable, then use the \code{manual.chr} and \code{manual.dist} parameters
#' @param drop.chrom = c("chrY", "chrM")
#' @param manual.chr = NA Specify a vector of chromosome names in the data
#' @param manual.dist = NA Specify a same length vector as manual.chr with the 
#' chromosomal distances corresponding to each element in the manual.chr
#' @param BPPARAM = bpparam() Parameters to pass to bplapply
#'
#' @return sparseHiCdatum of all chromosomal interactions pairwise
#'
#' @examples
#' # matrix.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' # bed.files <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' resolutions <- "1000000"
#' sampleName = "hESC1"
#' # x <- import.HiCPro.full(matrix.file, bed.file, res, sampleName = "hESC1", genomeBuild)
#' 
#' @export
setGeneric(name = "import.HiCPro.full",
          def = function(matrix.files, bed.files, resolutions, sampleName, genomeBuild = NA,
                         drop.chrom = c("chrY", "chrM"), manual.chr = NA, manual.dist = NA,
                         BPPARAM = BiocParallel::bpparam())
              standardGeneric("import.HiCPro.full"))

#' @rdname import.HiCPro.full
setMethod(f = "import.HiCPro.full",
          def = function(matrix.files, bed.files, resolutions, sampleName, genomeBuild = NA,
                         drop.chrom = c("chrY", "chrM"), manual.chr = NA, manual.dist = NA,
                         BPPARAM = BiocParallel::bpparam()){
              
    # Early fails; more in the chrDistBuild function that will serve globally
    stopifnot(is.character(resolutions))
    stopifnot(length(matrix.files) == length(bed.files) & length(bed.files) == length(resolutions))
    
    # Configure some of the user parameters
    dist <- .chrDistBuild(genomeBuild, manual.chr, manual.dist)
    dist <- dist[!(names(dist) %in% drop.chrom)]
    
    # Import data for each supplied resolution
    collectedRes <- lapply(1:length(resolutions), function(i){
        bed.GRanges <- GRanges(data.frame(
            read_tsv(bed.files[i], col_names = c("chr", "start", "stop", "region"))))
        matrix.file <- matrix.files[i]
        dat.long <-read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))
        
        olists <- bplapply(names(dist), function(chr1) {
            chr2s <- names(dist[seq.int(match(chr1, names(dist)), length(dist), 1)])
            ilist <- lapply(chr2s, function(chr2) {
                mat.chrom <- .matrixBuild.inter(chr1, chr2, bed.GRanges, dat.long, resolutions[i], dist)
                Matrix(as.matrix(mat.chrom))
            })
            names(ilist) <- paste0(chr1, "-", chr2s)
            ilist
        }, BPPARAM = BPPARAM)
        
        list.dat <- unlist(olists, recursive = FALSE)
    })
    names(collectedRes) <- resolutions
    if(is.na(genomeBuild)) genomeBuild <- "custom"
    md <- data.frame(genomeBuild)
    obj <- new("sparseHiCdatum", sampleName = sampleName, resolutionNamedList = collectedRes, metaData = md)
    return(obj)
})
