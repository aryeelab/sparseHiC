#' @include core.R
NULL

#' Compress Hi-C Pro output data
#'
#' \code{sparseCompress.HiCPro} takes output data from HiC-Pro and compresses
#' it. Only intrachromosomal interactions will be preserved. 
#'
#' There are two main modes that can be indicated by the save logic variable.
#' If save = FALSE, this function will return a list of sparse matrices where
#' the names of the matrices correspond to the chromosomes. If save = TRUE,
#' one (if list = TRUE) or many (list = FALSE) .rds files will be created
#' of the compressed Hi-C interactions in sparse matrices. 
#' 
#' @param matrix.file Path to .matrix file from Hi-C Pro Output
#' @param bed.file Path to .bed file from Hi-C Pro Output
#' @param res Character of the resolution of the Hi-C output
#' @param genomeBuild = NA Can specify one of c("hg38", "hg19", "hg18", "mm10", "mm9", "mm8")
#' that are built-in options for Hi-C chromosomes and distances. If not of these options
#' are suitable, then use the \code{manual.chr} and \code{manual.dist} parameters
#' @param n Number of off-diagonal rows to retain features to retain. 
#' If n = 0, retain the full data. 
#' @param out.pre = NA Prefix required if saving to disk
#' @param drop.chrom = c("chrY", "chrM") Chromosomes dropped from compression
#' @param list = TRUE If the user wants to split based on chromosome, must specify
#' that save is TRUE
#' @param save = FALSE Save to .rds in lieu of returning the value? Useful when splitting. 
#' @param dir.create If list is FALSE and save is TRUE, this allows a new folder to be
#' created of the resulting .rds files
#' @param compress if dir.create is TRUE, compresses the directory (to .tgz) and removes
#' the raw data.
#' @param manual.chr = NA Specify a vector of chromosome names in the data
#' @param manual.dist = NA Specify a same length vector as manual.chr with the 
#' chromosomal distances corresponding to each element in the manual.chr
#' @param BPPARAM = bpparam() Parameters to pass to bplapply
#'
#' @return Either .rds files szved on the local disk or a list of sparse matrices
#'
#' @examples
#' # matrix.file <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' # bed.file <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' res <- "1000000"
#' # x <- sparseCompress.HiCPro(matrix.file, bed.file, res, genomeBuild)
#' # x <- sparseCompress.HiCPro(matrix.file, bed.file, res, genomeBuild, list = FALSE,
#' #          save = TRUE, out.pre = "LOL", dir.create = TRUE, compress = TRUE, n = 0)
#' 
#' @import GenomicRanges
#' @importFrom utils read.table combn tar
#' @import Matrix
#' @importFrom BiocParallel bplapply
#' @importFrom readr read_tsv
#' @importFrom reshape2 dcast

#' @export
setGeneric(name = "sparseCompress.HiCPro",
          def = function(matrix.file, bed.file, res, genomeBuild = NA, n = 40,
                         out.pre = NA, drop.chrom = c("chrY", "chrM"), list = TRUE,
                         save = FALSE, dir.create = FALSE, compress = FALSE,
                         manual.chr = NA, manual.dist = NA, BPPARAM = BiocParallel::bpparam())
               standardGeneric("sparseCompress.HiCPro"))

#' @rdname sparseCompress.HiCPro
setMethod(f = "sparseCompress.HiCPro",
          def = function(matrix.file, bed.file, res, genomeBuild = NA, n = 40,
                         out.pre = NA, drop.chrom = c("chrY", "chrM"), list = TRUE,
                         save = FALSE, dir.create = FALSE, compress = FALSE,
                         manual.chr = NA, manual.dist = NA, BPPARAM = BiocParallel::bpparam()){
              
    # Early fails; more in the chrDistBuild function that will serve globally
    if( is.na(out.pre) & save ) stop("Specifiy out.pre for the saved files")
    stopifnot(is.character(res))
    stopifnot(!list & save)
    stopifnot(n >= 0 )
    
    # Configure some of the user parameters
    dist <- chrDistBuild(genomeBuild, manual.chr, manual.dist)
    dist <- dist[!(names(dist) %in% drop.chrom)]
    op <- paste0(out.pre, "_", res)
    fs <- ""
    if(dir.create){
        dir.create(op)
        fs <- paste0(op, "/")
    }
    bed.GRanges <- GRanges(data.frame(
        read_tsv(bed.file, col_names = c("chr", "start", "stop", "region"))))
    dat.long <-read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))

    # Create either a list or save split files
    list.dat <- bplapply(names(dist), function(chr) {
        mat.chrom <- matrixBuild.inter(chr, bed.GRanges, dat.long, res, dist, n)
        if(save & !list){
            saveRDS(Matrix(as.matrix(mat.chrom)),
                    file = paste0(fs, op, "-", chr, ".rds"))
        } else {
            Matrix(as.matrix(mat.chrom))   
        }
    }, BPPARAM = BPPARAM)
    
    if(list) { 
        names(list.dat) <- names(dist)
        if(save) saveRDS(list.dat, file = paste0(fs, op, ".rds"))
        return(list.dat)
    }
    
    if(compress & dir.create) {
        tar(paste0(op,".tgz"), op, compression='gzip')
        unlink(out.pre, recursive = TRUE)
        }
    
})

#' Compress all Hi-C Pro output data, including intrachromosomal interactions
#'
#' \code{sparseCompress.HiCPro.full} takes output data from HiC-Pro and compresses
#' it. All interactions will be preserved including intrachromosomal. No values 
#' can be zeroed out like in the \code{sparseCompress.HiCPro} function. 
#'
#' The option to return a list is not supported since the intrachromosomal data
#' will likely be larger and the number of pairwise chromosome-chromosome interactions
#' scales to O(n^2). 
#' 
#' @param matrix.file Path to .matrix file from Hi-C Pro Output
#' @param bed.file Path to .bed file from Hi-C Pro Output
#' @param res Character of the resolution of the Hi-C output
#' @param genomeBuild = NA Can specify one of c("hg38", "hg19", "hg18", "mm10", "mm9", "mm8")
#' that are built-in options for Hi-C chromosomes and distances. If not of these options
#' are suitable, then use the \code{manual.chr} and \code{manual.dist} parameters
#' @param out.pre File name prefix; will also be name of the new directory created. 
#' "_" + res will be appended to the output file name. 
#' @param drop.chrom = c("chrY", "chrM")
#' @param compress Compress the resulting directory? 
#' @param dir.create Create directory for resulting .rds files? 
#' @param manual.chr = NA Specify a vector of chromosome names in the data
#' @param manual.dist = NA Specify a same length vector as manual.chr with the 
#' chromosomal distances corresponding to each element in the manual.chr
#' @param BPPARAM = bpparam() Parameters to pass to bplapply
#'
#' @return Pairwise .rds files will be saved to the disk of each chromosome interaction, 
#' including interchromosomal interactions
#'
#' @examples
#' # matrix.file <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' # bed.file <- paste(system.file("extdata", package = "processedHiCdata"),
#' # "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' res <- "1000000"
#' # x <- sparseCompress.HiCPro.full(matrix.file, bed.file, res, genomeBuild, "hESC_rep1")
#' 
#' @export
setGeneric(name = "sparseCompress.HiCPro.full",
          def = function(matrix.file, bed.file, res, genomeBuild = NA,
                         out.pre, drop.chrom = c("chrY", "chrM"),
                         dir.create = TRUE, compress = TRUE, 
                         manual.chr = NA, manual.dist = NA,
                         BPPARAM = BiocParallel::bpparam())
               standardGeneric("sparseCompress.HiCPro.full"))

#' @rdname sparseCompress.HiCPro.full
setMethod(f = "sparseCompress.HiCPro.full",
          def = function(matrix.file, bed.file, res, genomeBuild = NA,
                         out.pre, drop.chrom = c("chrY", "chrM"),
                         dir.create = TRUE,compress = TRUE, 
                         manual.chr = NA, manual.dist = NA,
                         BPPARAM = BiocParallel::bpparam()){
              
    # Early fails; more in the chrDistBuild function that will serve globally
    stopifnot(is.character(res))
    
    # Configure some of the user parameters
    dist <- chrDistBuild(genomeBuild, manual.chr, manual.dist)
    dist <- dist[!(names(dist) %in% drop.chrom)]
    op <- paste0(out.pre, "_", res)
    fs <- ""
    if(dir.create){
        dir.create(op)
        fs <- paste0(op, "/")
    }
    bed.GRanges <- GRanges(data.frame(
        read_tsv(bed.file, col_names = c("chr", "start", "stop", "region"))))
    dat.long <-read_tsv(matrix.file, col_names = c("idx1", "idx2", "region"))

    # Create either a list or save split files
    junk <- bplapply(names(dist), function(chr1) {
        lapply(names(dist[seq.int(match(chr1, names(dist)), length(dist), 1)]), function(chr2) {
            mat.chrom <- matrixBuild.intra(chr1, chr2, bed.GRanges, dat.long, res, dist)
            saveRDS(Matrix(as.matrix(mat.chrom)),
                    file = paste0(fs, op, "-", chr1, "-", chr2, ".rds"))
        })
    }, BPPARAM = BPPARAM)
    
    
    if(compress & dir.create) {
        tar(paste0(op,".tgz"), op, compression='gzip')
        unlink(out.pre, recursive = TRUE)
        }
    
})