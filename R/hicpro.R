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
#' @param n Number of off-diagonal rows to retain features to retain. 
#' If n = 0, retain the full data. 
#' @param genomeBuild = NA Can specify one of c("hg19", "hg18", "mm9") that are
#' built-in options for Hi-C chromosomes and distances. If not of these options
#' are suitable, then use the \code{manual.chr} and \code{manual.dist} parameters
#' @param out.pre = NA Prefix required if 
#' @param drop.chrom = c("chrY", "chrM")
#' @param list = TRUE
#' @param save = FALSE
#' @param compress if dir.create is TRUE, compresses the directory (to .tgz) and removes
#' the raw data. 
#' @param dir.create If list is FALSE and save is TRUE
#' @param manual.chr = NA Specify a vector of chromosome names in the data
#' @param manual.dist = NA Specify a same length vector as manual.chr with the 
#' chromosomal distances corresponding to each element in the manual.chr
#' @param BPPARAM = bpparam() Parameters to pass to bplapply
#'
#' @return Either .rds files szved on the local disk or a list of sparse matrices
#'
#' @examples
#' matrix.file <- paste(system.file("extdata", package = "processedHiCData"), "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
#' bed.file <- paste(system.file("extdata", package = "processedHiCData"), "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
#' genomeBuild <- "hg19"
#' res <- "1000000"
#' # x <- sparseCompress.HiCPro(matrix.file, bed.file, res, genomeBuild)
#' 
#' @import GenomicRanges
#' @import Matrix
#' @importFrom BiocParallel bplapply
#' @importFrom readr read_tsv
#' @importFrom reshape2 dcast

#' @export
setGeneric(name = "sparseCompress.HiCPro",
          def = function(matrix.file, bed.file, res, n = 40, genomeBuild = NA,
                         out.pre = NA, drop.chrom = c("chrY", "chrM"), list = TRUE,
                         save = FALSE, compress = FALSE, dir.create = FALSE,
                         manual.chr = NA, manual.dist = NA, BPPARAM = bpparam())
               standardGeneric("sparseCompress.HiCPro"))

#' @rdname sparseCompress.HiCPro
setMethod(f = "sparseCompress.HiCPro",
          def = function(matrix.file, bed.file, res, n = 40, genomeBuild = NA,
                         out.pre = NA, drop.chrom = c("chrY", "chrM"), list = TRUE,
                         save = FALSE, compress = FALSE, dir.create = FALSE,
                         manual.chr = NA, manual.dist = NA, BPPARAM = bpparam()){
              
    # Early fails; more in the chrDistBuild function that will serve globally
    if( is.na(out.pre) & save ) stop("Specifiy out.pre for the outputted .rds files")
    stopifnot(is.character(res))
    stopifnot(is.integer(n), n >= 0 )
    
    # Configure some of the user parameters
    dist <- chrDistBuild(genomeBuild, manusal.chr, manual.dist)
    dist <- dist[!(names(dist) %in% drop.chrom)]
    op <- paste0(out.pre, "_", res)
    if(dir.create) dir.create(op)
    
    # Create either a list or save split files
    list.dat <- bplapply(names(dist), function(chr) {
        mat.chrom <- matrixBuild(chr, bed.GRanges, dat.long, dist, n = 41)
        if(save & !list){
            saveRDS(Matrix(as.matrix(mat.chrom)),
                    file = paste0(op, "/", op, "-", chr, ".rds"))
        } else {
            Matrix(as.matrix(mat.chrom))   
        }
    }, BPPARAM = param)
    
    if(list) { 
        names(list.dat) <- names(dist)
        if(save) saveRDS(list.dat, file = paste0(op, ".rds"))
        return(list.dat)
    }
    
    if(compress & dir.create) {
        tar(paste0(op,".tgz"), op, compression='gzip')
        unlink(out.pre, recursive = TRUE)
        }
    
})