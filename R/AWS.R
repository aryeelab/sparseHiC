#' @include toDNAlandscapeR.R
NULL

#' Import HiC Sample from AWS Bucket
#'
#' \code{getSampleFromBucket} returns a \code{sparseHiCdatum} object
#' when the user specifies the sample name and the bucket that the object
#' is housed in. If the sample is split by chromosomes (e.g. GM12878 in
#' the DNAlandscapeR bucket), then the function collates those into one
#' object. However, the user can specify only certain chromosomes be 
#' imported using the chr paramter. As a side note, all AWS buckets are
#' lower case letters only, so that may be a source of error when trying
#' to import the bucket. 
#'
#' @param sample The name of the sample (e.g. GM12878 or BT142-rep1)
#' @param bucket Name of the AWS bucket. 
#' @param chr = NA Specify which chromosomes to pull down; if left with
#' NA, then by default, all chromosomes will be pulled down
#' @param organism = "human" or "mouse" likely supported
#' @return A sparseHiCdatum object
#'
#' @examples
#' # Import chromosome 4 from DNAlandscapeR bucket
#' samp <- "hESC1"
#' #hesc <- getSampleFromBucket(samp, "dnalandscaper2", c("chr1", "chr4"))

#' @export
setGeneric(name = "getSampleFromBucket", def = function(sample, bucket, chr = NA, organism = "human")
    standardGeneric("getSampleFromBucket"))

#' @rdname getSampleFromBucket
setMethod("getSampleFromBucket", def= function(sample, bucket, chr = NA, organism = "human") {
    base <- paste0("https://s3.amazonaws.com/", bucket, "/data/", organism, "/hic/", sample, "-HiC")
    resFile <- paste0(base, ".resolutions.txt")
    chrSplit <- as.logical(strsplit(as.character(read.table(resFile, skip = 2)[1,]), split = "=")[[1]][2])
    
    if(chrSplit){ # Separate chromosome
        
        if(!all(is.na(chr))){ # Return specific chromosomes
            md <- data.frame()
            dat <- lapply(chr, function(ct){
                xcon <- gzcon(url(paste0(base, "/", sample, "-HiC", "-", ct, ".rds")))
                z <- readRDS(xcon); close(xcon)
                md <- z@metaData
                z@resolutionNamedList
            })   
            datflat <- unlist(dat)
            ndf <- names(datflat)
            res <- unique(sapply(strsplit(ndf, split = "\\."), "[[", 1))
            lo  <- lapply(res, function(r){
                d <- datflat[grepl(paste0(r, "."), ndf)]
                names <- sapply(strsplit(ndf[grepl(paste0(r, "."), ndf)], split = "\\."), "[[", 2)
                names(d) <- names
                d
            })
            names(lo) <- res
            obj <- new("sparseHiCdatum", sampleName = sample, resolutionNamedList = lo, metaData = md)
        }
        
    } else { # Bound together
        
        xcon <- gzcon(url(paste0(base, ".rds")))
        z <- readRDS(xcon); close(xcon)
        if(!all(is.na(chr))){ # Return specific chromosomes
            nd <- names(z@resolutionNamedList)
            lo <- lapply(z@resolutionNamedList, function(d){
                d[chr]
            })
            names(lo) <- nd
            obj <- new("sparseHiCdatum", sampleName = sample, resolutionNamedList = lo, metaData = z@metaData)
        }
    }
    return(z)
})

#' @include toDNAlandscapeR.R
NULL

#' List Available Hi-C Samples from an AWS Bucket
#'
#' \code{samplesInBucket} returns a character vector of the samples
#' that can be imported from the specified bucket. Should also specify
#' the organism.
#'
#' @param bucket name of the AWS bucket (all lower case)
#' @importFrom aws.s3 get_bucket
#'
#' @examples
#' samplesInBucket("dnalandscaper")

#' @export
setGeneric(name = "samplesInBucket", def = function(bucket)
    standardGeneric("samplesInBucket"))

#' @rdname samplesInBucket
setMethod("samplesInBucket", signature("character"),
          definition = function(bucket) {
              t <- unlist(get_bucket(bucket = bucket))
              full <- t[grepl("hic", t) & grepl("resolutions.txt", t)]
              t2 <- unlist(strsplit(full, split = "-HiC"))
              return(sapply(strsplit(t2[!grepl("resolutions.txt", t2)], "/"), function(sample){
                  v <- sample[4]
                  names(v) <- sample[2]
                  v
              }))
          })
