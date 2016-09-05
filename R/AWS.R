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
#' samp <- "IMR90"
#' # hesc <- getSampleFromBucket(samp, "dnalandscaper", c("chr1", "chr4"))

#' @export
setGeneric(name = "getSampleFromBucket", def = function(sample, bucket, chr = NA, organism = "human")
    standardGeneric("getSampleFromBucket"))

#' @rdname getSampleFromBucket
setMethod("getSampleFromBucket", def= function(sample, bucket, chr = NA, organism = "human") {
    if(length(sample) > 1) stop("use getSamplesFromBucket (plural) to get multiple samples from the same bucket")
    extURL <- paste0("data/", organism, "/hic/", sample, "-HiC/", sample, "-HiC")
    base <- paste0("https://s3.amazonaws.com/", bucket, "/", extURL)
    resFile <- paste0(base, ".sparseHiC.meta")
    chrSplit <- as.logical(strsplit(as.character(read.table(resFile, skip = 2)[1,]), split = "=")[[1]][2])
    
    if(chrSplit){ # Separate chromosome
        
        # User didn't define chromosomes, so import them all
        if(all(is.na(chr))) {
            t <- gsubfn::strapplyc(RCurl::getURL((paste0("https://s3.amazonaws.com/", bucket))), "<Key>(.*?)</Key>", simplify = c)
            chr <- sapply(strsplit(t[grepl(extURL,t) & grepl(".rds", t) & grepl("chr",t)], split = extURL),
                   function(hit){ strsplit(strsplit(hit[[2]], ".rds")[[1]], "-")[[1]][2] })
            chr <- .chrOrder(chr)
        }
            
        md <- data.frame()
        dat <- lapply(chr, function(ct){
            xcon <- gzcon(url(paste0(base, "-", ct, ".rds")))
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
    return(obj)
})

#' Import multiple Hi-C samples from AWS Bucket
#'
#' \code{getSamplesFromBucket} returns a \code{sparseHiCdata} object
#' when the user specifies multiple sample names present in the bucket
#'
#' @param samples The names of the samples 
#' @param bucket Name of the AWS bucket. 
#' @param chr = NA Specify which chromosomes to pull down; if left with
#' NA, then by default, all chromosomes will be pulled down
#' @param organism = "human" or "mouse" likely supported
#' @return A sparseHiCdata object
#'
#' @examples
#' # Import chromosome 4 from DNAlandscapeR bucket
#' samples <- c("GM12878", "IMR90")
#' # two <- getSamplesFromBucket(samples, "dnalandscaper", c("chr1", "chr4"))

#' @export
setGeneric(name = "getSamplesFromBucket", def = function(samples, bucket, chr = NA, organism = "human")
    standardGeneric("getSamplesFromBucket"))

#' @rdname getSamplesFromBucket
setMethod("getSamplesFromBucket", def= function(samples, bucket, chr = NA, organism = "human") {
    lo <- lapply(samples, function(sample) getSampleFromBucket(sample, bucket, chr = chr, organism = organism)) 
    names(lo) <- samples
    mdTotal <- data.frame(bucket = rep(bucket, length(samples)))
    row.names(mdTotal) <- samples
    obj <- new("sparseHiCdata", HiCSamplesList = lo, metaData = mdTotal)
})


#' List Available Hi-C Samples from an AWS Bucket
#'
#' \code{samplesInBucket} returns a character vector of the samples
#' that can be imported from the specified bucket. Should also specify
#' the organism.
#'
#' @param bucket name of the AWS bucket (all lower case)
#' @importFrom gsubfn strapplyc
#' @importFrom RCurl getURL
#'
#' @examples
#' # samplesInBucket("dnalandscaper")

#' @export
setGeneric(name = "samplesInBucket", def = function(bucket)
    standardGeneric("samplesInBucket"))

#' @rdname samplesInBucket
setMethod("samplesInBucket", signature("character"),
          definition = function(bucket) {
              t <- gsubfn::strapplyc(RCurl::getURL((paste0("https://s3.amazonaws.com/", bucket))), "<Key>(.*?)</Key>", simplify = c)
              full <- t[grepl("hic", t) & grepl("sparseHiC.meta", t)]
              if(length(full) > 0){
                  t2 <- unlist(strsplit(full, split = "-HiC"))
                  return(sapply(strsplit(t2[grepl("data", t2)], "/"), function(sample){
                      v <- sample[4]
                      names(v) <- sample[2]
                      v
                  }))
              } else {
                  return("none")
              }
          })
