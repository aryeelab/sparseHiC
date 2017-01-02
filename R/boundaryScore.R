#' @include data.R
NULL

# Calculate boundary scores for a matrix (typically one chromosome)
.boundary_score_matrix <- function(x, num_bins=5) {
  x <- as.matrix(x) # TODO: Find another way to get diag()
  # Make sure rows and columns are aligned
  stopifnot(all(rownames(x) == colnames(x)))
  # These index the subsetted matrix
  a_idx <- 1:num_bins
  b_idx <- (num_bins+1):(num_bins*2)
  
  border_idxs <- (num_bins+1) : (ncol(x)-num_bins)
  df <- t(sapply(border_idxs, function(i) {
    idx <- (i-num_bins):(i+num_bins-1)
    xx <- x[idx, idx]
    min_diag <- min(diag(xx))
    diag(xx) <- 0
    a <- sum(xx[a_idx, a_idx])
    b <- sum(xx[b_idx, b_idx])
    c <- sum(xx[b_idx, a_idx])
    c(a, b, c, min_diag)
  }))
  df <- data.frame(pos=as.numeric(rownames(x)[border_idxs]), df)
  colnames(df) <- c("pos", "a", "b", "c", "min_diag")
  df <- subset(df, min_diag>0)
  df$within_vs_cross <- with(df, (a+b) / (c+5))
  return(df)
}

    
#' Compute boundary scores
#'
#' \code{computeBoundaryScores} returns a \code{data.frame} object
#' with boundary scores.
#'
#' @param hic A sparseHiC Object
#' @param sample The name of the sample (e.g. GM12878 or BT142-rep1)
#' @param chrs = paste0("chr", c(1:22))  Specify which chromosomes
#' to compute boundary scores for
#' @param res = "20000"
#' @param num_bins = "human" or "mouse" likely supported
#' @return A data.frame of boundary scores
#' 
#' @import foreach 
#' @examples
#' # Need to make the example better
#' samp <- "IMR90"
#' # hesc <- getSampleFromBucket(samp, "dnalandscaper", c("chr1", "chr4"))
#' 

#' @export
setGeneric(name = "computeBoundaryScores", def = function(hic, sample, chrs, res, num_bins)
    standardGeneric("computeBoundaryScores"))

#' @rdname computeBoundaryScores
setMethod("computeBoundaryScores",
          definition = function(hic, sample, chrs=paste0("chr", c(1:22)), res="20000", num_bins=5) {
  foreach(chr=chrs, .combine="rbind") %dopar% {
    x <- getHiCMatrix(hic, chr=chr, res=res, sampleName = sample)
    data.frame(chr=chr, .boundary_score_matrix(x, num_bins=num_bins))
  }
}
)


