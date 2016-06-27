#' @include sparseHiC.R
NULL

# Internal function that deals with the user-specified chromosome builds
chrDistBuild <- function(genomeBuild, manual.chr, manual.dist) {
    if(!is.na(genomeBuild)){
        if(genomeBuild == "hg19"){
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_hg19.sizes", sep = "/")
        } else if(genomeBuild == "hg18") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_hg18.sizes", sep = "/")
        } else if(genomeBuild == "mm9") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_mm9.sizes", sep = "/")
        } else { 
            stop("Specified genomeBuild is not supported; try manually inputted chromosome names and distances")
        }
        sdf <- read.table(sizeFile)
        dist <- sdf$V2
        names(dist) <- sdf$V1
        
    } else {
        stopifnot(!is.na(manual.chr), !is.na(manual.dist), length(manual.dist) == length(manual.chr))
        dist <- manual.dist
        names(dist) <- manual.chr
    }
    return(dist)
}

matrixBuild <- function(chr, bed.GRanges, dat.long, res, dist, n){
    options(scipen=999)
    cur.chrom <- bed.GRanges[seqnames(bed.GRanges) == chr]
    vals <-  mcols(cur.chrom)$region 
    dat.chrom <- dat.long[dat.long$idx1 %in% vals & dat.long$idx2 %in% vals, ]
    starts <- start(cur.chrom)
    names(starts) <- vals
    dat.chrom$idx1 <- starts[as.character(dat.chrom$idx1)]
    dat.chrom$idx2 <- starts[as.character(dat.chrom$idx2)]
    
    # Create zeroes matrix to handle missing data
    bins <- seq(0, as.numeric(dist[chr]), as.numeric(res))
    zeros.long <- cbind(t(combn(bins, 2)), 0)
    zeros.long <- rbind(zeros.long, cbind(bins, bins, 0))
    colnames(zeros.long) <- c("idx1", "idx2", "region")

    mat.chrom <- dcast(data = rbind(dat.chrom, zeros.long), formula = idx2 ~ idx1,
                       value.var = "region", fill = 0, fun.aggregate = sum)
    row.names(mat.chrom) <- as.character(format(mat.chrom[, 1], scientific = FALSE))
    mat.chrom <- mat.chrom[, -1]
    
    i <- dim(mat.chrom)[1]
    j <- dim(mat.chrom)[2]
    if (n > 0){for (k in 1:i) mat.chrom[(k + n):j, k] <- 0} 
    mat.chrom <- mat.chrom[1:i, 1:j]
}