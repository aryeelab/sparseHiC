#' @include coreHelp.R
NULL

# Internal function that deals with the user-specified chromosome builds
.chrDistBuild <- function(genomeBuild, manual.chr, manual.dist) {
    if(!is.na(genomeBuild)){
        if(genomeBuild == "hg19"){
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_hg19.sizes", sep = "/")
        } else if(genomeBuild == "hg18") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_hg18.sizes", sep = "/")
        } else if(genomeBuild == "hg38") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_hg38.sizes", sep = "/")
        } else if(genomeBuild == "mm9") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_mm9.sizes", sep = "/")
        } else if(genomeBuild == "mm8") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_mm8.sizes", sep = "/")
        } else if(genomeBuild == "mm10") {
            sizeFile <- paste(system.file("extdata", package = "sparseHiC"), 
                              "chrom_mm10.sizes", sep = "/")
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

# Computes pairwise 
.matrixBuild.intra <- function(chr, bed.GRanges, dat.long, res, dist, n){
    options(scipen=999)
    cur.chrom <- bed.GRanges[seqnames(bed.GRanges) == chr]
    vals <-  mcols(cur.chrom)$region 
    dat.chrom <- dat.long[dat.long$idx1 %in% vals & dat.long$idx2 %in% vals, ]
    dat.chrom$region[abs(dat.chrom$idx1 - dat.chrom$idx2) > n] <- 0
    starts <- start(cur.chrom)
    names(starts) <- vals
    dat.chrom$idx1 <- starts[as.character(dat.chrom$idx1)]
    dat.chrom$idx2 <- starts[as.character(dat.chrom$idx2)]
    
    # Create zeroes matrix to handle missing data
    bins <- seq(0, as.numeric(dist[chr]), as.numeric(res))
    zeros.long <- cbind(t(combn(bins, 2)), 0)
    zeros.long <- rbind(zeros.long, cbind(bins, bins, 0))
    colnames(zeros.long) <- c("idx1", "idx2", "region")

    return(reshape2::acast(data = rbind(dat.chrom, zeros.long), formula = idx2 ~ idx1,
                       value.var = "region", fill = 0, fun.aggregate = sum))
    
}

# Computes pairwise (faster and less memory)
.matrixBuild.intra2 <- function(chr, bed.GRanges, dat.long, res, dist, n){
    options(scipen=999)
    cur.chrom <- bed.GRanges[seqnames(bed.GRanges) == chr]
    vals <-  mcols(cur.chrom)$region
    dat.chrom <- dat.long[dat.long$idx1 %in% vals & dat.long$idx2 %in% vals, ]
    if( n > 0 ){
        dat.chrom$region[abs(dat.chrom$idx1 - dat.chrom$idx2) > n] <- 0
    }
    cur.chrom$matIdx <- seq_along(cur.chrom)
    matIdx <- cur.chrom$matIdx
    names( matIdx ) <- cur.chrom$region
    dat.chrom$matIdx1 <- matIdx[as.character( dat.chrom$idx1 )]
    dat.chrom$matIdx2 <- matIdx[as.character( dat.chrom$idx2 )]
    numBins <- length( cur.chrom )
    matr <- matrix( 0, nrow=numBins, ncol=numBins )
    for( j in seq_len( numBins ) ){
        ww <- dat.chrom$matIdx1 == j
        if( sum( ww ) > 0 ){
            matr[j,dat.chrom$matIdx2[ww]] <- dat.chrom$region[ww]
        }
    }
    rownames(matr) <- start( cur.chrom )
    colnames(matr) <- start( cur.chrom )
    matr
}


# Modification of function above to do two chromosomes instead of just one
.matrixBuild.inter <- function(chr1, chr2, bed.GRanges, dat.long, res, dist){
    options(scipen=999)
    cur.chrom1 <- bed.GRanges[seqnames(bed.GRanges) == chr1]
    vals1 <-  mcols(cur.chrom1)$region 
    cur.chrom2 <- bed.GRanges[seqnames(bed.GRanges) == chr2]
    vals2 <-  mcols(cur.chrom2)$region
    dat.chrom <- dat.long[dat.long$idx1 %in% vals1 & dat.long$idx2 %in% vals2, ]
    starts1 <- start(cur.chrom1)
    names(starts1) <- vals1
    starts2 <- start(cur.chrom2)
    names(starts2) <- vals2
    dat.chrom$idx1 <- starts1[as.character(dat.chrom$idx1)]
    dat.chrom$idx2 <- starts2[as.character(dat.chrom$idx2)]
    
    # Create zeroes matrix to handle missing data
    bins1 <- seq(0, as.numeric(dist[chr1]), as.numeric(res))
    bins2 <- seq(0, as.numeric(dist[chr2]), as.numeric(res))
    zeros.long <- cbind(expand.grid(bins1,bins2), 0)
    colnames(zeros.long) <- c("idx1", "idx2", "region")

    return(acast(data = rbind(dat.chrom, zeros.long), formula = idx2 ~ idx1,
                       value.var = "region", fill = 0, fun.aggregate = sum))
}
