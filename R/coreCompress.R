#' @include sparseHiC-class.R
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

# Internal function for library normalization
# Takes a list of sparse matricies and returns the same list

.libraryNormHiC <- function(losm){
    options(scipen=999)
    
    # Set up long matrix to get the differences
    all <- summary(Reduce("+", losm))
    counts <- sapply(losm, function(m){ as.matrix(m[cbind(all$i, all$j)])})
    long <- data.matrix(cbind(as.numeric(colnames(losm[[1]])[all$i]), as.numeric(colnames(losm[[1]])[all$j]), counts))
    diff <- abs(long[,1] - long[,2])
    longdiff <- cbind(long, diff)
    colnames(longdiff) <- c("idx1", "idx2", paste0("s", seq(1, length(losm), 1)), "diff")

    # Infer resolution and max size
    res <- min(diff[diff > 0])
    ma <- res * (dim(losm[[1]])[1]-1)
    
    # Aggregate and append means
    m <- aggregate(x = longdiff, by = list(longdiff[,6]), FUN = "mean")
    scaled <- m[,4:(3 + length(losm))]/rowMeans(m[,4:(3 + length(losm))])
    mlo <- data.frame(cbind(diff = m[,4 + length(losm)], scaled))
    ldm <- merge(longdiff, mlo, by.x = c("diff"), by.y = c("diff"))
    
    # Make long zeros matrix
    bins <- seq(0, ma, res)
    zeros.long <- cbind(t(combn(bins, 2)), 0)
    zeros.long <- rbind(zeros.long, cbind(bins, bins, 0))
    colnames(zeros.long) <- c("idx1", "idx2", "val")
    
    # Apply transform and make new list
    dat <- lapply(1:(length(losm)), function(i){
        a <- data.frame(cbind(idx1=ldm[,2], idx2=ldm[,3], val=ldm[,i+3]/ldm[,i+3+length(losm)]))
        Matrix(reshape2::acast(rbind(a, zeros.long), formula = idx2 ~ idx1, value.var = "val", fill = 0, fun.aggregate = sum))
    })
    return(dat)
}