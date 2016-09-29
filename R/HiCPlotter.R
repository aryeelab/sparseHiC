#' @include primitive.R
NULL

# Internal function for color mapping
.hicColors <- function(p) {
    if(p == 1) return(colorRampPalette(colorRamps::matlab.like2(100)))
    if(p == 2) return(colorRampPalette(c("#ffffff", "#40826D")))
    if(p == 3) return(colorRampPalette(c("#ffffff", "#9D9A96")))
    if(p == 4) return(colorRampPalette(c("#ffffff", "#2956B2")))
    if(p == 5) return(colorRampPalette(c("#ffffff", "#E34234")))
    if(p == 6) return(colorRampPalette(c("#ffffff", "#E6E6FA")))
    if(p == 7) return(colorRampPalette(c("#ffffff", "#ACE1AF")))
    if(p == 8) return(colorRampPalette(c("#ffffff", "#FF0080")))
    if(p == 9) return(colorRampPalette(c("#ffffff", "#FF9933")))
    if(p == 10) return(colorRampPalette(c("#ffffff", "#E34234")))
    if(p == 11) return(colorRampPalette(c("#ffffff", "#4B0082")))
    if(p == 12) return(colorRampPalette(c("black","blue","#1E90FF","orange","#FF8C00")))
    if(p == 13) return(colorRampPalette(c("black","blue","#1E90FF","#00BFFF","#B0E2FF")))
    if(p == 14) return(colorRampPalette(grDevices::heat.colors(100)))
    if(p == 15) return(colorRampPalette(grDevices::topo.colors(100)))
    if(p == 16) return(colorRampPalette(colorRamps::blue2red(100)))
}

#' Get a color ramp for plotting
#'
#' \code{getHiCcolor} returns a color ramp from one of the 
#' prespecified color options. This can either be accessed
#' from the name or the index (p) of the ramp
#' @param p = Integer index of the plot color them listed below
#' @param name = Character name of the plot color theme listed below
#' @importFrom grDevices colorRampPalette heat.colors topo.colors
#' @importFrom colorRamps matlab.like2 blue2red
#' @import Sushi
#' @examples
#' color.choices <- list(
#'    "Pallet" = 1,
#'    "Viridian" = 2,
#'    "Pewter" = 3,
#'    "Cerulean" = 4,
#'    "Vermillion" = 5,
#'    "Lavendar" = 6,
#'    "Celadon" = 7,
#'    "Fuchsia" = 8,
#'    "Saffron" = 9,
#'    "Cinnabar" = 10,
#'    "Indigo" = 11,
#'    "Master" = 12,
#'    "Black and Blue" = 13,
#'    "Heat" = 14,
#'    "Topology" = 15,
#'    "Blue to Red" = 16)
#'    
#' @export
setGeneric(name = "getHiCcolor", def = function(p = character(0), name = character(0))
    standardGeneric("getHiCcolor"))

#' @rdname getHiCcolor
setMethod("getHiCcolor", def= function(p = character(0), name = character(0)) {
    stopifnot(sum(length(p), length(name)) == 1)
    
    color.choices <- list(
        "Pallet" = 1,
        "Viridian" = 2,
        "Pewter" = 3,
        "Cerulean" = 4,
        "Vermillion" = 5,
        "Lavendar" = 6,
        "Celadon" = 7,
        "Fuchsia" = 8,
        "Saffron" = 9,
        "Cinnabar" = 10,
        "Indigo" = 11,
        "Master" = 12,
        "Black and Blue" = 13,
        "Heat" = 14,
        "Topology" = 15,
        "Blue to Red" = 16)
    
    if(length(name) == 1) p <- color.choices[[name]]
    return(.hicColors(as.integer(p)))
})


#' Plot HiC samples
#'
#' \code{plotHiC} prints a base R image of the region specified for 
#' the sparseHiCdata supplied. 
#'
#' @param x A sparseHiCdata object with samples to be visualized
#' @param y A GRanges object containing region of interest
#' @param res The resolution to be plotted
#' @param libraryNorm = FALSE When plotting multiple samples,
#' normalize for library complexity?
#' @param log2trans = TRUE log2transform the HiC data?
#' @param quantiles = c(5, 95). A vector of length 2 OR 
#' the boolean FALSE. Every value below the first element of the
#' vector (e.g. less than the 5th percentile) will be plotted as the 
#' minimum and everything bigger than the second element of the vector
#' (e.g. greater than the 95th percentile) will be plotted as the maximum.
#' @param color A colorRamp defined by the
#' getHiCcolor function. Default is 13 from the getHiCcolor function
#' @param missingco = "min" By default, missing pixels are plotted
#' with the minimum color, but users can specify other color options
#' so long as they exist in base R (e.g. 'red', 'black', etc.)
#' @param showLegend = TRUE display legend on plots?
#' @param organism = 'h';  'h' for human or 'm' for mouse supported
#' for viewing the gene annotation. Only matters if showGenes is TRUE
#' @param showGenes = FALSE; show gene annotation at bottom panel.
#' 
#' @importFrom graphics plot polygon
#' @importFrom stats quantile
#' @importFrom reshape2 melt
#' @importFrom grDevices recordPlot
#' @importFrom graphics mtext par
#' 
#' @examples
#' 
#' library(GenomicRanges)
#' region <- GRanges(seqnames=c("22"),ranges=IRanges(start=c(0),end=c(9000000)))
#' rdsA<-paste(system.file('rds',package='sparseHiC'),'hESCdata.rds',sep='/')
#' hESCdata <- readRDS(rdsA)
#' region <- GRanges(seqnames=c("chr21"),ranges=IRanges(start=c(23000000),end=c(40000000)))
#' res <- "1000000"
#' a <- plotHiC(hESCdata, region, res, showGenes = FALSE)
#' 
#' @export
setGeneric(name = "plotHiC", def = function(x, y, res, libraryNorm = FALSE, log2trans = TRUE, quantiles = c(5, 95), 
                color = getHiCcolor(16), missingco = "min", showLegend = TRUE, organism = 'h', 
                showGenes = FALSE) standardGeneric("plotHiC"))

#' @rdname plotHiC
setMethod("plotHiC", def = function(x, y, res, libraryNorm = TRUE, log2trans = TRUE, quantiles = c(5, 95), 
                color = getHiCcolor(16), missingco = "min", showLegend = TRUE, organism = 'h', 
                showGenes = FALSE) {
    
    # Data typing; also deal with single sample
    stopifnot(class(x) == "sparseHiCdatum" | class(x) == "sparseHiCdata")
    if(class(x) == "sparseHiCdatum") x <- .as.sparseHiCdata(x)
    
    # Deal with quantiles
    if(length(quantiles) != 2) quantiles <- c(0,0) # no quantile cutoffs
    stopifnot(quantiles[1] <= quantiles[2])
    
    # Handle Region information
    region <- y
    chr <- as.character(seqnames(region))
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    
    sampleNames <- names(x@HiCSamplesList)
    
    # Get the data; check for normalization
    hicdatalist <- lapply(sampleNames, function(sampleName){ getHiCMatrix(x, chr, res, sampleName) })
    if(libraryNorm) hicdatalist <- .libraryNormHiC(hicdatalist)
    names(hicdatalist) <- sampleNames
    
    # FINALLY: Plotting
    hplot <- recordPlot()
    par(mfrow = c(length(sampleNames) + as.numeric(showGenes), 1), mar = c(3, 1, 1, 1), oma = c(0, 0, 3, 0))
    j1 <- sapply(sampleNames, function(sample){
        hicdata <- hicdatalist[[sample]]
        # Call the plot
        .hic.plot(hicdata, region, sample, color = color, log2trans = log2trans,
              missingco = missingco, showlegend = showLegend, Qmin = quantiles[1], Qmax = quantiles[2])
    })
    
    # Add gene annotation if specified
    geneinfo <- ""
    if(showGenes){
        if(organism == "h") {
            rda <- paste(system.file("rda", package = "sparseHiC"), 
                         "geneinfo.h.rda", sep = "/")
        } else if (organism == "m") {
            rda <- paste(system.file("rda", package = "sparseHiC"), 
                         "geneinfo.m.rda", sep = "/")
        }
        load(rda)
        chrom <- gsub("chr", "", chr) 
        geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start - 10000 & geneinfo$stop < end + 10000,]
        
        if(dim(geneinfo)[1] == 0){ #Dummy plot
            plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, 
                      plottype = "loops", heights = 0, lwdrange = c(0, 0), 
                      main = "", adj=0)
        } else {
            pg <- plotGenes(geneinfo = geneinfo, chrom = chrom, chromstart = start, 
                            chromend = end, bheight = 0.1, plotgenetype = "box", 
                            bentline = FALSE, labeloffset = 0.4, fontsize = 1, arrowlength = 0.025, 
                            labeltext = TRUE)
        }
    } 
    mtext(paste0("Region: ", chr, ":", start, "-", end), outer = TRUE, line = 1)
    return(hplot)
})

# Internal method that does the plotting
.hic.plot <- function(hicdata, region, sample, color, log2trans, missingco, showlegend, Qmin, Qmax){
    
    # Set up region
    chromchr <- as.character(seqnames(region))
    start <- as.integer(start(ranges(range(region))))
    end <- as.integer(end(ranges(range(region))))
    palette <- color
    
    rows <- as.numeric(rownames(hicdata))
    cols <- as.numeric(colnames(hicdata))
    
    hicregion <- as.matrix(hicdata[which(rows >= start & rows <= end), which(cols >= start & cols <= end), drop=FALSE])

    if(log2trans) {hicregion <- log2(hicregion); hicregion[is.infinite(hicregion)] <- 0}
    
    if(dim(hicregion)[1]==0 | dim(hicregion)[2]==0){ #Nothing comes up from subsetting
        hicregion <- matrix(0)  
        colnames(hicregion) <- as.character(as.integer(start))
        rownames(hicregion) <- as.character(as.integer(end))
    }

    # determine number of bins
    rvs <- as.numeric(rownames(hicregion))
    cvs <- as.numeric(colnames(hicregion))
    min_bp <-  min(c(rvs, cvs))
    max_bp <-  max(c(rvs, cvs))
    if(length(c(diff(rvs), diff(cvs))) == 0){
        resolution <- 0
    } else {
        resolution <-  min(c(diff(rvs), diff(cvs)))
    }
    if(is.infinite(resolution)){ resolution <- max(rvs,cvs) - min(rvs,cvs)} #1x1 matrix 
    
    if(resolution != 0) {  nbins <- (max_bp-min_bp)/resolution } else { nbins <- 1 }
    
    stepsize <- abs(start - end)/(2 * nbins)
    max_z <- max(hicregion, na.rm = TRUE)
    min_z <- min(hicregion[hicregion != 0], na.rm = TRUE)    
    if(is.infinite(max_z) | is.na(max_z) | is.nan(max_z) | max_z == 0) max_z <- 10000
    if(is.infinite(min_z) | is.na(min_z) | is.nan(min_z)) min_z <- 0.0000001
    
    
    if((as.numeric(Qmax) > as.numeric(Qmin))){
        mreg <- reshape2::melt(hicregion)
        mreg.subset <- mreg[mreg[,3] > 0  & (mreg[,1] != mreg[,2]), ]
        if(dim(mreg.subset)[1] == 0){
            max_z <- as.numeric(hicregion)
            min_z <- as.numeric(hicregion)
        } else {
            max_z <- quantile(mreg.subset[,3], as.numeric(Qmax)*0.01)
            min_z <- quantile(mreg.subset[,3], as.numeric(Qmin)*0.01)
        }
    }    
    
    # map to colors
    breaks <- seq(min_z, max_z, length.out = 100) - 0.001
    cols <- palette(length(unique(breaks)))
    
    if(missingco == "min") { cols <- c(cols[1], cols) } else { cols <- c(missingco, cols) }
    if(length(cols) == 2 ){ cols[2] <- palette(100)[100]}
    
    hicmcol <- matrix(as.character(cut(hicregion, c(-Inf, unique(breaks), Inf), labels = cols)), nrow = nrow(hicregion))
    
    f <- 1; ylim <- c(0, 20); side <- 1

    # initialize plot
    plot(1, 1, xlim = c(start, end), ylim = ylim, type = "n", xaxs = "i", yaxs = "i",
         bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = sample, adj = 0)
    if(dim(hicmcol)[1] != 1) {
    # fill plot
    h <- 20/min(40, dim(hicregion)[2]) * f
    for (rownum in (1:nrow(hicregion))) {
        y = -1*h
        x = start + (rownum * 2 * stepsize) - (stepsize * 3)
        for (colnum in (rownum:ncol(hicregion))) {
            x = x + stepsize
            y = y + h
            if((y <= 20 & f == 1) | (y >= -20 & f == -1)){
                if(colnum != rownum & y!=20*f){ # Square
                    xs = c(x - stepsize, x, x + stepsize, x, x - stepsize)
                    ys = c(y, y + h, y, y - h, y)
                } else if(y == 20 | y == -20){ #upside down triangle; at top
                    xs = c(x - stepsize, x, x + stepsize)
                    ys = c(y, y - h, y)
                } else { #basic triangle
                    xs = c(x - stepsize, x, x + stepsize)
                    ys = c(y, y + h, y)
                }
                if(rownum <= dim(hicmcol)[2] & colnum <= dim(hicmcol)[1]){
                    col <- hicmcol[colnum, rownum]
                } else {col <- cols[1]}
                polygon(xs, ys, border = NA, col = col)
            }
        }
    }
    } else {
        xs = c(start, start+stepsize, end)
        ys = c(0, f*20, 0)
        polygon(xs, ys, border = NA, col = hicmcol[1, 1])
    }
    
    if(min_z == max_z) min_z <- 0
    if(showlegend){
        addlegend(c(min_z, max_z), palette = palette, title="", side="right",
            bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
            width=0.025, title.offset=0.035, labels.digits=1)
    }
    labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
}

