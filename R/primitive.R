#' @include coreCompress.R
NULL

# Simply a resolution named list slot
# such that a minimal data hierarchy is present
.cleanResolutionNamedList <- function(rnl){
    datflat <- unlist(rnl)
    ndf <- names(datflat)
    
    # Remove duplicates
    d <- duplicated(ndf)
    ndf <- ndf[!d]
    datflat <- datflat[!d]
    
    res <- unique(sapply(strsplit(ndf, split = "\\."), "[[", 1))
    lo  <- lapply(res, function(r){
        d <- datflat[grepl(paste0(r, "."), ndf)]
        names <- sapply(strsplit(ndf[grepl(paste0(r, "."), ndf)], split = "\\."), "[[", 2)
        names(d) <- names
        d
    })
    names(lo) <- res
    return(lo)
}

#' Combining together multiple Hi-C Samples
#' 
#' First, we check if the resolutions are different but the sample
#' is the same. If so, then we append the new resolution, returning a
#' sparseHiCdatum object. If not, then we append it as a new sample.
#' The metaData handling is designed for functionality, not inclusiveness.
#' 
#' @param e1 A sparseHiCdatum object
#' @param e2 A sparseHiCdatum object
#' @param x A sparseHiCdatum object
#' @param y A sparseHiCdatum object
#' 
#' @return A sparseHiCdatum/data object depending on if joining multiple samples
#' 
#' @importFrom BiocGenerics combine
#' 
#' @examples 
#' rdsA<-paste(system.file('rds',package='sparseHiC'),'hESCdatum1.rds',sep='/')
#' hESCdatum1 <- readRDS(rdsA)
#' rdsB<-paste(system.file('rds',package='sparseHiC'),'hESCdatum1.rds',sep='/')
#' hESCdatum2 <- readRDS(rdsB)
#' rdsC<-paste(system.file('rds',package='sparseHiC'),'IMR90datum1.rds',sep='/')
#' IMR90datum1 <- readRDS(rdsC)
#' hESCdata <- hESCdatum1 + hESCdatum2
#' threeA <- IMR90datum1 + hESCdata
#' threeB <- hESCdata + IMR90datum1
#' 
setMethod("+", signature(e1 = "sparseHiCdatum", e2 = "sparseHiCdatum"),
          definition = function(e1, e2) {
              
    #Basic Error Handling
    if((e1@sampleName == e2@sampleName) &&
       (names(e1@resolutionNamedList) == names(e2@resolutionNamedList)) &&
       (names(e1@resolutionNamedList[[1]]) == names(e2@resolutionNamedList[[1]]))){
            stop("No distinct sample names, resolutions, or chromosomes being combined together")
    }
              
    if(e1@sampleName == e2@sampleName){ # Adding resolutions to the same sample
        
        ss <- .cleanResolutionNamedList(append(e1@resolutionNamedList,  e2@resolutionNamedList))  
        md <- e1@metaData
        obj <- new("sparseHiCdatum", sampleName = e1@sampleName, resolutionNamedList = ss, metaData = md)
        
    } else { # Linking multiple samples together
        
        # Collect Samples
        sampleNames <- c(e1@sampleName, e2@sampleName)
        sampleList <- list(e1,e2)
        names(sampleList) <- sampleNames
        
        #Collect Meta Data
        cc <- intersect(colnames(e1@metaData), colnames(e2@metaData))
        mdTotal <- rbind(subset(e1@metaData, select = cc), subset(e2@metaData, select = cc))
        rownames(mdTotal) <- sampleNames
        
        obj <- new("sparseHiCdata", HiCSamplesList = sampleList, metaData = mdTotal)
    }
    return(obj)
})

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("+", signature(e1 = "sparseHiCdata", e2 = "sparseHiCdata"), definition = function(e1, e2) {
    ev1 <- unlist(lapply(e1@HiCSamplesList, function(sample){unlist(sample@resolutionNamedList)}))
    ev2 <- unlist(lapply(e2@HiCSamplesList, function(sample){unlist(sample@resolutionNamedList)}))
    
    datflat <- unlist(append(ev1,ev2))
    ndf <- names(datflat)
    
    # Remove duplicates
    d <- duplicated(ndf)
    ndf <- ndf[!d]
    datflat <- datflat[!d]
    
    res <- unique(sapply(strsplit(ndf, split = "\\."), "[[", 2))
    samples <- unique(sapply(strsplit(ndf, split = "\\."), "[[", 1))
    lo  <- lapply(samples, function(s){
        li <- lapply(res, function(r){
            d <- datflat[grepl(paste0(s, ".", r, "."), ndf)]
            names <- sapply(strsplit(ndf[grepl(paste0(s, ".", r, "."), ndf)], split = "\\."), "[[", 3)
            names(d) <- names
            d
        })
        names(li) <- res
        li
    })
    names(lo) <- samples
    
    #Collect Meta Data
    cc <- intersect(colnames(e1@metaData), colnames(e2@metaData))
    mdTotal <- rbind(subset(e1@metaData, select = cc), subset(e2@metaData, select = cc))
    rownames(mdTotal) <- samples

    obj <- new("sparseHiCdata", HiCSamplesList = lo, metaData = mdTotal)
    return(obj)
})

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("+", signature(e1 = "sparseHiCdatum", e2 = "sparseHiCdata"),
          definition = function(e1, e2) {
              return(.as.sparseHiCdata(e1) + e2)
          })

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("+", signature(e1 = "sparseHiCdata", e2 = "sparseHiCdatum"),
          definition = function(e1, e2) {
              return(e1 + .as.sparseHiCdata(e2))
          })

# Alias function call `combine`

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("combine", signature(x = "sparseHiCdatum", y = "sparseHiCdatum"),
          definition = function(x, y) {
              return(x + y)
          })

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("combine", signature(x = "sparseHiCdata", y = "sparseHiCdata"),
          definition = function(x, y) {
              return(x + y)
          })

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("combine", signature(x = "sparseHiCdatum", y = "sparseHiCdata"),
          definition = function(x, y) {
              return(.as.sparseHiCdata(x) + y)
          })

#' @rdname plus-sparseHiCdatum-sparseHiCdatum-method
setMethod("combine", signature(x = "sparseHiCdata", y = "sparseHiCdatum"),
          definition = function(x, y) {
              return(x + .as.sparseHiCdata(y))
          })
