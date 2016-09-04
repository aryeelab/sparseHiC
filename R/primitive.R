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

#' @include coreCompress.R
NULL

#' Combining together multiple Hi-C Samples
#' 
#' First, we check if the resolutions are different but the sample
#' is the same. If so, then we append the new resolution, returning a
#' sparseHiCdatum object. If not, then we append it as a new sample.
#' The metaData handling is designed for functionality, not inclusiveness.
#' 
#' @param e1 A sparseHiCdatum object
#' @param e2 A sparseHiCdatum object
#' 
#' @return A sparseHiCdatum/data object depending on if joining multiple samples
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
        
        ss <- append(e1@resolutionNamedList,  e2@resolutionNamedList)     
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
