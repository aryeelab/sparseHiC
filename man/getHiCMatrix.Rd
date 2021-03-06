% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/access.R
\docType{methods}
\name{getHiCMatrix}
\alias{getHiCMatrix}
\alias{getHiCMatrix,sparseHiCdata-method}
\alias{getHiCMatrix,sparseHiCdatum-method}
\title{Grab a specific sparse Hi-C matrix}
\usage{
getHiCMatrix(obj, chr, res, sampleName)

\S4method{getHiCMatrix}{sparseHiCdatum}(obj, chr, res, sampleName)

\S4method{getHiCMatrix}{sparseHiCdata}(obj, chr, res, sampleName)
}
\arguments{
\item{obj}{A \code{sparseHiCdatum} object}

\item{chr}{Chromosome desired}

\item{res}{Resolution desired}

\item{sampleName}{When accessing a matrix with multiple
samples in the object}
}
\value{
Returns a Hi-C Matrix
}
\description{
\code{getHiCMatrix} takes a \code{sparseHiCdatum} object
and returns a single matrix. Thus, the length of chr
and res must be 1
}
\examples{
chr <- "chr21"
res <- "1000000"

#Single sample
rds<-paste(system.file('rds',package='sparseHiC'),'hESCdatum1.rds',sep='/')
hESCdatum1 <- readRDS(rds)
r <- getHiCMatrix(hESCdatum1,chr = chr, res = res) 

rdsA<-paste(system.file('rds',package='sparseHiC'),'hESCdata.rds',sep='/')
hESCdata <- readRDS(rdsA)
getHiCMatrix(hESCdata, chr = chr, res = res, sampleName = "hESC2") 
}

