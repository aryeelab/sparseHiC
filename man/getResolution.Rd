% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/access.R
\docType{methods}
\name{getResolution}
\alias{getResolution}
\alias{getResolution,sparseHiCdatum-method}
\title{Grab a sample's specific resolution}
\usage{
getResolution(obj, res, a.list = TRUE)

\S4method{getResolution}{sparseHiCdatum}(obj, res, a.list = TRUE)
}
\arguments{
\item{obj}{A \code{sparseHiCdatum} object}

\item{res}{Resolution(s) desired}

\item{a.list}{= TRUE Return a list of matrices? If FALSE, return another
\code{sparseHiCdatum} object}
}
\value{
Returns object subsetted by resolution
}
\description{
\code{getResolution} takes a \code{sparseHiCdatum} object
and returns a list of sparse Hi-C matrices. If a list is 
returned, then the names are the chromosomes.
}
\examples{
resolutions <- "1000000"
rds<-paste(system.file('rds',package='sparseHiC'),'hESCdatum1.rds',sep='/')
hESCdatum1 <- readRDS(rds)
r <- getResolution(hESCdatum1, resolutions, a.list = FALSE) #trivial
}

