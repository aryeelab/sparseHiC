
context("Basic Hi-C Pro impot works")

test_that("We get a list of 23 chromosomes and the matrix is square for chromosome 1", {
    matrix.file <- paste(system.file("extdata", package = "processedHiCdata"),
                          "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
    bed.file <- paste(system.file("extdata", package = "processedHiCdata"),
                      "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
    genomeBuild <- "hg19"
    res <- "1000000"
    sampleName <- "hESC1"
    x <- import.HiCPro(matrix.file, bed.file, res, sampleName, genomeBuild)
    expect_equal(length(x@resolutionNamedList), 1)
    expect_equal(length(x@resolutionNamedList[[1]]), 23)
    expect_equal(dim(x@resolutionNamedList[[1]][[1]])[1], dim(x@resolutionNamedList[[1]][[1]])[2])
    expect_equal(dim(x@resolutionNamedList[[1]][[1]])[1], 250)
})

