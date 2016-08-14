
context("Basic Hi-C Pro compression works")

test_that("We get a list of 23 chromosomes and the matrix is square for chromosome 1", {
    matrix.file <- paste(system.file("extdata", package = "processedHiCdata"),
                          "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_iced.matrix", sep = "/")
    bed.file <- paste(system.file("extdata", package = "processedHiCdata"),
                      "HiC-Pro/hESC_Rep1/hESC_Rep1_1000000_abs.bed", sep = "/")
    genomeBuild <- "hg19"
    res <- "1000000"
    x <- sparseCompress.HiCPro(matrix.file, bed.file, res, genomeBuild)
    expect_equal(length(x), 23)
    expect_equal(dim(x[[1]])[1], dim(x[[1]])[2])
    expect_equal(dim(x[[1]])[1], 250)
})

