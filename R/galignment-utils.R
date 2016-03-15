##
## The following is extracted from the GenomicAlignments package
##
##  Once BioC adds a generic for isFirstSegment, we'll probably just be able to
##  import the generic and methods
.isFirstSegment.matrix <- function(x)
{
    is_paired <- as.logical(x[ , "isPaired"])
    is_first0 <- as.logical(x[ , "isFirstMateRead"])
    is_last0 <- as.logical(x[ , "isSecondMateRead"])
    ## According to SAM Spec, bits 0x40 (isFirstMateRead) and 0x80
    ## (isSecondMateRead) can both be set or unset, even when bit 0x1
    ## (isPaired) is set. However we are not interested in those situations
    ## (which have a special meaning).
    is_paired & is_first0 & (!is_last0)
}

.isFirstSegment.integer <- function(flag)
{
    bitnames <- c("isPaired", "isFirstMateRead", "isSecondMateRead")
    .isFirstSegment.matrix(bamFlagAsBitMatrix(flag, bitnames=bitnames))
}

.isFirstSegment.GAlignments <- function(x)
  .isFirstSegment.integer(mcols(x)$flag)

.isLastSegment.matrix <- function(x)
{
    is_paired <- as.logical(x[ , "isPaired"])
    is_first0 <- as.logical(x[ , "isFirstMateRead"])
    is_last0 <- as.logical(x[ , "isSecondMateRead"])
    ## According to SAM Spec, bits 0x40 (isFirstMateRead) and 0x80
    ## (isSecondMateRead) can both be set or unset, even when bit 0x1
    ## (isPaired) is set. However we are not interested in those situations
    ## (which have a special meaning).
    is_paired & is_last0 & (!is_first0)
}

.isLastSegment.integer <- function(flag)
{
    bitnames <- c("isPaired", "isFirstMateRead", "isSecondMateRead")
    .isLastSegment.matrix(bamFlagAsBitMatrix(flag, bitnames=bitnames))
}

.isLastSegment.GAlignments <- function(x)
    .isLastSegment.integer(mcols(x)$flag)
