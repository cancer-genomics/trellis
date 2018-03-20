#' Compute allele frequencies at germline heterozygous positions
#
#' Given a set of genomic coordinates this function identifies heterozygous positions in a bam file 
#' given by \code{normalBam} and retuns the minor allele frequency at these positions in a bam file 
#' specified by \code{tumorBam}.  
#'  
#' @param normalBam The path to a bam file 
#' @param tumorBam The path to a bam file
#' @param genome The genome build of \code{tumorBam} and \code{normalBam}.  Possible values are "hg38", "hg19", and "hg18". 
#' @param positions A \code{GRanges} object consisting of the genomic regions of interest.  The default set of positions is 
#' the \code{snps} object from \code{svfilters.hg38}, \code{svfilters.hg19} or \code{svilters.hg18} depending on the user-specified \code{genome} argument.  The \code{snps}
#' object contains 1,000,000 positions that are frequently seen as heterozygous, therefore fewer positions are needed to find a sufficient 
#' number of heterozygous positions.  The \code{snps} object contains more than enough positions for tumor ploidy/ploidy analysis on well-covered
#' WGS bam files, although for bam files genertated from targeted capture-based sequencing data it is recommended to use all or most the \code{dbsnp150_snps}
#' object in \code{svfilters.hg38}, \code{svfilters.hg19} or \code{svfilters.hg18} to find a sufficient number of heterozygous positions as these objects contain
#' over 12,000,000 positions.  
#' @param region If \code{region} is specified only SNPs in \code{position} that fall in \code{region} will 
#' be used.  This argument is useful for capture-based sequencing technologies (e.g. Whole-exome sequencing) where we tend to only have high enough
#' coverage to accurately call SNPs in the targeted regions.  By supplying a \code{GRanges} object containing the targeted 
#' regions this function avoids calculating coverage metrics and allele frequencies for low coverage off-target regions. Alternatively, 
#' users can provide the full path to a bed file in BED format.
#' @param n The number of positions to use for pileup in \code{normalBam}.  The \code{snps} objects 
#' in \code{svfilters.hg38}, \code{svfilters.hg19} and \code{svfilters.hg18} contain 1 million frequently heterozygous positions.  By specifying the \code{n} argument
#' a random sample of the positions object of length \code{n} will be used.  For tumor purity/ploidy
#' analysis 10,000 heterozygous positions well spread across the genome is typically plenty.  The 
#' default value of \code{n = 50000} is generally sufficient to achieve this on a 30X WGS bam file.  For sequencing data
#' from targeted capture protocols, it is recommended to use the full set of SNPs in \code{dbsnp150_snp} in \code{svfilters.hg38}, 
#' \code{svfilters.hg19} or \code{svfilters.hg18} to identify a suffient number of SNPs with high enough coverage.
#' @param minCovNormal The minimum coverage of a position in \code{normalBam} to be considered.
#' @param minCovTumor The minimum coverage of a position in \code{tumorBam} to be considered.
#' @param max_depth The maximum number of alignments considered for each position.
#' @param minMafNormal The minimum minor allele frequency (MAF) in order to consider a position 
#' as heterozygous in \code{normalBam}.
#' @param minMafTumor The minimum minor allele frequency (MAF) in order to consider a position 
#' as heterozygous in \code{tumorBam}.  It is recommended to set this value to at least 0.05 
#' when setting \code{normalBam = NULL} to avoid outputting homozygous positions.  
#' 
#' @details If using this function to generate allele frequencies in a tumor 
#' sample at germline heterozygous positions identified in a matched normal sample then
#' \code{tumorBam} should point to the bam file for the tumor sample and \code{normalBam} 
#' should point to the bam file for its matched normal.
#' 
#' @examples 
#' extdir <- system.file("extdata", package="svbams")
#' bam1 <- file.path(extdir, "cgov10t.bam")
#' bam2 <- file.path(extdir, "cgov44t_revised.bam")
#' 
#' data(snps, package = "svfilters.hg19")
#' snps <- keepSeqlevels(snps, c("chr3", "chr5"), pruning.mode = "coarse")
#' \dontrun{
#' svAF(normalBam=bam1, 
#'      tumorBam=bam2, 
#'      genome="hg19", 
#'      positions = snps, 
#'      n = 1000, 
#'      minCovNormal = 10, 
#'      minCovTumor = 10, 
#'      max_depth = 1e3, 
#'      minMafNormal = 0.3, 
#'      minMafTumor = 0)
#'}
#' @return
#' A data.frame with the following columns:
#' 
#' \strong{Chrom}: The chromosome of the event \cr
#' \strong{Pos}: The coordinate of the event \cr
#' \strong{RefBase}: The base in the reference genome (corresponds to the refUCSC column in dbSNP build 150) \cr
#' \strong{AltBase}: The other observed base \cr
#' \strong{Normal.Mut.Count}: The coverage of AltBase in normalBam \cr
#' \strong{Normal.Coverage}: The distinct coverage of RefBase + AltBase in normalBam \cr
#' \strong{Tumor.Mut.Count}: The coverage of AltBase in tumorBam \cr
#' \strong{Tumor.Coverage}: The distinct coverage of RefBase + AltBase in tumorBam \cr
#' \strong{Tumor.MAF}: The minor allele frequency of the event in tumorBam
#' 
#' @export
#' 
svAF <- function(normalBam, 
                 tumorBam, 
                 genome, 
                 positions,
                 region,
                 n = 50000, 
                 minCovNormal = 20, 
                 minCovTumor = 20, 
                 max_depth = 1e3, 
                 minMafNormal = 0.3,
                 minMafTumor = 0) {
  
  if (!is.null(normalBam)) {
    if (!file.exists(normalBam)) {
      stop(paste0(normalBam), " doesn't exist!")
    }
  }
  
  if (!is.null(tumorBam)) {
    if (!file.exists(tumorBam)) {
      stop(paste0(tumorBam), " doesn't exist!")
    }
  }
  
  if (minCovNormal < 0) {
    stop("minCovNormal cannot take on a negative value.")
  }
  
  if (minCovTumor < 0) {
    stop("minCovTumor cannot take on a negative value.")
  }
  
  if ((minMafNormal < 0) | (minMafNormal > 0.5)) {
    stop("minMafNormal must be in the range [0,0.5]")
  }
  
  if ((minMafTumor < 0) | (minMafTumor > 0.5)) {
    stop("minMafTumor must be in the range [0,0.5]")
  }
  
  if (missing(positions)) {
    data(snps, package = paste0("svfilters.", genome), envir = environment())
    SNPs <- snps
  } else {
    SNPs <- positions
  }
  
  if (max_depth <= 0) {
    stop("The value of 'max_depth' must be greater than 0")
  }
  
  if (n < 1) {
    stop("The value of 'n' must be at least 1")
  }
  
  if (n > length(SNPs)) {
    stop("The value of 'n' must be less than or equal to the number of SNPs in 'positions'")
  }
  
  if (!(genome %in% c("hg18", "hg19", "hg38"))) {
    stop(paste0(genome, " is not a possible value for 'genome'.  Possible values include 'hg19' and 'hg18'"))
  }
  
  if (!missing(region)) {
    if (!is(region, "GRanges")) {
      if (file.exists(region)) {
        region <- bed2gr(bedPath = region)
      } else {
        stop("The file path given as an argument to 'region' does not exist.  Make sure that the full path to
             the bed file is specified e.g. '/Users/XSVuser/Documents/capture-region.bed'")
      }
      }
    }
  
  out.df <- data.frame(Chrom = character(0),
                       Pos = character(0),
                       RefBase = character(0),
                       AltBase = character(0),
                       Normal.Mut.Count = character(0),
                       Normal.Coverage = character(0),
                       Tumor.Mut.Count = character(0),
                       Tumor.Coverage = character(0),
                       Tumor.MAF = character(0))
  
  if (!missing(region)) {
    SNPs <- subsetByOverlaps(SNPs, region)
    if (length(SNPs) > 0) {
      message(paste0(length(SNPs), " SNPs in 'positions' overlap with 'region'"))
    } else {
      stop(paste0(length(SNPs), " SNPs in 'positions' overlap with 'region'"))
    }
  }
  
  if (n > length(SNPs)) {
    n <- length(SNPs)
    message("Scaling down 'n' to ", n)
  }
  
  SNPs <- SNPs[sample(1:length(SNPs), size = n, replace = FALSE)]
  
  if (!is.null(normalBam) & !is.null(tumorBam)) {
    message(paste0("Computing pileup at ", n, " position(s) in ", normalBam))
    querySNPs <- SNPs
    normalPU <- pileup(normalBam,
                       scanBamParam=ScanBamParam(which=querySNPs,
                                                 flag = scanBamFlag(isDuplicate = FALSE)),
                       pileupParam=PileupParam(max_depth=max_depth,
                                               distinguish_strands=FALSE,
                                               include_deletions=FALSE))
    
    if (nrow(normalPU) == 0) {
      warning(paste0("0 coverage was found at every specified position in", normalBam), call. = FALSE)
      return(out.df)
    }
    
    message("Identifying heterozygous positions")
    normalSNPs <- filterSNPs(pu = normalPU, SNPs = querySNPs, min.cov = minCovNormal, min.maf = minMafNormal, keepSingles = FALSE)
    
    if (length(normalSNPs) == 0) {
      warning(paste0("0 heterozygous positions were found in", normalBam), call. = FALSE)
      return(out.df)
    }
    
    message(paste0("Computing pileup at ", length(normalSNPs), " position(s) in ", tumorBam))
    tumorPU <- pileup(tumorBam,
                      scanBamParam=ScanBamParam(which=normalSNPs,
                                                flag = scanBamFlag(isDuplicate = FALSE)),
                      pileupParam=PileupParam(max_depth=max_depth,
                                              distinguish_strands=FALSE,
                                              include_deletions=FALSE))
    
    if (nrow(tumorPU) == 0) {
      warning(paste0("0 coverage was found at every identified germline heterozygous position in", tumorBam), call. = FALSE)
      return(out.df)
    }
    
    tumorSNPs <- filterSNPs(pu = tumorPU, SNPs = querySNPs, min.cov = minCovTumor, min.maf = minMafTumor, keepSingles = TRUE)
    alleleFreqs <- calcAlleleFreq(tumorSNPs = tumorSNPs, normalSNPs = normalSNPs)
    message(paste0("Retuning allele frequencies for ", nrow(alleleFreqs), " position(s)"))
    
    return(alleleFreqs)
  }
  
  # When tumorBam is specified and normalBam isn't, do pileup in tumorBam with minMafTumor
  if (is.null(normalBam) & !is.null(tumorBam)) {
    querySNPs <- SNPs
    message(paste0("Computing pileup at ", length(querySNPs), " position(s) in ", tumorBam))
    tumorPU <- pileup(tumorBam,
                      scanBamParam=ScanBamParam(which=querySNPs,
                                                flag = scanBamFlag(isDuplicate = FALSE)),
                      pileupParam=PileupParam(max_depth=max_depth,
                                              distinguish_strands=FALSE,
                                              include_deletions=FALSE))
    
    if (nrow(tumorPU) == 0) {
      warning(paste0("0 coverage was found at every identified germline heterozygous position in", tumorBam), call. = FALSE)
      return(out.df)
    }
    
    message("Identifying heterozygous positions")
    tumorSNPs <- filterSNPs(pu = tumorPU, SNPs = querySNPs, min.cov = minCovTumor, min.maf = minMafTumor, keepSingles = FALSE)
    
    if (length(tumorSNPs) == 0) {
      warning(paste0("0 heterozygous positions were found in", tumorBam), call. = FALSE)
      return(out.df)
    }
    
    alleleFreqs <- calcAlleleFreq(tumorSNPs = tumorSNPs, normalSNPs = NULL)
    message(paste0("Retuning allele frequencies for ", nrow(alleleFreqs), " position(s)"))
    
    return(alleleFreqs)
  }
  
  }



# Takes as input a pileup output, returns a
# GRanges of heterozygous positions if called with keepSingles = FALSE.
# If called with keepSingles = TRUE, then keeps homozygous positions
# as long as they match the SNP database.  This is to detect high
# purity LOH (e.g. LOH in a cell line).  The high purity LOH
# will be missed in a tumor-only analysis because we can't
# distingiush it from being germline homozygous.
filterSNPs <- function(pu, SNPs, min.cov, min.maf, keepSingles) {
  
  pu$position <- paste(pu$seqnames, pu$pos, sep = ":")
  
  # If a position has more than 2 bases counted (becaue of errors, mutations, ...) then
  # only keep the 2 most frequent alleles
  potentialHet <- as.data.frame(
    pu %>% 
      group_by(position) %>% 
      top_n(2, count)
  )
  
  # If keepSingles == FALSE then single-allele positions are removed.  This will be the
  # case when identifying hets in a matchedNormal, or in a tumor-only analysis to remove homozygotes.  
  # These positions are kept in the tumor in T/N analysis because of LOH in a high purity sample
  if (keepSingles == FALSE) {
    potentialHet <- as.data.frame(
      potentialHet %>%
        dplyr::group_by(position) %>%
        dplyr::filter(n() > 1)
    )
  }
  
  if (nrow(potentialHet) == 0) {
    return(GRanges())
  }
  
  # Remove position with coverage < min.cov
  # Remove positions with alt allele frequency < min.maf
  summed <- potentialHet %>%
    group_by(position) %>%
    summarize(tot = sum(count)) %>%
    filter(tot < min.cov)
  
  mafs <- potentialHet %>%
    group_by(position) %>%
    mutate(maf = count / sum(count)) %>%
    filter(maf < min.maf)
  
  rem <- c(summed$position, mafs$position)
  potentialHet <- subset(potentialHet, !(position %in% rem))
  
  if (nrow(potentialHet) == 0) {
    return(GRanges())
  }
  
  ids <- as.character(unique(potentialHet$position))
  keepers <- GRanges()
  for (i in ids) {
    pos <- potentialHet[which(potentialHet$position == i),]
    
    if (nrow(pos) == 2) {
      genotype <- paste(pos$nucleotide, collapse = "/")
    } else if (nrow(pos) == 1) {
      genotype <- as.character(pos$nucleotide)
    }
    
    chrom <- strsplit(i, split = ":")[[1]][1]
    start <- as.integer(strsplit(i, split = ":")[[1]][2])
    end <- start
    coord <- GRanges(seqnames = chrom,
                     ranges = IRanges(start = start,
                                      end = end))
    snp <- subsetByOverlaps(SNPs, coord)
    
    # Only include positions that match the snp database for both reference and alt alleles if 2 alleles are present
    # If one allele is present, only that allele must match the SNP database
    if (nrow(pos) == 2) {
      if ((length(grep(snp$refUCSC, genotype)) == 0) || (length(grep(gsub("/", "|", snp$altAllele), genotype)) == 0)) {
        next
      }
    } else if (nrow(pos == 1)) {
      if ((length(grep(snp$refUCSC, genotype)) == 0) && (length(grep(snp$altAllele, genotype)) == 0)) {
        next
      }
    }
    
    if (length(which(pos$nucleotide == snp$refUCSC)) == 1) {
      snp$refCount <- pos$count[which(pos$nucleotide == snp$refUCSC)]
    } else {
      snp$refCount <- 0
    }
    
    altAllele <- strsplit(snp$altAllele, split = "/")[[1]]
    ind <- which(pos$nucleotide %in% altAllele)
    if (length(ind) == 1) {
      snp$altCount <- pos$count[ind]
      snp$altAllele <- pos$nucleotide[ind]
    } else if (length(ind) == 0) {
      snp$altCount <- 0
    }
    
    keepers <- c(keepers, snp)
  }
  return(keepers)
}



calcAlleleFreq <- function(tumorSNPs, normalSNPs = NULL, min.cov) {
  message("Computing allele frequencies")
  
  if (!is.null(normalSNPs)) {
    # Matching the indicies of tumorSNPs and normalSNPs 
    normalSNPs <- subsetByOverlaps(normalSNPs, tumorSNPs)
    normalSNPs <- sort(normalSNPs)
  }
  
  tumorSNPs <- sort(tumorSNPs)
  
  maf <- tumorSNPs$refCount / (tumorSNPs$refCount + tumorSNPs$altCount)
  maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]
  
  if (!is.null(normalSNPs)) {
    out.df <- data.frame(Chrom = as.character(seqnames(tumorSNPs)),
                         Pos = start(ranges(tumorSNPs)),
                         RefBase = tumorSNPs$refUCSC,
                         AltBase = tumorSNPs$altAllele,
                         Normal.Mut.Count = normalSNPs$altCount,
                         Normal.Coverage = (normalSNPs$refCount + normalSNPs$altCount),
                         Tumor.Mut.Count = tumorSNPs$altCount,
                         Tumor.Coverage = (tumorSNPs$altCount + tumorSNPs$refCount),
                         Tumor.MAF = round(maf, digits = 2))  
  }
  
  if (is.null(normalSNPs)) {
    out.df <- data.frame(Chrom = as.character(seqnames(tumorSNPs)),
                         Pos = start(ranges(tumorSNPs)),
                         RefBase = tumorSNPs$refUCSC,
                         AltBase = tumorSNPs$altAllele,
                         Normal.Mut.Count = NA,
                         Normal.Coverage = NA,
                         Tumor.Mut.Count = tumorSNPs$altCount,
                         Tumor.Coverage = (tumorSNPs$altCount + tumorSNPs$refCount),
                         Tumor.MAF = round(maf, digits = 2))
  }
  return(out.df)
}


# Input path to bed file, returns GRanges of regions
# It is expected that bedPath points to a tab-delimited text
# file in 3-column bed format.  If more than 3 columns are present, 
# it is expected that:
# Column 1: Chromosome
# Column 2: Start Coordinate
# Column 3 : End Coordinate
bed2gr <- function(bedPath) {
  bed <- read.table(bedPath, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  bed.gr <- NULL
  try({
    bed.gr <- GRanges(seqnames = bed$V1, 
                      ranges = IRanges(start = bed$V2, 
                                       end = bed$V3))
  }, silent = TRUE)
  
  if (is.null(bed.gr)) {
    stop("The supplied bed file cannot be coerced to a GRanges object.  
          Check that the bed file follows BED formatting rules (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
          Note that only the first three BED fields are required.")
  } else {
    return(bed.gr) 
  }
}
