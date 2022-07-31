#' generate random GRanges
#'
#' create random granges objects based on supplied chromosome lenghts and desired region widths
#'
#' @name random_granges
#' @param n number of regions to create
#' @param txdb txdb object containing the sizes of chromosomes
#' @param width sizes of the regions, vector of integers from which is sampled
#' @export
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
#' random_granges(n = 30, txdb = txdb, width = pmax(1, as.integer(rnorm(30, mean = 50, sd = 30))))

random_granges <- function(n, txdb, width=1){
  chr_sizes <- seqlengths(txdb)
  chr_sizes = chr_sizes[grepl('^chr[0-9XY]+$', names(chr_sizes))]
  #chr_gr = keepStandardChromosomes(GRanges(seqnames=names(chr_sizes), ranges = IRanges(start = 1, width = 1)), pruning.mode = 'coarse')
  random_chr <- sample(x=names(chr_sizes), size=n, prob=chr_sizes, replace=T)
  random_pos <- sapply(random_chr, function(chrTmp){sample(chr_sizes[names(chr_sizes)==chrTmp],1)})
  random_widhts = sample(width, n, replace = T)
  gr <- GRanges(seqnames = random_chr, ranges = IRanges(start = random_pos, width = random_widhts))
  return(gr)
}