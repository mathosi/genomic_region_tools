#' compare the percent overlap between two granges
#'
#' returns the regions in a_gr that have at least minOverlap with the regions in b_gr (other width of overlap comparisons can be achieved by setting `olap_direction`
#' as well as the percentage of a_gr that overlaps b_gr
#'
#' @name granges_overlap
#' @param a_gr Granges object
#' @param b_gr Granges object to compare with a_gr
#' @param minOverlap minimum required overlap of the region width in a_gr, value between 0-1 (0=minimum 1bp overlap)
#' @param olap_direction Direction for which minOverlap is applied. Should regions in a_gr have minOverlap with regions in b_gr (choose 'a'), vice versa ('b'), in both directions ('both') or any direction is fine ('any')
#' @param return_type Should a data.frame with pairs of overlapping regions in a_gr and b_gr ('pairs'), a logical vector indicating regions in a_gr with minOverlap in b_gr ('logical'), a vector of the region strings from a_gr with minOverlap with b_gr ('strings') be returned, or a message summarizing the percentage overlap be returned
#' @returns Depending on `return_type` returns a data.frame, logical or character vector, or just a message
#' @export
#' @examples
#' #example1: generate random granges 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
#' gr1 = random_granges(n = 20000, txdb = txdb, chroms = 'chrY', width = pmax(1, as.integer(rnorm(20000, mean = 300, sd = 30))))
#' gr2 = random_granges(n = 5000, txdb = txdb, chroms = 'chrY', width = pmax(1, as.integer(rnorm(5000, mean = 500, sd = 100))))
#' 
#' #get the overlap percentage summarized
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'a', return_type = 'message', verbose = F)
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'both', return_type = 'message', verbose = F)
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'any', return_type = 'message', verbose = F)
#' #get the matched regions in a_gr and b_gr as a data.frame
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'any')
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'both', verbose = F)
#' #return as logical vector that can be used to subset a_gr
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'both',return_type = 'logical', verbose = F)
#' #or directly return a_gr regions as region strings
#' granges_overlap(a_gr = gr1, b_gr = gr2, minOverlap = 0.4, olap_direction = 'both',return_type = 'strings', verbose = F)
#' 
#' #example2: calculate overlap of hic regions with atac regions
#' hic_gr <- GRanges(Rle(c("chr1", "chr1", "chrX")),
#'                   IRanges(start = c(3, 22, 1), end=c(5, 25, 10)),
#'                   name= paste0('hic_', 1:3))
#' names(hic_gr) = hic_gr$name
#' atac_gr <- GRanges(Rle("chr1", 5),
#'                    IRanges(start = c(4,5,10,20,25), end=c(30,8,22,24,27)),
#'                    name = paste0('atac_', 1:5))
#' names(atac_gr) = atac_gr$name
#' 
#' #3 atac regions have minimum overlap of 20% of their length with a hic region
#' granges_overlap(a_gr = atac_gr,
#'                 b_gr = hic_gr,
#'                 minOverlap = 0.2,
#'                 olap_direction = 'a')
#' 
#' #visualize the result (only looking at region on chromosome 1)
#' plot_granges_overlap(a_gr = atac_gr, 
#'               b_gr = hic_gr, 
#'               region_gr = Signac::StringToGRanges('chr1-1-100'),
#'               minOverlap = 0.2,
#'               olap_direction = 'a')

granges_overlap = function(a_gr, b_gr, minOverlap = 0, olap_direction = 'a',
                           return_type = 'pairs', verbose = TRUE){
  
  #returns the regions in a_gr (column selectedRegions) that have at least minOverlap with the regions in b_gr (shown in the column selected coordinates)
  #not vice versa in contrast to getOverlapRegionsFromCoordinates
  stopifnot(minOverlap >=0 & minOverlap <=1)
  olap_direction = match.arg(arg = olap_direction, choices = c('a','b','both','any')) 
  return_type = match.arg(arg = return_type, choices = c('pairs','logical','strings', 'message'))
  if(is.null(names(a_gr))){
    #if(verbose) warning("a_gr must have names set. Setting names(a_gr) <- paste(seqnames(a_gr), ranges(a_gr), sep='-')")
    names(a_gr) <- paste(seqnames(a_gr), ranges(a_gr), sep='-')
  }
  if(is.null(names(b_gr))){
    #if(verbose) warning("b_gr must have names set. Setting names(b_gr) <- paste(seqnames(b_gr), ranges(b_gr), sep='-')")
    names(b_gr) <- paste(seqnames(b_gr), ranges(b_gr), sep='-')
  }
  a_gr$name = names(a_gr)
  b_gr$name = names(b_gr)
  names(a_gr) = seq_along(a_gr)
  names(b_gr) = seq_along(b_gr)
  
  a_grolap = a_regions_that_have_minOverlap_with_b_regions(a_gr, b_gr, minOverlap = 0, verbose = F)
  colnames(a_grolap) = c('a_index', 'b_index', 'pct_a_olap_b')
  
  b_grolap = a_regions_that_have_minOverlap_with_b_regions(b_gr, a_gr, minOverlap = 0, verbose = F)
  colnames(b_grolap) = c('b_index', 'a_index', 'pct_b_olap_a')
  
  pair_df = suppressMessages(dplyr::full_join(a_grolap, b_grolap))
  pair_df$a_name = a_gr$name[as.integer(pair_df$a_index)]
  pair_df$b_name = b_gr$name[as.integer(pair_df$b_index)]
  pair_df = pair_df %>% dplyr::select(a_index, b_index, a_name, b_name, pct_a_olap_b, pct_b_olap_a)
  
  if(olap_direction == 'a'){
    pair_df = pair_df[pair_df$pct_a_olap_b > minOverlap, ]
  }else if(olap_direction == 'b'){
    pair_df = pair_df[pair_df$pct_b_olap_a > minOverlap, ]
  }else if(olap_direction == 'any'){
    pair_df = pair_df[pair_df$pct_a_olap_b > minOverlap | pair_df$pct_b_olap_a > minOverlap, ]
  }else if(olap_direction == 'both'){
    pair_df = pair_df[pair_df$pct_a_olap_b > minOverlap & pair_df$pct_b_olap_a > minOverlap, ]
  }
  
  if(return_type == 'pairs'){
    return(pair_df)
  }else if(return_type == 'strings'){
    region_matches = unique(pair_df[, 'a_name'])
    return(region_matches) 
  }else if(return_type == 'logical'){
    return(names(a_gr) %in% pair_df$a_index)
  }else if(return_type == 'message'){
    la_gr = length(a_gr) 
    lb_gr = length(b_gr) 
    lolap1 = length(unique(pair_df$a_index)) 
    lolap2 = length(unique(pair_df$b_index))
    lolap1pct = lolap1 / la_gr *100 #fraction of a_gr regions overlapping region in b_gr
    lolap2pct = lolap2 / lb_gr *100 #fraction of b_gr regions overlapping region in a_gr
    message(sprintf('Minimum overlap of %f%% in direction %s\nFraction of a_gr regions: %i/%i = %.3f%%\nFraction of b_gr regions: %i/%i = %.3f%%',
                    minOverlap*100, olap_direction, lolap1, la_gr, lolap1pct, lolap2, lb_gr, lolap2pct))
  }
}


a_regions_that_have_minOverlap_with_b_regions <- function(
    a_gr,
    b_gr,
    minOverlap=0,
    verbose = TRUE,
    ...)
{
  stopifnot(minOverlap >=0 & minOverlap <=1)
  if(is.null(names(a_gr))){
    if(verbose) warning("a_gr must have names set. Setting names(a_gr) <- paste(seqnames(a_gr), ranges(a_gr), sep='-')")
    names(a_gr) <- paste(seqnames(a_gr), ranges(a_gr), sep='-')
  } 
  if(is.null(names(b_gr))){
    if(verbose) warning("b_gr must have names set. Setting names(b_gr) <- paste(seqnames(b_gr), ranges(b_gr), sep='-')")
    names(b_gr) <- paste(seqnames(b_gr), ranges(b_gr), sep='-')
  } 
  
  dbRegionsOverlap <- GenomicRanges::findOverlaps(a_gr, b_gr, type='any', select="all", ignore.strand=TRUE, ...)
  if(length(dbRegionsOverlap)==0){
    if(verbose) warning('No overlaps found.')
    selectedMapped <- as.data.frame(matrix(,nrow = 0, ncol=3))
    colnames(selectedMapped) = c('selectedRegions', 'selectedCoordinates', 'percentOverlapCoordinates')
    return(selectedMapped)
  }
  
  overlaps <- pintersect(a_gr[queryHits(dbRegionsOverlap)], b_gr[subjectHits(dbRegionsOverlap)])
  percentOverlapCoordinates <- width(overlaps) / width(a_gr[queryHits(dbRegionsOverlap)])
  dbRegionsOverlap <- dbRegionsOverlap[percentOverlapCoordinates >= minOverlap]
  if(length(dbRegionsOverlap)==0){
    if(verbose) warning(sprintf('No overlaps with >= %f %% overlap found.', minOverlap*100))
    selectedMapped <- as.data.frame(matrix(,nrow = 0, ncol=3))
    colnames(selectedMapped) = c('selectedRegions', 'selectedCoordinates', 'percentOverlapCoordinates')
    return(selectedMapped)
  }
  percentOverlapCoordinates <- percentOverlapCoordinates[which(percentOverlapCoordinates >= minOverlap)]
  
  selectedRegions <- a_gr[queryHits(dbRegionsOverlap)]
  selectedRegions <- names(selectedRegions) #paste(as.vector(seqnames(selectedRegions)), '-', as.vector(start(selectedRegions)), '-', as.vector(end(selectedRegions)), sep='')
  selectedCoordinates <- names(b_gr[subjectHits(dbRegionsOverlap)])
  selectedMapped <- data.frame(selectedRegions, selectedCoordinates, percentOverlapCoordinates, row.names=NULL)
  
  return(selectedMapped)
}
