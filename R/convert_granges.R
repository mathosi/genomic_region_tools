#' GRanges object conversion
#'
#' convert GRanges to strings or data.frame/strings to GRanges
#'
#' @name convert_granges
#' @param obj number of regions to create
#' @param ... arguments passed to `makeGRangesFromDataFrame`
#' @returns if obj is character vector or data.frame, it is converted to a GRanges object. If it is a GRanges object, it is converted to a character vector
#' @export
#' @examples
#' region_vec = c("chr1-105442114-105442778", "chr1-105453942-105454629", "chr1-106231358-106231675")
#' #convert region strings to GRanges
#' (region_gr = convert_granges(region_vec))
#' #convert GRanges to region strings
#' convert_granges(region_gr)
#' #convert data.frame to GRanges
#' region_df = data.frame(seqnames=c('chr1', 'chr2', 'chr2'),
#'                        start=c(1,100,345),
#'                        end = c(306, 3409,909),
#'                        strand = c('*','+', '-'), 
#'                        symbol=c('gene_1', 'gene_2', 'gene_3'))
#' convert_granges(region_df, keep.extra.columns = TRUE)

convert_granges = function(obj, ...){
  library(GenomicRanges)
  if(any(class(obj) == 'character')){
    obj = gsub(':', '-', obj)
    res = do.call('rbind', strsplit(obj, '-'))
    res = data.frame(seqnames = res[, 1], start = as.integer(res[, 2]), end = as.integer(res[, 3]))
    res = makeGRangesFromDataFrame(res, ...)
  }else if(any(class(obj) == 'GRanges')){
    res = paste(seqnames(obj), ranges(obj), sep = '-')
  }else if(any(class(obj) %in% c('data.frame', 'DataFrame'))){
    res = makeGRangesFromDataFrame(obj, ...)
  }else{
    stop('obj should be one of character vector, data.frame or Granges')
  }
  res
}
