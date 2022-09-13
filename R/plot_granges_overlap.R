#' plot overlaps between two granges objects
#'
#'
#' @name plot_granges_overlap
#' @param a_gr Granges object
#' @param b_gr Granges object to compare with a_gr
#' @param minOverlap minimum required overlap of the region width in a_gr, value between 0-1 (0=minimum 1bp overlap)
#' @param region_gr Granges of length 1 specifying the region that should be plotted
#' @param olap_direction regions in a_gr are highlighted if there is a minOverlap in the specified direction. Currently, regions in b_gr are not highlighted.
#' @export
#' @examples
#' hic_gr <- GRanges(Rle(c("chr1", "chr1", "chrX")),
#'                   IRanges(start = c(3, 22, 1), end=c(5, 25, 10)),
#'                   name= paste0('hic_', 1:3))
#' names(hic_gr) = hic_gr$name
#' atac_gr <- GRanges(Rle("chr1", 5),
#'                    IRanges(start = c(4,5,10,20,25), end=c(30,8,22,24,27)),
#'                    name = paste0('atac_', 1:5))
#' names(atac_gr) = atac_gr$name
#' 
#' #visualize the result (only looking at region on chromosome 1)
#' plot_granges_overlap(a_gr = atac_gr,
#'               b_gr = hic_gr,
#'               region_gr = Signac::StringToGRanges('chr1-1-100'),
#'               minOverlap = 0.2,
#'               olap_direction = 'a')
#' plot_granges_overlap(a_gr = atac_gr,
#'               b_gr = hic_gr,
#'               region_gr = Signac::StringToGRanges('chr1-1-100'),
#'               minOverlap = 0.4,
#'               olap_direction = 'b')
#' plot_granges_overlap(a_gr = unname(atac_gr),
#'               b_gr = unname(hic_gr),
#'               region_gr = Signac::StringToGRanges('chr1-1-100'),
#'               minOverlap = 0.4,
#'               olap_direction = 'both')

plot_granges_overlap = function(a_gr, b_gr, region_gr, minOverlap=0, olap_direction='a'){
  #subset a_gr and b_gr to regions overlapping with `region_gr`
  stopifnot(length(region_gr) == 1)
  a_gr = a_gr[granges_overlap(a_gr, region_gr, minOverlap = 0, return_type = 'logical', olap_direction = 'a', verbose = F)]
  b_gr = b_gr[granges_overlap(b_gr, region_gr, minOverlap = 0, return_type = 'logical', olap_direction = 'a', verbose = F)]
  
  
  a_grolap = granges_overlap(a_gr, b_gr, minOverlap = minOverlap, return_type = 'pairs', olap_direction = olap_direction, verbose=F)
  a_df = as.data.frame(a_gr)
  if(is.null(names(a_gr))){
    a_df$name = paste(seqnames(a_gr), ranges(a_gr), sep= '-')
  }else{
    a_df$name = names(a_gr)
  }
  
  b_df = as.data.frame(b_gr)
  if(is.null(names(b_gr))){
    b_df$name = paste(seqnames(b_gr), ranges(b_gr), sep= '-')
  }else{
    b_df$name = names(b_gr)
  }
  a_df$pos = 1:nrow(a_df)
  a_df$min_overlapping = a_df$name %in% a_grolap$a_gr
  #a_df$min_overlapping = factor(a_df$min_overlapping, levels = c('FALSE','TRUE'))
  b_df$pos = (nrow(a_df)+1):(nrow(b_df)+nrow(a_df))
  #b_df$min_overlapping = b_df$name %in% a_grolap$a_gr

  ggplot() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_segment(data = a_df[a_df$min_overlapping, ], aes(x = start-0.5, xend = end+0.5, y = pos, yend = pos, color='blue'),
                       size = 5) +
    geom_segment(data = a_df[!a_df$min_overlapping, ], aes(x = start-0.5, xend = end+0.5, y = pos, yend = pos, color='lightblue'),
                 size = 5) +
    #to add border color
          #geom_segment(data = a_df, aes(x = start-0.3, xend = end+0.3, y = pos, yend = pos),
          #       size = 4, color='grey') +
          geom_text(data=a_df, aes(x=start-0.5, y=pos, label = name), hjust=1) +
          geom_segment(data = b_df, aes(x = start-0.5, xend = end+0.5, y = pos, yend = pos),
                       size = 5, color='red') +
          geom_text(data=b_df, aes(x=start-0.5, y=pos,  label = name), hjust=1) + 
    scale_colour_manual(name = 'minimum overlap', 
                        values =c('lightblue'='lightblue','blue'='blue'), labels = c('FALSE','TRUE')) + 
    labs(x=seqnames(region_gr), y = '')
}

plot_gene_and_peaks = function(peaks, genes, region){
  #visualize peak and gene positions in one track
  #region: region to plot (may be extended if features intersect borders), eg. chr1-5-10000
  #peaks: strings of peak positions, eg. rownames(seurat_atac)
  #genes can be 'mouse', 'human' or a custom gene genomic ranges object (must contain `symbol` metacolumn)
  if(is.character(genes)){
    if(genes == 'mouse'){
      genes = read_rds('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-09-11-mouse_normal_markers/archr_mm10_gene_annot.RDS')[[1]]
    }else if(genes=='human'){
      genes = read_rds('/omics/groups/OE0436/data2/simonma/projects/scATAC/analysis/2020-09-11-mouse_normal_markers/archr_hg19_gene_annot.RDS')[[1]]
    }
  }
  names(genes) =  paste(seqnames(genes), ranges(genes), sep = '-')
  region_plot = StringToGRanges(region)
  names(region_plot) = region
  peaks = StringToGRanges(peaks)
  names(peaks) = paste(seqnames(peaks), ranges(peaks), sep = '-')
  gene_olap = a_regions_that_have_minOverlap_with_b_regions(genes, region_plot, minOverlap = 0.00000001, return_a_regions_only = T)
  peak_olap = a_regions_that_have_minOverlap_with_b_regions(peaks, region_plot, minOverlap = 0.00000001, return_a_regions_only = T)
  gene_df = data.frame(StringToGRanges(gene_olap))
  peaks_df = data.frame(StringToGRanges(peak_olap))
  gene_df$name = paste(genes$symbol[match(gene_olap, names(genes))], 
                       paste(gene_df$start, gene_df$end, sep=','), 
                       sep = ': ')
  peaks_df$name = paste(peaks_df$start, peaks_df$end, sep = '-')
  gene_df$pos = 1:nrow(gene_df)
  peaks_df$pos = nrow(gene_df)+1
  
  min_use = min(start(region_plot),peaks_df$start, gene_df$start)
  max_use = max(end(region_plot), peaks_df$end, gene_df$end)
  ggplot() +  
    ylim(c(0.5,nrow(gene_df)+1.5)) + 
    xlim(min_use, max_use) +
    labs(x='position', y='genes (blue) / peaks (red)', title = region) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_segment(data = gene_df, aes(x = start, xend = end, y = pos, yend = pos),
                 size = 40, color='blue') +
    geom_segment(data = peaks_df, aes(x = start, xend = end, y = pos, yend = pos),
                 size = 20, color='red') + 
    geom_label(data=gene_df, aes(x=start, y=pos,  label = name), hjust=0) + #, 
    #color='grey', nudge_x = c(-300, -200, -200, -100))
    ggrepel::geom_label_repel(data = peaks_df,
                              aes(x = start, y = pos, label = name),
                              #segment.color = "grey", color = 'black', fill='brightgrey',
                              nudge_x = 0.3, nudge_y=0.4, size = 3, direction = "both")
}
