#' signature enrichment/depletetion p-value
#'
#' Given a peak signature, (binary) scATAC peakmatrix and a cell-to-cluster annotation, calculate an enrichment or
#' depletion p-value for each cluster using bias-matched background peak sets
#'
#' @name get_signature_pvalue
#' @param pmat binary peak matrix (rows = regions in the format 'chr1-10-100', columns = cells)
#' @param group_vector vector of group annotation for the cells in the peak matrix (e.g. clustering result)
#' @param signature_peaks_gr GRanges object of the peak signature
#' @param genome genome for selection of background sequences, e.g. BSgenome.Mmusculus.UCSC.mm10
#' @param n_background_sets number of background sets to compare against, higher values give more precise estimates
#' @keywords ggplot,plot_grobs,plot_grid
#' @export
#' @examples
#' #bad example since there are no real clusters in this dataset, but the function still works
#' signac_atac = Signac::atac_small
#' peakmat = BinarizeCounts(GetAssayData(signac_atac))
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' sig_gr = Signac::StringToGRanges(rownames(FindMarkers(signac_atac, assay = 'peaks', ident.1 = '0')))
#' res = get_signature_pvalue(pmat = peakmat,
#'                      group_vector = signac_atac$seurat_clusters, 
#'                      signature_peaks_gr = sig_gr,
#'                      genome = BSgenome.Hsapiens.UCSC.hg38)

get_signature_pvalue  = function(pmat, group_vector, signature_peaks_gr, genome, n_background_sets = 300){
  stopifnot(identical(ncol(pmat), length(group_vector)))
  library(SummarizedExperiment)
  library(chromVAR)
  library(jj)
  se = SummarizedExperiment(assays = SimpleList(counts = pmat), 
                            rowRanges = StringToGRanges(rownames(pmat)),
                            colData = DataFrame(group = as.character(group_vector)))
  #library(BSgenome.Mmusculus.UCSC.mm10)
  se <- addGCBias(se, genome = genome)
  se = se[!is.na(rowData(se)$bias), ]
  #draw background peaks from the same peakset for each peak
  message('Getting ', n_background_sets, ' sets of background peaks')
  bgpeaks <- getBackgroundPeaks(se,
                                bias = rowData(se)$bias,
                                niterations = n_background_sets, 
                                w = 0.1, bs = 50)
  message('Calculating peak matrix overlap with signature')
  signature_olap = granges_overlap(rowRanges(se), signature_peaks_gr, return_type = 'logical')
  message(sprintf('%i/%i peaks are overlapping with the signature.',sum(signature_olap), nrow(pmat)))
  bgpeaks = bgpeaks[signature_olap, ]
  
  message('Calculating signature peak set scores per group')
  mean_cluster_acc = jj_summarize_sparse_mat(pmat[signature_olap, ], 
                                             summarize_by_vec = se$group,
                                             method = 'mean')
  #returns matrix with features in rows and groups in columns
  #mean_cluster_acc = scale(t(mean_cluster_acc))
  #mean_cluster_acc = apply(mean_cluster_acc, 1, rank)
  #mean_cluster_acc = rowMeans(mean_cluster_acc) #returns groups in rows and features in columns
  mean_cluster_acc = colMeans(mean_cluster_acc)
  print(sort(mean_cluster_acc, decreasing=T))
  
  message('Calculating background peak set scores per group')
  bg_mat = jj_initialize_df(ncol = length(unique(se$group)), 
                            nrow = n_background_sets, 
                            init = NA, 
                            col.names = names(mean_cluster_acc),
                            return_matrix = T)
  
  for(i in 1:n_background_sets){
    if(i %% 25 == 0) message(i,'/',n_background_sets)
    peaks_use = bgpeaks[, i]
    mean_bg_cluster_acc = jj_summarize_sparse_mat(pmat[peaks_use, ], 
                                                  summarize_by_vec = se$group,
                                                  method = 'mean')
    #mean_bg_cluster_acc = scale(t(mean_bg_cluster_acc))
    #mean_bg_cluster_acc = apply(mean_bg_cluster_acc, 1, rank)
    bg_mat[i, ] = colMeans(mean_bg_cluster_acc) #rowMeans(mean_bg_cluster_acc)
  }
  stopifnot(identical(colnames(bg_mat), names(mean_cluster_acc)))
  enr_score = sapply(1:ncol(bg_mat), function(x) mean(mean_cluster_acc[x] / bg_mat[, x], na.rm = T))
  names(enr_score) = names(mean_cluster_acc)
  enr_score = sort(enr_score, decreasing = T)
  
  message('Computing p values')
  #enrichment
  pvals_bigger = sapply(seq_along(unique(group_vector)), function(x) sum(bg_mat[, x] > mean_cluster_acc[x]) / n_background_sets)
  names(pvals_bigger) = names(mean_cluster_acc)
  padjusted_bigger = sort(p.adjust(pvals_bigger, method = 'BH'))
  
  #depletion
  pvals_smaller = sapply(seq_along(unique(group_vector)), function(x) sum(bg_mat[, x] < mean_cluster_acc[x]) / n_background_sets)
  names(pvals_smaller) = names(mean_cluster_acc)
  
  padjusted = p.adjust(c(pvals_bigger, pvals_smaller), method = 'BH')
  padj_enriched = sort(padjusted[1:length(pvals_bigger)])
  padj_depleted = sort(padjusted[(length(pvals_bigger)+1):length(padjusted)])
  
  # result_list = list(
  #   #p_enriched = pvals_bigger,
  #   #p_depleted = pvals_smaller,
  #   padj_enriched = padj_enriched,
  #   padj_depleted = padj_depleted,
  #   rank_score = avfc #rank of signature / rank of background
  # )
  
  pval_df = data.frame(group = gtools::mixedsort(names(padj_enriched)))
  pval_df$padj_enriched = padj_enriched[match(pval_df$group, names(padj_enriched))]
  pval_df$padj_depleted = padj_depleted[match(pval_df$group, names(padj_depleted))]
  pval_df$enr_score = enr_score[match(pval_df$group, names(enr_score))]
  
  return(pval_df)
}
