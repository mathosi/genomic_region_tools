#' percentage signature overlap
#' 
#' For each cell calculate how high the percentage of signature peaks covered by the peaks of that cell is.
#' To avoid bias due to count depth, perform an aggregation to above a count threshold among nearest neighbours
#'
#' @name get_percentage_overlap
#' @param peak_matrix binary matrix with peak regions in rows (format chr1-1-99) and cells in columns
#' @param proj ArchR project to use
#' @param reduction dimensionality reduction from ArchR proj, which is used to search nearest neighbours
#' @param reduced_dim_df data.frame with cells in rows and coordinates from dimension reduction in columns
#' @param signature_gr peak signature to quantify as GRanges object
#' @param nFrags_vec vector with the fragment count per cell (same length as ncol in peak matrix)
#' @param count_thres number of fragments to aggregate cells to (nearest neighbour counts are added until the threshold is reached)
#' @param k number of nearest neighbours to search for each cell
#' @param min_overlap minimum overlap of the peak regions with signature regions to count as overlapping feature
#' @param workers number of workers to use
#' @param verbose print progress messages
#' @export
#' @examples
#' 
#' 
# sig_df = get_percentage_overlap(peak_matrix = BinarizeCounts(peak_matrix), 
#                                 reduced_dim_df = Embeddings(SObj, 'pca'),
#                                 signature_gr = peaks_gr, 
#                                 nFrags_vec = SObj$nCount_RNA, 
#                                 count_thres = 50000)


get_percentage_overlap = function(peak_matrix, reduced_dim_df, signature_gr, nFrags_vec, 
                                 count_thres = 3e5, k=100, min_overlap=0.5, workers=1, unify_peak_matrix = TRUE, 
                                 verbose = FALSE){
  library(BiocParallel)
  library(GenomicRanges)

  multicoreParam <- MulticoreParam(workers = workers, progressbar = verbose)
  stopifnot(max(peak_matrix) == 1)
  stopifnot(identical(colnames(peak_matrix), rownames(reduced_dim_df)))
  
  names(signature_gr) <- paste(seqnames(signature_gr), ranges(signature_gr), sep='-')
  coordinates = StringToGRanges(rownames(peak_matrix))
  names(coordinates) <- paste(seqnames(coordinates), ranges(coordinates), sep='-')
  
  stopifnot(class(coordinates)=='GRanges' & class(signature_gr)=='GRanges')
  #return peaks from seurat that overlap with signature peaks
  if(verbose) message('Determining signature - peak matrix overlap.')
  regionsSignature <- granges_overlap(coordinates,
                                      signature_gr, 
                                      minOverlap=min_overlap,
                                      return_type = 'strings')
  if(length(regionsSignature) == 0){
    stop('No overlap between peaks and signature found.')
  }
  if(verbose) message(sprintf('\nSignature contains %i regions, %i regions in the peak matrix overlap with them.',
                  length(signature_gr), length(regionsSignature)))
  
  ##
  if(unify_peak_matrix){
    if(verbose) message('Unifying peak matrix and signature regions')
    dummy_mat = matrix(dimnames = list(names(signature_gr), c('A','B')), 
                       nrow = length(signature_gr), ncol = 2, data = 0)
    peak_matrix = unify_peak_matrices(list(peak_matrix, dummy_mat), type = 'union')[[1]]
    peak_matrix[peak_matrix>1] = 1
    coordinates = convert_granges(rownames(peak_matrix))
    regionsSignature <- granges_overlap(coordinates,
                                        signature_gr, 
                                        minOverlap=min_overlap,
                                        return_type = 'strings')
    if(verbose) message(sprintf('\nSignature contains %i regions, %i regions in the unified peak matrix overlap with them.\nMaximum possible overlap: %.2f%%',
                                length(signature_gr), length(regionsSignature), length(regionsSignature)/length(signature_gr)*100))
  }
  ##
  
  if(verbose) message(sprintf('Finding %i nearest neighbours for each cell.', k))
  nn_map2 <- FNN::get.knnx(data=reduced_dim_df, query=reduced_dim_df, k = k)
  n_ID = t(nn_map2$nn.index)
  
  if(verbose) message('Determining number of neighbours to aggregate for each cell.')
  cum_sum_mat = apply(n_ID, 2, function(x) cumsum(nFrags_vec[x])) #each column is a cell, each row adds another neighbour cell to the counts
  n_ID_cumsum = apply(cum_sum_mat, 2, function(x) which(x > count_thres)[1]) #find the number of cells to aggregate to reach the count threshold
  ncount_after = sapply(1:ncol(cum_sum_mat), function(x) cum_sum_mat[n_ID_cumsum[x], x]) #fragment number after aggregation
  if(anyNA(n_ID_cumsum)){
    stop('NA in ID list. Set k higher or count_thres lower.')
  } #if na, then too high thres or too little neighbours
  #get the column IDs of cells to aggregate for each cell as a list
  ids_use_list = lapply(1:ncol(peak_matrix), function(x) n_ID[1:n_ID_cumsum[x], x])
  
  #message('Subsetting peak matrix.')
  #get subset of counts matrix with peaks overlapping signature peaks
  pmat_use <- peak_matrix[rownames(peak_matrix) %in% regionsSignature, ]
  
  if(verbose) message('Counting number of signature peaks overlapped by each aggregated cell.')
  vals_vec = pmat_use@i
  ids_vec = rep(1:ncol(pmat_use), times = Matrix::colSums(pmat_use))
  # observed_peaks = sapply(1:ncol(pmat_use), function(x){
  #   if(x %% 1000 == 0) message(sprintf('%i/%i cells calculated.', x, ncol(pmat_use)))
  #   olap = length(unique(vals_vec[ids_vec %in% ids_use_list[[x]]]))
  #   return(olap)
  # }) 
  
  observed_peaks = unlist(bplapply(1:ncol(pmat_use), function(x){
    #if(x %% 1000 == 0) message(sprintf('%i/%i cells calculated.', x, ncol(pmat_use)))
    olap = length(unique(vals_vec[ids_vec %in% ids_use_list[[x]]]))
    return(olap)
  }, BPPARAM = multicoreParam))
  
  # pmat_use2 = pmat_use
  # for(i in 1:ncol(pmat_use)){
  #   if(i %% 1000 == 0) message(sprintf('%i/%i cells calculated.', i, ncol(pmat_use)))
  #   pmat_use2[, i] =  pmin(1, Matrix::rowSums(pmat_use[, ids_use_list[[i]], drop=F]))
  # }
  # stopifnot(max(pmat_use2) == 1)
  
  #observed_peaks <- Matrix::colSums(pmat_use2)
  signature_pct_overlap = observed_peaks / length(signature_gr) *100 #percent of full signature that is overlapped
  dataset_max_pct_overlap = observed_peaks / length(regionsSignature) *100  #percent of overlapping peaks in the dataset that have >0 fragments
  sig_df = data.frame(signature_pct_overlap = signature_pct_overlap, 
                      dataset_max_pct_overlap = dataset_max_pct_overlap,
                      signature_n_overlap = observed_peaks,
                      n_cells_aggregated = n_ID_cumsum,
                      nFrags_aggregated =  ncount_after)
  return(sig_df)
}

#' @rdname get_percentage_overlap
#' @export
get_signature_overlap_archr = function(proj, signature_gr, reduction='UMAP', count_thres = 3e5, k=100, min_overlap=0.5){
  dr_df = get_reduction_coords(proj, reduction, embedding_only = T)
  
  matse = getMatrixFromProject(proj, 'PeakMatrix', binarize = T)
  pmat = assays(matse)[[1]]
  pmat = pmat[, match(proj$cellNames, colnames(pmat))]
  stopifnot(identical(colnames(pmat), proj$cellNames))
  
  stopifnot(identical(rownames(dr_df), colnames(pmat)))
  feat_avail = ArchR::getPeakSet(proj)#[sample_use2]
  rownames(pmat) = paste(seqnames(feat_avail), ranges(feat_avail), sep='-')
  
  sig_df = get_percentage_overlap(peak_matrix = pmat, reduced_dim_df = dr_df, signature_gr = signature_gr,
                                 nFrags_vec = proj$nFrags, count_thres = count_thres, k = k, min_overlap = min_overlap)
  return(sig_df)
}

#' @rdname get_percentage_overlap
#' @export
get_signature_overlap_seurat = function(seurat_obj, signature_gr, nFrags_vec,
                                        assay ='scATAC_raw', slot='counts', reduction='umap',
                                        count_thres = 3e5, k=100, min_overlap=0.5){
  #for seurat
  dr_df = Embeddings(seurat_obj, reduction)
  pmat = GetAssayData(seurat_obj, assay = assay, slot=slot)
  sig_df = get_percentage_overlap(peak_matrix = pmat, 
                                 reduced_dim_df = dr_df, 
                                 signature_gr = signature_gr,
                                 nFrags_vec = nFrags_vec, count_thres = count_thres, k = k, min_overlap = min_overlap)
  return(sig_df)
}