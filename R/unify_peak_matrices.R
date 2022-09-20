#' reduce two peak matrices to a common feature set
#'
#' Filters two peak matrices for their overlapping features and returns them with new feature names,
#' which are identical between the two matrices.
#'
#' @name unify_peak_matrices
#' @param pmat_list list containing two peak matrices, where each matrix has region strings (e.g. chr-110-324) as rows and samples as columns
#' @param minOlap minimum overlap of regions in both directions required to keep them as a common feature
#' @returns list with the two peak matrices, but with identical features as rows
#' @export
#' @examples
#' #get somee random regions as features
#' set.seed(1)
#' gr1 = random_granges(n = 30, chr_annot = structure(5000, names='chr1'), width = sample(1:100,30,replace=T))
#' gr2 = random_granges(n = 20, chr_annot = structure(5000, names='chr1'), width = sample(1:100,30,replace=T))
#' # show the overlap between the two granges
#' granges_overlap(gr1, gr2, return_type = 'message')
#' # intialize peak matrices
#' pmat1 = matrix(1:90, nrow=30, dimnames=list(convert_granges(gr1), paste0('dataset1_cell_', 1:3)))
#' pmat2 = matrix(1:40, nrow=20, dimnames=list(convert_granges(gr2), paste0('dataset2_cell_', 1:2)))
#' # return the matrices with unified features
#' unify_peak_matrices(list(pmat1, pmat2))

unify_peak_matrices = function(pmat_list, minOverlap = 0){
  if(length(pmat_list) != 2){
    stop('Currently only works for unifying two peak matrices, therefore length(pmat_list) must be 2.')
  }
  peaks_gr1 = convert_granges(rownames(pmat_list[[1]]))
  peaks_gr2 = convert_granges(rownames(pmat_list[[2]]))
  
  olap_df = granges_overlap(peaks_gr1, peaks_gr2, minOverlap = minOverlap, olap_direction = 'both')
  #for now just keep the first reported overlap although this is not an optimal solution
  olap_df = olap_df[!(duplicated(olap_df$a_index) | duplicated(olap_df$b_index)), ]
  if(nrow(olap_df) == 0){
    stop('No overlapping regions between the matrices were found.')
  }
  olap_df$new_name = paste0('peak_', 1:nrow(olap_df))
  pmat_list[[1]] = pmat_list[[1]][as.integer(olap_df$a_index), ]
  pmat_list[[2]] = pmat_list[[2]][as.integer(olap_df$b_index), ]
  rownames(pmat_list[[1]]) = rownames(pmat_list[[2]]) = olap_df$new_name
  
  pmat_list
}
