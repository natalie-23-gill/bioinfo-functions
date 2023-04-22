#'
NULL

#' downsample_seurat: downsample the Seurat object to the minimum number of cells per sample
#'
#' @param seurat_obj A Seurat object
#' @param sample_metadata_slot The name of the metadata slot containing the sample IDs
#' @param target_cells The minimum number of cells per sample
#' @return A Seurat object with the minimum number of cells per sample
#' 
#' @export
#' 
downsample_seurat <- function(seurat_obj, sample_metadata_slot, target_cells) {
  # Get the unique sample IDs from the metadata
  sample_ids <- unique(seurat_obj@meta.data[[sample_metadata_slot]])
  cells_to_keep <- c()
  for (sam in sample_ids) {
    # Get the cells for this sample
    cells <- Cells(seurat_obj)[which(seurat_obj[[sample_metadata_slot]] == sam)]
    if (length(cells) > target_cells) {
      # If there are more cells than the target, downsample
      keep_cells <- sample(cells, target_cells)
      cells_to_keep <- c(cells_to_keep, keep_cells)
    } else {
      # If there are fewer cells than the target, keep all
      cells_to_keep <- c(cells_to_keep, cells)
    }
  }
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  print("Downsampled Seurat object to minimum number of cells per sample, cells per sample:")
  print(table(seurat_obj@meta.data[[sample_metadata_slot]]))
  return(seurat_obj)
}
