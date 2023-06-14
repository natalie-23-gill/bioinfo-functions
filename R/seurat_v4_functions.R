#'
NULL

#' downsample_seurat: downsample the Seurat object to the minimum number of cells per
#' metadata category (sample, batch etc.)
#' 
#' @param seurat_obj A Seurat object
#' @param metadata_slot The name of the metadata slot 
#' @param target_cells The minimum number of cells per sample
#' @return A Seurat object with the minimum number of cells per sample
#' 
#' @export
#' 
downsample_seurat <- function(seurat_obj, metadata_slot, target_cells) {
  # Get the unique sample IDs from the metadata
  ids <- unique(seurat_obj@meta.data[[metadata_slot]])
  cells_to_keep <- c()
  for (sam in ids) {
    # Get the cells for this sample
    cells <- Cells(seurat_obj)[which(seurat_obj[[metadata_slot]] == sam)]
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
  print(table(seurat_obj@meta.data[[metadata_slot]]))
  return(seurat_obj)
}
#' refine_metadata_levels: removes unused metadata levels from a Seurat object,
#' useful for removing unused levels from the metadata after subsetting. Credit to
#' github user michael-kotliar for this function 
#' (https://github.com/satijalab/seurat/issues/5069#issuecomment-1372599846).
#' 
#' @param seurat_obj A Seurat object
#' @return A Seurat object with unused metadata levels removed
#' 
#' @export
#' 
refine_metadata_levels <- function(seurat_data){
    for (i in base::colnames(seurat_data@meta.data)){
        if (base::is.factor(seurat_data@meta.data[[i]])){
            base::print(base::paste("Re-evaluating levels for a factor column", i))
            base::print(
                base::paste(
                    "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
            seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
            base::print(
                base::paste(
                    "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
        }
    }
    return (seurat_data)
}