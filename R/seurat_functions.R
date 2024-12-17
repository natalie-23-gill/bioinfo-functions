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
  ids <- as.character(unique(seurat_obj@meta.data[[metadata_slot]]))
  
  # Define a helper function to downsample cells for a given sample
  downsample_cells <- function(sam) {
    # Get the cells for this sample
    cells <- Cells(seurat_obj)[which(seurat_obj[[metadata_slot]] == sam)]
    if (length(cells) > target_cells) {
      # If there are more cells than the target, downsample
      keep_cells <- sample(cells, target_cells)
    } else {
      # If there are fewer cells than the target, keep all
      keep_cells <- cells
    }
    return(keep_cells)
  }
  
  # Use lapply to iterate over sample IDs and downsample cells
  cells_to_keep <- unlist(lapply(ids, downsample_cells))
  
  # Subset the seurat_obj using the cells to keep
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  message("Downsampled Seurat object to minimum number of cells per sample, cells per sample:")
  message(table(seurat_obj@meta.data[[metadata_slot]]))
  
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
            base::message(base::paste("Re-evaluating levels for a factor column", i))
            base::message(
                base::paste(
                    "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
            seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
            base::message(
                base::paste(
                    "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
        }
    }
    return (seurat_data)
}

#' plot_metadata: Plots all of the categorical metadata in grouped
#' and split UMAPs, autosizes based on the number of unique values
#' in the metadata. Plots are output as PDFs to a directory. Does not
#' work with coord_fixed() so not suitable for final publication plots.
#'
#' @param seurat_obj A Seurat object
#' @param output_dir Output directory
#' @param pt.size Point size for the UMAP plots (default is 1)
#' @param plot_only A character vector of metadata columns to plot
#' @return Writes out plots to directory
#'
#' @export
#'
plot_metadata <- function(seurat_obj,
                          output_dir,
                          pt.size = 1,
                          plot_only = NULL) {
  ## Plots all of the categorical metadata in grouped and split UMAPs,
  ## autosizes based on the number of unique values in the metadata
  if (!(dir.exists(output_dir))) {
    dir.create(output_dir, recursive = T)
  }
  # get the metadata columns that are not numeric
  meta_cols <- names(seurat_obj@meta.data)[!sapply(seurat_obj@meta.data, is.numeric)]
  # Exclude "orig.ident"
  meta_cols <- meta_cols[meta_cols != "orig.ident"]

  # If plot_only is specified, limit the meta_cols to those specified
  if (!is.null(plot_only)) {
    meta_cols <- intersect(meta_cols, plot_only)
  }
  for (meta in meta_cols) {
    message(paste0("***** Plotting ", meta, " *****"))
    pdf(file.path(paste0(output_dir, "/", meta, "_UMAP.pdf")),
      width = 10,
      height = 10
    )
    print(DimPlot(seurat_obj,
      raster = FALSE,
      order = TRUE,
      label = FALSE,
      group.by = meta,
      pt.size = pt.size
    )) 
    dev.off()
    h <- ifelse(length(unique(seurat_obj@meta.data[[meta]])) == 2,
      10,
      3 * ceiling(length(unique(seurat_obj@meta.data[[meta]])) / 2)
    )
    pdf(file.path(paste0(output_dir, "/", meta, "_split_UMAP.pdf")),
      width = 15,
      height = h
    )
    print(DimPlot(seurat_obj,
      raster = FALSE,
      order = TRUE,
      label = FALSE,
      group.by = meta,
      split.by = meta,
      ncol = 2,
      pt.size = pt.size # Set point size here
    ) + NoLegend())
    dev.off()
  }
}

#' plot_metadata_numeric: Plots all of the numeric metadata
#' using feature plot UMAPs. Plots are output as PDFs to a 
#' directory.
#' 
#' @param seurat_obj A Seurat object
#' @param output_dir Output directory
#' @return Writes out plots to directory
#' 
#' @export
#' 
plot_metadata_numeric <- function(seurat_obj, output_dir) {
  # Plot UMAPs for all numeric metadata columns
  nonchar_cols <- names(seurat_obj@meta.data)[sapply(seurat_obj@meta.data, is.numeric)]
  for (meta in nonchar_cols) {
    message(paste0("***** Plotting ", meta, " *****"))
    pdf(file.path(paste0(output_dir, "/", meta, "_UMAP.pdf")),
      width = 10,
      height = 10
    )
    print(FeaturePlot(seurat_obj,
      features = meta,
      raster = FALSE,
      order = TRUE
    ))
    dev.off()
  }
}

#' Convert a Seurat object to BPCells on-disk format
#'
#' Converts the counts matrix of specified assays in a Seurat object to BPCells
#' format, saving the matrices on disk and updating the Seurat object to use
#' these on-disk matrices. Only works for single count layers.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_dir Directory where the BPCells matrices will be saved. Defaults
#' to a subdirectory in the system's temporary directory named after the Seurat object.
#' @param assays Character vector of assays to convert. Defaults to "RNA".
#' @return The updated Seurat object using on-disk matrices.
#' @export
#' @examples
#' \dontrun{
#' # Example usage
#' seurat_obj <- readRDS("/path/to/your/seurat_object.rds")
#' seurat_obj <- convert_seurat_to_bpcells(seurat_obj)
#' }
convert_seurat_to_bpcells <- function(seurat_obj, output_dir = NULL,
                                      assays = "RNA") {
  # Derive the name of the seurat object
  obj_name <- deparse(substitute(seurat_obj))
  
  # Set default output directory to TMPDIR using the object's name
  if (is.null(output_dir)) {
    output_dir <- file.path(tempdir(), obj_name)
  }
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Iterate over each specified assay in the Seurat object
  for (assay_name in assays) {
    if (!assay_name %in% names(seurat_obj@assays)) {
      warning(paste(
        "Assay", assay_name,
        "not found in the Seurat object. Skipping."
      ))
      next
    }
    
    # Convert to v5 assay
    seurat_obj[[assay_name]] <- as(object = seurat_obj[[assay_name]], Class = "Assay5")
    # Check if the counts matrix is already in BPCells format to avoid reprocessing
    if (inherits(seurat_obj[[assay_name]]@layers$counts, "BPMatrix")) {
      message(paste(
        "Counts matrix for assay",
        assay_name, "is already in BPCells format. Skipping."
      ))
      next
    }
    
    # Write counts matrix to BPCells format
    counts_dir <- file.path(output_dir, paste0(assay_name, "_counts"))
    BPCells::write_matrix_dir(mat = seurat_obj[[assay_name]]@layers$counts, dir = counts_dir)
    
    # Update the counts matrix to on-disk BPCells matrix
    seurat_obj[[assay_name]]@layers$counts<- BPCells::open_matrix_dir(dir = counts_dir)
  }
  
  # Return the updated Seurat object
  return(seurat_obj)
}