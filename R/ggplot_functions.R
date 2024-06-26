#'
NULL

#' multipage_plot: Creates and saves a multipage PDF of a list of ggplot
#' objects, with a specified number of plots per page. Plots are
#' automatically sized to fit the page.
#' 
#' @param plot_list a list of ggplot objects
#' @param per_page the number of plots per page
#' @param filename the name of the output file
#' @param ncol the number of columns per page
#' @param page_width the width of the page (default = 8.5 inches)
#' @param page_height the height of the page (default = 11 inches)
#' @return a multipage PDF of the plots written to the filename specified
#' 
#' @export
#' 
multipage_plot <- function(plot_list,
                           per_page,
                           filename,
                           ncol = 2,
                           page_width = 8.5,
                           page_height = 11) {
  # Create a PDF file to save the plots
  pdf(file = filename, height = page_height, width = page_width)
  
  # Split the plots into groups of per_page per page
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / per_page)
  plot_indices <- split(seq_len(num_plots),
                        rep(seq_len(num_pages),
                            each = per_page)) |>
    suppressWarnings()
  
  # Plot the plots on each page
  for (i in seq_len(num_pages)) {
    if (length(plot_indices[[i]]) < per_page) {
      # Fill in the rest of the last page with blank plots
      this_index <- max(plot_indices[[i]])
      for (j in seq_len(per_page - length(plot_indices[[i]]))) {
        this_index <- this_index + 1
        plot_indices[[this_index]] <- ggplot() +
          geom_blank()
        plot_indices[[i]] <- c(plot_indices[[i]], this_index)
      }
    }
    plots <- ggpubr::ggarrange(
      plotlist = plot_list[plot_indices[[i]]],
      ncol = ncol, nrow = ceiling(length(plot_indices[[i]]) / ncol)
    )
    plots_grob <- ggplotGrob(plots)
    grid::grid.newpage()
    grid::grid.draw(plots_grob)
  }
  # Close the PDF file
  dev.off()
}
