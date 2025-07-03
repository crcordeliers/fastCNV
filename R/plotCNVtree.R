#' Plot an Annotated Phylogenetic Tree with CNV Events
#'
#' This function generates a plot of an annotated phylogenetic tree using `ggtree`.
#' It displays tip labels, tip points, and labels for CNV events associated with
#' each node.
#'
#' @param tree_data A data frame containing tree structure and annotations,
#' typically produced by `annotateCNVtree`.
#' @param clone_cols a color palette to color the clones. If NULL, points are
#' not colored. If TRUE, clones are colored using default color palette. If a
#' palette is given, clones are colored following the palette, with
#' values passed to [`scale_color_manual`](ggplot::scale_color_manual).
#'
#' @return A `ggplot` object representing the annotated phylogenetic tree.
#'
#' @examples
#' cnv_matrix <- structure(c(0.2, 0.4, 0, 0, 0.1, 0, 0.1, 0.2, 0.2), dim = c(3L,
#' 3L), dimnames = list(c("Clone 1", "Clone 2", "Clone 3"), c("Region 1",
#'                                                           "Region 2", "Region 3")))
#' tree <- CNVtree(cnv_matrix)
#' tree_data <- annotateCNVtree(tree, cnv_matrix)
#' plotCNVtree(tree_data)
#'
#' @importFrom ggtree ggtree geom_tiplab geom_tippoint theme_tree
#' @importFrom ggplot2 aes geom_label scale_x_continuous scale_color_manual guides
#' @export

plotCNVtree <- function(tree_data, clone_cols = NULL) {
  tree_plot <- ggtree(tree_data) +
    geom_tiplab(hjust = -0.2) +
    geom_label(aes(label = .data$events), size = 3, vjust = -1, hjust = 1, na.rm = TRUE) +
    scale_x_continuous(expand = c(0, 0.2)) +
    guides(color = "none") +
    theme_tree()
  if(is.null(clone_cols)) {
    tree_plot <- tree_plot +
      geom_tippoint()
  } else {
    tree_plot <- tree_plot +
      geom_tippoint(aes(color = label), , size = 3)
    if(length(clone_cols) == sum(tree_data$isTip)) {
      tree_plot <- tree_plot +
        scale_color_manual(values = clone_cols)
    }
  }
  return(tree_plot)
}
