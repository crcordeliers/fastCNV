#' Plot an Annotated Phylogenetic Tree with CNV Events
#'
#' This function generates a plot of an annotated phylogenetic tree using `ggtree`.
#' It displays tip labels, tip points, and labels for CNV events associated with
#' each node.
#'
#' @param tree_data A data frame containing tree structure and annotations,
#' typically produced by `annotateCNVtree`.
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
#' @importFrom ggplot2 aes geom_label scale_x_continuous
#' @export

plotCNVtree <- function(tree_data) {
  tree_plot <- ggtree(tree_data) +
    geom_tiplab(hjust = -0.2) +
    geom_tippoint() +
    geom_label(aes(label = .data$events), size = 3, vjust = -1, hjust = 1) +
    scale_x_continuous(expand = c(0, 0.2)) +
    theme_tree()
  return(tree_plot)
}
