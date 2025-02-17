# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#' @param omim_data A data frame containing OMIM disease-gene associations.
#' @param degs_file A CSV file containing a list of significant DEGs (one column).
#' @return A list containing the network analysis results, including hubs and degree of nodes.
#' @examples
#' construct_comorbid_network(omim_data, "degs.csv")
#' @export
library(igraph)
library(dplyr)
library(circlize)
library(tools)
# Load required libraries
# Load required libraries


# Load required libraries
if(!require(igraph)) install.packages("igraph", dependencies = TRUE)
library(igraph)
library(tools)

# Function to load and process gene expression data from a CSV file
load_gene_expression_data <- function(file_path) {
  data <- read.csv(file_path)
  upregulated <- data$Gene[data$LogFC > 1]
  downregulated <- data$Gene[data$LogFC < -1]
  return(list(upregulated = upregulated, downregulated = downregulated))
}

# Function to calculate shared upregulated and downregulated genes separately
calculate_shared_genes_by_type <- function(central_genes, other_diseases_genes) {
  shared_upregulated <- list()
  shared_downregulated <- list()

  for (disease in names(other_diseases_genes)) {
    shared_up <- intersect(central_genes$upregulated, other_diseases_genes[[disease]]$upregulated)
    shared_down <- intersect(central_genes$downregulated, other_diseases_genes[[disease]]$downregulated)

    shared_upregulated[[disease]] <- shared_up
    shared_downregulated[[disease]] <- shared_down
  }

  return(list(upregulated = shared_upregulated, downregulated = shared_downregulated))
}

# Function to construct the comorbidities network with gene nodes
construct_comorbidities_network_with_genes <- function(central_disease, shared_genes, type = c("upregulated", "downregulated")) {
  type <- match.arg(type)
  network <- data.frame(Node1 = character(), Node2 = character(), stringsAsFactors = FALSE)

  for (disease in names(shared_genes[[type]])) {
    for (gene in shared_genes[[type]][[disease]]) {
      # Create two edges: Central Disease -> Gene, Gene -> Comorbid Disease
      network <- rbind(network, data.frame(Node1 = central_disease, Node2 = gene))  # Central disease to shared gene
      network <- rbind(network, data.frame(Node1 = gene, Node2 = disease))         # Shared gene to comorbid disease
    }
  }

  return(network)
}

# Updated function to visualize the network with shared genes as intermediary nodes
visualize_diseasome_network_with_genes <- function(network, network_type) {
  g <- graph_from_data_frame(network, directed = FALSE)

  # Color the nodes: diseases in green, genes in blue
  V(g)$color <- ifelse(V(g)$name %in% unique(network$Node2), "lightblue", "lightgreen")

  # Plot the graph with gene names inside nodes
  plot(g,
       vertex.size = 30,
       vertex.label.cex = 0.8,
       vertex.label.color = "black",
       edge.color = "gray",
       vertex.shape = "circle",
       main = paste(network_type, "Diseasome Network with Shared Genes"))
}

# Function to export network for Cytoscape
export_for_cytoscape_with_genes <- function(network, file_path) {
  write.table(network, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Main function to load data, construct, analyze, and visualize networks with genes as nodes
ConstructAnalyzeComorbiditiesNetwork <- function(central_file, comorbid_files) {
  # Get central disease name (without .csv)
  central_disease <- file_path_sans_ext(basename(central_file))

  # Load central disease data
  central_genes <- load_gene_expression_data(central_file)

  # Load comorbid/risk factor diseases data
  comorbid_genes <- lapply(comorbid_files, load_gene_expression_data)

  # Assign names to comorbid diseases using the file names (without .csv extension)
  names(comorbid_genes) <- sapply(comorbid_files, function(x) file_path_sans_ext(basename(x)))

  # Find shared upregulated and downregulated genes separately
  shared_genes <- calculate_shared_genes_by_type(central_genes, comorbid_genes)

  # Construct comorbidities network for upregulated genes with genes as nodes
  upregulated_network <- construct_comorbidities_network_with_genes(central_disease, shared_genes, type = "upregulated")

  # Construct comorbidities network for downregulated genes with genes as nodes
  downregulated_network <- construct_comorbidities_network_with_genes(central_disease, shared_genes, type = "downregulated")

  # Visualize the networks with shared genes as intermediary nodes
  visualize_diseasome_network_with_genes(upregulated_network, "Upregulated")
  visualize_diseasome_network_with_genes(downregulated_network, "Downregulated")

  # Export both networks for Cytoscape using the central disease name in the file name
  export_for_cytoscape_with_genes(upregulated_network, paste0(central_disease, "_upregulated_comorbidities_network_with_genes.txt"))
  export_for_cytoscape_with_genes(downregulated_network, paste0(central_disease, "_downregulated_comorbidities_network_with_genes.txt"))

  return(list(upregulated_network = upregulated_network, downregulated_network = downregulated_network))
}

# Example Usage
# Replace these paths with your actual file paths
#central_file <- "path/to/central_disease.csv"
#comorbid_files <- list("path/to/disease1.csv", "path/to/disease2.csv", "path/to/disease3.csv")

# Construct and analyze networks
#results <- ConstructAnalyzeComorbiditiesNetwork(central_file, comorbid_files)

# Example Usage
# Replace these paths with your actual file paths
#central_file <- "central_disease.csv"
#comorbid_files <- list("disease1.csv", "disease2.csv", "disease3.csv")

# Construct and analyze networks
#comorbid_files <- list("SMK.csv", "covid19.csv")
# central_file <- "LC.csv"
#results <- ConstructAnalyzeComorbiditiesNetwork(central_file, comorbid_files)

hello <- function() {
  print("Hello, world!")
}
