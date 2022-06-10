#adding a line to the control panel!
library(igraph)
library(Biobase)
library(umap)
library("ggplot2")
library("purrr")
library("tibble")
library("dplyr")
library("tidyr")
library("stringr")
library("readr")
library("forcats")
library("magrittr")
library(reticulate)
library(tidygraph)
library(markdown)
library(dbscan)
library(class)
library(tictoc)
library(ggraph)
library(R.utils)
library(plotly)
library(clusterProfiler)
library(enrichplot)
#base_dir <-  "~/Documents/R/netbid_app/netbid_analysis_pipeline_light/"
data_dir <- "DATA"
transformation_factor = 0.00000000001
color_dir <- "color_palettes"
categorical_palette <- "cat_blackzero.RDS"
#"aeb_colors.rds"
#"out_pal_33.RDS"
categorical_color_file <- file.path(data_dir,color_dir,categorical_palette)


act_file <- "BRCA_TCGA_Lee_act.RDS"
exp_file <-  "BRCA_TCGA_Lee_exp.RDS"
net_file <-  "BRCA_TCGA_Lee_net.RDS"
drl_options <- drl_defaults$default
#default, coarsen, coarsest, refine, final
#making this change on experimental
layout_parameters <- list(
  drl =                   list(
    use.seed = FALSE,
  #  seed = matrix(runif(vcount(graph) * 2), ncol = 2),
    options = drl_defaults$default,
    #weights = E(graph)$weight,
    fixed = NULL, 
    dim = 2
  )
)

rand_sample <- FALSE 
sample_n <- 5000
id_col = "id"

layout_used <- "drl"
pathway_files <- c("rp.csv","mrp.csv")
gene_mapping_file <- "gene_mapping_sig_tf_20200624_ensembl_changed.RDS"
pathway_file <- "rp.csv"   #what file we are coloring by based on ancestor/target
pathway <- read_csv(file.path(".","pathway_lists",pathway_file),col_names = FALSE)
collapse_edges <- TRUE
filter_edge_type <- FALSE
calculate_alternative_distance <- FALSE


pathway_file_1 <- "mrp.csv"   #what file we are coloring by based on ancestor/target
pathway_1 <- read_csv(file.path(".","pathway_lists",pathway_file_1),col_names = FALSE)


pathway_file_2 <- "rp.csv"   #what file we are coloring by based on ancestor/target
pathway_2 <- read_csv(file.path(".","pathway_lists",pathway_file_2),col_names = FALSE)
cluster_algorithm <- "hdbscan_knn"

cluster_parameters <- list(
  edge_betweennes = list(
    directed = TRUE, edge.betweenness = TRUE, merges = TRUE,
    bridges = TRUE, modularity = TRUE, membership = TRUE
  ),
  fast_greedy  = list(
    merges = TRUE, modularity = TRUE,
    membership = TRUE #, weights = E(graph)$weight
  ), 
  infomap =list(
    e.weights = NULL, v.weights = NULL,
    nb.trials = 10, modularity = TRUE  
  ), 
  label_prop = list( weights = NULL, initial = NULL,
                     fixed = NULL),
  leading_eigen = list(steps = -1, weights = NULL,
                       start = NULL, options = arpack_defaults, callback = NULL,
                       extra = NULL#, env = parent.frame()
  ), 
  louvain = list(
    weights = NULL
  ),
  optimal = list(
    weights = NULL
  ),
  
  spinglass = list(
    weights = NULL, vertex = NULL, spins = 25,
    parupdate = FALSE, start.temp = 1, stop.temp = 0.01,
    cool.fact = 0.99, update.rule = c("config", "random", "simple"),
    gamma = 1, implementation = c("orig", "neg"), gamma.minus = 1
  ),
  walktrap = list(
    steps = 4,
    merges = TRUE, modularity = TRUE, membership = TRUE  #, weights = E(graph)$weight
  ),
  hdbscan_knn = list(
    coordinate_variables = c("x","y"),
    id_variable = "id",
    jitter_factor = NULL,
    jitter_amount = 1,
    jitter_cluster = TRUE,
    hdbscan_min_pts = 150,
    force_hard_cluster = FALSE,
    k_for_knn = 50 
  )
)

parameters_used <- as.list(.GlobalEnv)
