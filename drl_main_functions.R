
#' Loads the basic data for the analysis to the global environment
#' @param dc_in the general data containder
#' @param parameters a list with parameters for the analysis
#' @param base_dir the base directory for the analysis
#' @param data_dir directory with where the data is
#' @param act_file a file with activty values
#' @param exp_file a file with expression values
#' @param net_file the file with the network table

load_data <- function(dc_in,
                      parameters = dc_in$parameters,
                      base_dir = parameters$base_dir,
                      data_dir = parameters$data_dir,
                      act_file = parameters$act_file,
                      exp_file = parameters$exp_file,
                      net_file = parameters$net_file){
  #######
  gene_mapping <- readRDS(gene_mapping_file)
  act <- readRDS(file.path(base_dir,data_dir, act_file))
  exp <- readRDS( file.path(base_dir,data_dir, exp_file))
  net <- readRDS( file.path(base_dir,data_dir, net_file)) %>% as_tibble()
 #######adding a comment 
 #assiging main data 
  dc_in$data$input$activity <- act
  dc_in$data$input$expression <- exp
  dc_in$data$input$network_table <- net
  dc_in$data$input$gene_mapping <- gene_mapping
  
  return(dc_in)
  
}

# a comment here can't hurt
preprocess_data <- function(dc_in, 
                            parameters = dc_in$parameters, 
                            net = dc_in$data$input$network_table, 
                            gene_mapping = dc_in$data$input$gene_mapping,
                            collapse_edges = parameters$collapse_edges,
                            filter_edge_type = parameters$filter_edge_type,
                            calculate_alternative_distance = parameters$calculate_alternative_distance,
                            act =  dc_in$data$input$activity ,
                              exp = dc_in$data$input$expression) {
  
  
  tidy_activity <- act %>% exprs() %>% as_tibble(rownames = "gene") %>% separate(gene, into = c("sym", "type"), sep = "_") %>% pivot_longer(cols = -c(sym,type), names_to = "sample", values_to = "activity")
  
  tidy_expression <- act %>% exprs() %>% as_tibble(rownames = "sym")  %>% pivot_longer(cols = -c(sym), names_to = "sample", values_to = "expression")
  
  
  sample_info <- act %>% pData() %>% as_tibble(rownames = "sample")
  
  
  if(collapse_edges){ 
    dist_df <-  net %>%
      group_by(source.symbol, target) %>%
      summarize_at(c("MI","pearson","spearman","rho","p.value"),mean) %>% #grouping and summarizing multiple edges in genes
      ungroup()
  } else if (filter_edge_type){
  
    dist_df <-  net %>%  
      separate(source, into = c("sym","type"), sep = "_") %>% 
      filter(type == "TF") %>% 
      dplyr::select(c("source.symbol" ,"target","MI","pearson","spearman","rho","p.value"))
  
    } else {
    #browser()
    dist_df <- net %>% dplyr::select(c("source.symbol" ,"target", "MI","pearson","spearman","rho","p.value"))
  }
  
  if(calculate_alternative_distance){
    dist_df <- dist_df %>% mutate(transformed_p = -log(p.value + transformation_factor)) %>% 
      mutate(ad_hoc_weight = transformed_p * MI)  
  }
  
  dist_df <- dist_df %>%    dplyr::filter((source.symbol %in% gene_mapping$hgnc_symbol) & (target %in% gene_mapping$hgnc_symbol)) 
  
  genes_in_network <- levels(factor(c( dist_df$target, dist_df$source.symbol)))  #
  gene_mapping <- gene_mapping %>% dplyr::filter(hgnc_symbol %in% genes_in_network )  
  
  dc_in$data$input$gene_mapping <- gene_mapping
  
  dc_in$data$input$dist_df <- dist_df
  dc_in$data$input$tidy_activity <- tidy_activity
  dc_in$data$input$tidy_expression <- tidy_expression
  dc_in$data$input$sample_info <- sample_info
  return(dc_in)
  
}



make_tidygraph <- function(dc_in, 
                           parameters = dc_in$parameters,  
                           dist_df = dc_in$data$input$dist_df,
                           gene_mapping  = dc_in$data$input$gene_mapping,
                           rand_sample = parameters$rand_sample,
                           sample_n = parameters$sample_n) {
  
  
  
  
  links1 <- dist_df %>% dplyr::select(source.symbol, target, MI)  %>% 
    dplyr::rename(from = source.symbol, to = target, weight = MI) #what we are using for the dimensionality; can be MI, ad_hoc_weighht, transformed_p
  
  nodes1 <- gene_mapping %>% group_by(hgnc_symbol) %>% slice_head(n = 1) %>% ungroup() %>% dplyr::rename({{id_col}} := hgnc_symbol ) %>% dplyr::select(.data[[id_col]],everything())  #getting the nodes to make the network
  
  
  
  network <- tbl_graph(edges = links1, nodes = nodes1, directed = TRUE, node_key = id_col)
  
  if(rand_sample){
    
    network  <- network %>% activate(nodes) %>%  slice(sample.int(n = dim(as_tibble(network)), size = sample_n))
  }
  
  dc_in$network <- network
  
  
  
  
  return(dc_in)
}






#dh takes too long with parameeters maxiter 10, fineiter 20
#lgl got bad results with default parameters, timing might be acceptanle
#kk took too long   with maxiter 10, kkconst = 1000
#fr got bad results with default parameters
#got poor results with gem la dee dadd
#stress performed pooly with defaults and took 900 seconds
#graphopt was taking a long time
#backbone was taking a long time, could be worth an overnight test; have to preprocdss network with network <-  network %>% tidygraph_to_igraph() %>% igraph::as.undirected() %>% igraph::simplify() %>% as_tbl_graph()
#eigen can be worth trying
#gem takes too long
layout_network_tidygraph <- function(dc_in,
                                     network = dc_in$network,
                                     parameters = dc_in$parameters,
                                     layout_used =  parameters$layout_used,
                                     options_used = parameters$layout_parameters[[layout_used]]
) {
  
  network <- network %>% enforce_tidygraph() %>% activate(nodes) 
  options_used <- dc_in$parameters$drl_parameters
  
  
  layout <- do.call(
    what = create_layout, 
    args = c(list(graph = network,layout = layout_used  ), options_used)
  )  %>%
    as_tibble()
  
  
  dc_in$network <- dc_in$network %>% left_join(layout[,c(id_col,"x","y")], by = .data[[id_col]])
  dc_out <- dc_in                                       
  
  return(dc_in)
}



cluster_network_igraph <- function(dc_in,network = dc_in$network,
                                   parameters = dc_in$parameters,
                                   clust_alg_used = parameters$cluster_algorithm,
                                   cluster_parameters = parameters$cluster_parameters[[clust_alg_used]]
) {
  network <-   enforce_igraph(network)
  
  dc_in$network <- apply_clust_alg_igraph(network = network,clust_alg = clust_alg_used,parameters = cluster_parameters) %>%
    enforce_tidygraph()
  
  dc_out <- dc_in
  return(dc_out)
  
}


