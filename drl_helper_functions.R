make_plot <- function(in_dr){
  layout <- in_dr %>% data.frame() %>% as_tibble() %>% ggplot(aes(x = X1, y = X2)) + geom_point(size = .1) 
  return(layout)
}



cluster_hdbscan_knn <- function(graph,
                                coordinate_variables = c("x","y"),
                                id_variable = "id",
                                jitter_factor = NULL,
                                jitter_amount = 1,
                                jitter_cluster = TRUE,
                                hdbscan_min_pts = 150,
                                force_hard_cluster = TRUE,
                                k_for_knn = 5){ 
  graph <- enforce_igraph(graph)
#  browser()
  wanted_attributes <-  as.list(coordinate_variables) %>% setNames(coordinate_variables)
  
  
  vertex_names <- set_names(
    
    
    list(names(V(graph))),
    id_variable
  )
  coordinates <-  pmap(.f = get.vertex.attribute, .l = list(name = wanted_attributes), graph = graph) 
  to_clust <- c( vertex_names, coordinates  ) %>% as_tibble()
  
  
  
  if(jitter_cluster){ 
    
    to_clust <- to_clust %>% dplyr::mutate(across(where(is.numeric),jitter,amount = jitter_amount, factor = jitter_factor))# map_df(jitter,amount = jitter_amount, factor = jitter_factor)
  } 
  
  
  #actually running hdbscan
  
  clusters <-  factor(hdbscan(select(to_clust, where(is.numeric)),minPts = hdbscan_min_pts)$cluster)
  
  
  to_clust$Cluster <- clusters
  
  if(force_hard_cluster){
    # use nearest_neighbors to force any "cluster 0" to join the nearest cluster 
    
    if(  any(to_clust$Cluster == 0)){
      
      
      
      to_clust <- to_clust %>% as_tibble()
      
      
      
      #separating successful clusters from what was labeled as noise points
      
      successful_cluster <- to_clust %>% dplyr::filter(Cluster != 0)
      cluster_zero <- to_clust %>% dplyr::filter(Cluster == 0)
      
      train <- successful_cluster %>% select_if(is.numeric) %>% as.matrix()
      rownames(train) <-  successful_cluster$id
      
      test <-  cluster_zero %>% select_if(is.numeric)    %>% as.matrix()
      rownames(test) <- cluster_zero$id
     
      cl <- successful_cluster$Cluster
      
      #actually running k-nearest neighbors
      
      assigned_cluster <- class::knn(train,test,cl,k = k_for_knn)
      names(assigned_cluster) <- rownames(test)
      
      clust_zero_indices <- match(names(assigned_cluster), to_clust$id)
      to_clust$Cluster[clust_zero_indices] <- assigned_cluster
     
      
      #assigning the results of knn to the cluster in the joined_dr tibble.  could probabably have done this and above steps with tidyverse
      clustered <- to_clust
      
    } else {
      clustered <- to_clust
    }
  } else {
    clustered <- to_clust 
  }
  
  
  clusters <- clustered$Cluster
  clusters <- as.numeric(levels(clusters))[clusters]
  clustered <- purrr::set_names(clusters, clustered$id) 
  
  
  return(clustered)
}




cluster_map_spinglass <- function(
  graph, 
  spinglass_parameters = list(weights = NULL,
                              vertex = NULL,
                              spins = 25,
                              parupdate = FALSE,
                              start.temp = 1,
                              stop.temp = 0.01,
                              cool.fact = 0.99,
                              update.rule = c("config", "random", "simple"),
                              gamma = 1,
                              implementation = c("orig", "neg"),
                              gamma.minus = 1)  
){
   
  apply_spinglass <- function(graph_in, parameters_in = spinglass_parameters){
   
    full_sg_params <- c(list(graph_in),parameters_in)
     return(do.call(what = "cluster_spinglass", args = full_sg_params))
    
  }  
    
  
  components <- igraph::decompose(graph)

   clustered <- purrr::pmap(.l = list(graph_in = components), .f = apply_spinglass) %>%
                       map(.f = igraph::membership)
   
   
   
   for(d in seq_along(clustered)){
     
     if(d == 1 ){
       clust_final <- clustered[[d]]
       to_add  <- 0
       } else {
         
         to_add <- max(clust_final)
      
          clust_final <- c(clust_final , (clustered[[d]] + to_add))
     }
   }
return(clust_final)
                      
    
}












tidygraph_to_igraph <- function(tidygraph_in,id_string = "id"){
  #tidygraph in is the tidygraph we are converting to igraph
  #id string is the column in the nodes that is the ID
  
  is_directed <- with_graph(tidygraph_in, graph_is_directed())
  nodes <- tidygraph_in %>% activate(nodes) %>% as_tibble()  
  edges <- tidygraph_in %>% activate(edges) %>% as_tibble()
  
  
  if((id_string %in% names(nodes)) ){
    nodes <- nodes %>%dplyr::relocate(.data[[id_string]], .before = everything())
    edges$from <- nodes[[id_string]][edges$from]
    edges$to   <- nodes[[id_string]][edges$to]
  }
  
  
  
  if((dim(nodes)[1] != 0) & (dim(nodes)[2] != 0)){
    return(graph_from_data_frame(d = edges, directed = is_directed, vertices = nodes))
  } else{
    return(graph_from_data_frame(d = edges, directed = is_directed))
    
  }
}

igraph_to_tidygraph <- function(igraph_in, id_string =  "id")  {
  return(
    as_tbl_graph(igraph_in)  %>% activate(nodes) %>%  dplyr::rename({{id_string}} := name) 
  )
}


enforce_tidygraph <- function(network_in){
  classes <- class(network_in)
  if(classes[[1]] != "tbl_graph"){
    return(igraph_to_tidygraph(network_in))
  } else(return(network_in))
}

enforce_igraph <- function(network_in){
  classes <- class(network_in)
  if(classes[[1]] != "igraph"){
    return(tidygraph_to_igraph(network_in))
  } else(return(network_in))
}



collapse_igraph <- function(graph_in){  igraph::simplify(graph_in)}
undirect_igraph <- function(graph_in){igraph::as.undirected(graph = graph_in, mode = "each" )}
collapse_undirect_igraph <- function(graph_in){ igraph::as.undirected(graph = graph_in, mode ="collapse")}
throw_error_graph <- function(graph_in){stop("network layout failed")}

apply_clust_alg_igraph <- function(network,clust_alg,parameters){#applies an arbitrary function to an igraph
  clust_alg <- str_c("cluster_",clust_alg)
  
  
  
  
  
  alg_in_igraph <- exists(clust_alg, where = "package:igraph", mode = "function")
  if (alg_in_igraph){
    
 
      all_params <- c(list(graph = network), parameters)
      communities_out <- do.call(what = clust_alg, args = all_params)  
  
    
    
    cluster_out <- memberhsip(communities_out)
    network <- set_vertex_attr(graph = network, name = "Cluster", value = cluster_out)
  } else {
    
    
    all_params <- c(list(graph = network), parameters)
    cluster_out <-  do.call(what = clust_alg, args = all_params)
    network <- set_vertex_attr(graph = network, name = "Cluster", value = cluster_out)
  }
  return(network)
}
