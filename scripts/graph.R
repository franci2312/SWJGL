
# circ.graph --------------------------------------------------------------

# @param mod fitted model where mod$theta contains the results
# @param idx index of the cluster 
# @param vertex_names labels
# @param weighted.edges bool indicating whether to plot a weighted graph or not

circ.graph <- function(mod, idx, vertex_names = gene_names, weighted.edges = T){
  
  require(igraph); require(ggokabeito)
  col_gg <- c('#F94144','#F9C74F','#43AA8B','#CC79A7','#90BE6D','#F8961E','#577590','#000') 
  
  adj <- (mod$theta[[idx]] != 0) + 0
  if(weighted.edges){
    adj <- mod$theta[[idx]]
  }
  graph.hat <- igraph::graph_from_adjacency_matrix(adjmatrix = adj, 
                         mode = 'undirected', diag = F, weighted = weighted.edges)  
  vertex_degrees <- degree(graph.hat)
  vertex_n <- c()
  for(i in 1:length(vertex_degrees)){
    vertex_n[i] <- ifelse(vertex_degrees[i] == 0, "", vertex_names[i])
  }
  
  min_size <- 3 
  max_size <- 15 
  vertex_sizes <- (vertex_degrees - min(vertex_degrees)) / (max(vertex_degrees) 
                        - min(vertex_degrees)) * (max_size - min_size) + min_size
  layout_circle <- layout_in_circle(graph.hat)
  
  
  # commands for label positions
  offset <- 1.25
  label_positions <- layout_circle * offset
  
  
  angles <- atan(layout_circle[, 2] / layout_circle[, 1])
  angles_deg <- angles * 180 / pi
  angles_deg <- ifelse(angles_deg < 0, angles_deg + 360, angles_deg)
  par(mar = c(4.1, 4.1, 4.1, 4.1)) 
  if(weighted.edges){
    edge_weights <- E(graph.hat)$weight
    edge_weights <- ifelse(is.na(edge_weights), 0, edge_weights) 
    min_weight <- min(edge_weights)
    max_weight <- max(edge_weights)
    edge_widths <- (edge_weights - min_weight) / (max_weight - min_weight) * (5 - 1) + 1  
    plt <- plot.igraph(x = graph.hat, layout = layout_in_circle(graph = graph.hat),
                       vertex.size = vertex_sizes, vertex.label = NA, 
                       vertex.frame.color = NA, edge.curved = -0.4, 
                       vertex.color = col_gg[idx], ylim = c(-1.5,1.5),
                       edge.width = edge_widths)
    
  } else {
    plt <- plot.igraph(x = graph.hat, layout = layout_in_circle(graph = graph.hat),
                       vertex.size = vertex_sizes, vertex.label = NA, 
                       vertex.frame.color = NA, edge.curved = 0.2, 
                       vertex.color = col_gg[idx], ylim = c(-1.5,1.5))
  }
  for (i in seq_along(vertex_names)) {
    text(x = label_positions[i, 1], y = label_positions[i, 2], 
         labels = vertex_n[i], cex = 0.7, col = "#000000", srt = angles_deg[i], 
         adj = c(0.5, 0.5))
  }
  
  return(plt)
  
}


