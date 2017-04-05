# This example shows how to plot subset of graph.
# and change the node's color based on the node's lable.
# Input data: yeast protein-protein interactions

library(dplyr)
library(igraph)


# read ppi, top 20000
ppi_yeast<-read.table("data_yeast_top_20000.txt")
dh_input<-read.table("yeast-all-datehubs.txt")



#convert ppi to edges for igraph input
ppi_edges_vector<-c()
for(i in 1:nrow(ppi_yeast)){
  ppi_edges_vector<-c(ppi_edges_vector, c(as.character(ppi_yeast[i,1]), as.character(ppi_yeast[i,2])))  
}

# create graph
ppi_graph <- graph(ppi_edges_vector,directed = FALSE)

# plot protein complex
complex1 <-as.numeric(V(ppi_graph)[c("YIL033C","YJL164C","YKL166C","YPL203W")])
complex1_graph <- induced.subgraph(graph=ppi_graph,vids = complex1)
plot(complex1_graph)

# Set different color for nodes based on node's lable
V(complex1_graph)$color <- ifelse((V(g2)$name %in% c("YIL033C","YJL164C")), "orange","yellow")
plot(complex1_graph)

#calculate density
complex1_dens<-edge_density(complex1_graph, loops = FALSE)