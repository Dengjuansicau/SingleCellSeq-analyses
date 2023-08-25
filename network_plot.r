###network plot###
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
mynet <- list()
phase <- c("em_E45_E60","em_E75_E105","em_E120_E135","non-ruminant","transition","ruminant")
s_short_levels <- c("PBC","BC1_FABP4+","BC2_IGFBP6+","BC3_CA1+","BC4_SELENBP1+","SC1_AGR2+","SC2_IGFBP6+","SC3_IGFBP2+","SC4_ALDH1A1+","SC5_KRT15+","SC6_FABP4+","SC7_KRT17+","GC1_KRTDAP+","GC2_S100A8+","MC","E","ECC","Stroma1_POSTN+","Stroma2_PTN+","M","Fib1_DLK+","Fib2_PLAC9+","Fib3_APOD+","Fib4_MFAP5+","NEC","Neu","ICC","SMC","B","T")
s_col <- c("skyblue2","deepskyblue2","dodgerblue2","mediumslateblue","mediumpurple2","turquoise3","turquoise","paleturquoise","seagreen4","seagreen2","palegreen4","palegreen","olivedrab2","yellowgreen","goldenrod1","goldenrod4","khaki2","darkgoldenrod2","darkorange","darkorange2","darkorange4","chocolate3","coral4","coral2","purple4","purple2","purple","orchid4","orchid2","plum4")
for(i in 1:6){
  setwd(paste(phase[[i]],"\\out",sep=""))
  mynet[[i]] <- read.table(paste(phase[[i]],"_s_count_network.txt",sep=""),header = T,sep="\t")
}

for (i in 1:6){
  mynet[[i]] <- subset(mynet[[i]], count>0)
  mynet[[i]]$SOURCE <- factor(mynet[[i]]$SOURCE, levels = s_short_levels)
  mynet[[i]] <- mynet[[i]][order(mynet[[i]]$SOURCE),]
  }
#########Set colors######
order <- unique(as.numeric(mynet[[1]]$SOURCE))
allcolour <- s_col[order]
index <- rank(order)
net<- graph_from_data_frame(mynet[[1]]) 
karate_groups <- cluster_optimal(net) 
coords <- layout_in_circle(net, order = index)
####Set the thickness of the edges###
E(net)$width  <- E(net)$count/10

edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

node <- unique(mynet[[1]]$SOURCE)
node <- as.vector(node) 
for (i in 1: length(node)){
  E(net)[map(node,function(x) {
    get.edge.ids(net,vp = c(node[i],x))
  }) %>% unlist()]$color <- allcolour[i]
}
E(net)$color=rgb(0.3,0.6,1,0.2)
#####Sets the size####
node_count <- read.table("interaction_count.txt",header = T,sep="\t")
node_count$X <- factor(node_count$X,levels=node)
size <- node_count[order(node_count$X),]
size$size <- size$all_sum/30

plot(net, edge.arrow.size=0, 
       edge.curved=0.2,
      # edge.label = vertex.color, 
     vertex.label.cex=.5,
    vertex.color=allcolour,
     vertex.size = size$size,
       vertex.frame.color=NA,
       vertex.label.color="black",
     vertex.label.dist=0,  
     vertex.label.degree=0, 
       layout = coords) 
