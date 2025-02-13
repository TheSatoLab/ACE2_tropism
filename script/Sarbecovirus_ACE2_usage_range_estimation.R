## set up
rm(list=ls())
setwd('~/Desktop/240503IA (Bat)/brown model/github')
library(phytools)

####READ IN DATA####

#phylogeny
tree<-read.nexus('imput/Sarbecovirus_RBD.fasta.tree')
tree$tip.label <- gsub("'", "", tree$tip.label)
tree$edge.length[which(tree$edge.length==0)] <- tree$edge.length[which(tree$edge.length==0)] + 0.000000001
tree_reroot <- reroot(tree, 79)
tree_rotate <- rotateNodes(tree_reroot, c(56:105))

#trait data
inf <- read.csv("imput/Sarbecovirus_infectivity_profile.csv", header = T, row.names = 1)

##log2 and scale
Inf.log2 <- log(inf+1,2)
Inf.log2.scale <- scale(Inf.log2)
Inf.median <- apply(Inf.log2.scale,1,median,row.names=1)

#estimate ancestral states
fit<-fastAnc(tree_reroot,Inf.median,vars=TRUE,CI=TRUE)

#create the figure
obj<-contMap(tree_rotate,Inf.median,plot=FALSE)
obj$cols[1:1001]<-colorRampPalette(c("blue","white","red"), space="Lab")(1001)
pdf("output/Sarbecovirus_ACE2_usage_range_estimation.pdf")
plot(obj,legend=0.7*max(nodeHeights(tree)),fsize=c(0.7,0.9))
dev.off()
