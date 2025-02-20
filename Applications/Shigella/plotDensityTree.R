library(ggplot2)
library(coda)
library(phangorn)
library(ggtree)
library(phytools)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


plasmids = c("NC_009345","NC_009346","NC_009347")
for (p in plasmids){
  print(p)
  system(paste("/Applications/BEAST\\ 2.6.7/bin/logcombiner -b 0 -log ./combined/Sonnei_equal_", p, ".trees -log ./combined/Sonnei_equal.", p, ".trees -o ./combined/SonneiComp.", p,".trees", sep=""))
  
  
  
  Nnet <- read.nexus(paste("./combined/SonneiComp.",p,".trees", sep=""))
  Nnet2 <- read.nexus(paste("./combined/Sonnei_equal_",p,".trees", sep=""))
  
  
  nrtreeseach = 25
  use_trees = sort(c(sample(length(Nnet2), nrtreeseach, replace = FALSE), sample(length(Nnet)-length(Nnet2), nrtreeseach, replace = FALSE)+length(Nnet2)))
  
  
  tree.plot <- Nnet
  color=list()
  cols = c("#023FA5", "#8E063B")
  
  for (i in seq(1,length(tree.plot))){
    if (i <=length(Nnet2)){
      color[i] = cols[1]
    }else{
      color[i] = cols[2]
    }
  }
  
  
  
  
  par(mar=c(1,1,1,1))
  
  pdf(paste("../../../Plasmids-Text/Figures/",  p, "_densitree.pdf", sep=""), width=5, height=5, onefile=FALSE)
  
  
  source("densitreecopy.R")
  densiTree(tree.plot, col=color, direction="leftwards", alpha=0.1, width=1, scale.bar=F, underscore=T,flipval=length(Nnet2), usetrees=use_trees)
  dev.off()
}


#system("/Applications/BEAST\\ 2.6.7/bin/treeannotator -b 0 -heights keep -target ./combined/SonFlex_weighted.LN624486.tree ./combined/SonFlex_LN624486.trees ./combined/SonFlex_LN624486.target.tree")


