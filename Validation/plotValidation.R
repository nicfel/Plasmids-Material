######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
# needed to calculate ESS values
library(coda)
library("methods")
library(colorblindr)
require(ggtree)
library(ggpubr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


t1 <- read.table("./out/simulate_serial5taxon5plasmids.log", header=TRUE, sep="\t")
t2 <- read.table("./out/simulate_serial5taxon10plasmids.log", header=TRUE, sep="\t")
t3 <- read.table("./out/testAll_serial5taxon5plasmids.log", header=TRUE, sep="\t")
t4 <- read.table("./out/testAll_serial5taxon10plasmids.log", header=TRUE, sep="\t")


dat = data.frame(method="simulated", noplasmids="5 plasmids", x=t1$network.height, value="height")
dat = rbind(dat, data.frame(method="simulated", noplasmids="10 plasmids", x=t2$network.height, value="height"))
dat = rbind(dat, data.frame(method="sampled", noplasmids="5 plasmids", x=t3$network.height, value="height"))
dat = rbind(dat, data.frame(method="sampled", noplasmids="10 plasmids", x=t4$network.height, value="height"))

dat = rbind(dat, data.frame(method="simulated", noplasmids="5 plasmids", x=t1$network.totalLength, value="length"))
dat = rbind(dat, data.frame(method="simulated", noplasmids="10 plasmids", x=t2$network.totalLength, value="length"))
dat = rbind(dat, data.frame(method="sampled", noplasmids="5 plasmids", x=t3$network.totalLength, value="length"))
dat = rbind(dat, data.frame(method="sampled", noplasmids="10 plasmids", x=t4$network.totalLength, value="length"))

dat = rbind(dat, data.frame(method="simulated", noplasmids="5 plasmids", x=t1$network.reassortmentNodeCount, value="no. plasmid transfer events"))
dat = rbind(dat, data.frame(method="simulated", noplasmids="10 plasmids", x=t2$network.reassortmentNodeCount, value="no. plasmid transfer events"))
dat = rbind(dat, data.frame(method="sampled", noplasmids="5 plasmids", x=t3$network.reassortmentNodeCount, value="no. plasmid transfer events"))
dat = rbind(dat, data.frame(method="sampled", noplasmids="10 plasmids", x=t4$network.reassortmentNodeCount, value="no. plasmid transfer events"))



p=ggplot(data=dat)+
  geom_density(aes(x=x,color=method, fill=method), adjust=1, alpha=0.2)+
  facet_grid(noplasmids~value, scales="free_x")+
  theme_minimal() +xlab("")
plot(p)

ggsave(plot=p,paste("../../Plasmids-Text/Figures/validation.pdf", sep=""),width=9, height=4)


