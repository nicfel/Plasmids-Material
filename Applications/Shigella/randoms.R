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


t = read.table("./out/Flexneri_rep0.log", header=T, sep="\t")
# t2 = read.table("./out/Flexneri_rep1.log", header=T, sep="\t")
# t3 = read.table("./out/Flexneri_rep2.log", header=T, sep="\t")

# t = rbind(t1,t2,t3)



ratio = t$NC_004337.tree.height/t$mutationRate.s.NC_004337
ratio = t$NC_004337.tree.treeLength*t$mutationRate.s.NC_004337


plot(ratio)

require("lubridate")


t = read.table("/Users/nmueller/Downloads/h3n2_ha_metadata.txt", sep="\t", header=T)
t$d = floor_date(as.Date(t$date), "year")
t = t[!is.na(t$d),]
p = ggplot(t, aes(x=d)) + geom_histogram(binwidth=365, color="white", fill="black") + theme_minimal() + ylab("number of influenza H3N2\nsequences on gisaid") + xlab("")+
   scale_x_date(limits=c(as.Date("1990-01-01"), as.Date("2022-12-31"))) 
plot(p)
ggsave(plot=p, "/Users/nmueller/Downloads/sequence.pdf", width=6, height=3)


