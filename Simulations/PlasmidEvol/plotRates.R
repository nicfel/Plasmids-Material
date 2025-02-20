######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("inf", ".*rep1.log", sep=""), full.names = TRUE)

dat = data.frame()
Ne = data.frame()

true.rates = read.table("rates.csv", header=T, sep=",")
snps = read.table("snpcounts.txt", header=T, sep="\t")


for (i in seq(1,length(log))){
  print(i)
  fname1 = log[[i]]
  # read in the log file
  t.1 <- read.table(fname1, header=TRUE, sep="\t")
  t.2 <- read.table(gsub("rep1","rep2",fname1), header=TRUE, sep="\t")
  t.3 <- read.table(gsub("rep1","rep3",fname1), header=TRUE, sep="\t")
  
  # take a 10% burnin
  t.1 <- t.1[-seq(1,ceiling(length(t.1$posterior)/5)), ]
  t.2 <- t.3[-seq(1,ceiling(length(t.2$posterior)/5)), ]
  t.3 <- t.3[-seq(1,ceiling(length(t.3$posterior)/5)), ]
  
  if (length(t.1$Sample)<10){
    print(fname1)
    das
  }else if (length(t.2$Sample)<10){
    print(gsub("rep1","rep2",fname1))
    dsa
    
  }else if (length(t.3$Sample)<10){
    print(gsub("rep1","rep3",fname1))
    das
  }
  
  
  
  # combine the replictes
  t = rbind(t.1,t.2,t.3)
  
  ess <- effectiveSize(t)
  if (min(ess[grepl("mutationRate",labels(ess))])>50){
    tmp = strsplit(log[i], split="_")[[1]]
    
    if (length(tmp)==3){
      hpd = HPDinterval(as.mcmc(t$mutationRate1))
      
      
      snpsno1 = snps[snps$run==tmp[[2]] & snps$plasmid=="core1", "snpcount" ]
      
      dat = rbind(dat, data.frame(run=tmp[[2]], name="core1", method="network",
                                  rate = median(t$mutationRate1), lower=hpd[1,"lower"], upper=hpd[1,"upper"], offset=1*runif(1), 
                                  snp=snps[snps$run==tmp[[2]] & snps$plasmid=="core1", "snpcount" ]))
      
      hpd = HPDinterval(as.mcmc(t$mutationRate2))
      
      dat = rbind(dat, data.frame(run=tmp[[2]], name="plasmid1", method="network",
                                  rate = median(t$mutationRate2), lower=hpd[1,"lower"], upper=hpd[1,"upper"], offset=1*runif(1),
                                  snp=snps[snps$run==tmp[[2]] & snps$plasmid=="plasmid1", "snpcount" ]))
      hpd = HPDinterval(as.mcmc(t$mutationRate3))
      
      dat = rbind(dat, data.frame(run=tmp[[2]], name="plasmid2", method="network",
                                  rate = median(t$mutationRate3), lower=hpd[1,"lower"], upper=hpd[1,"upper"], offset=1*runif(1),
                                  snp=snps[snps$run==tmp[[2]] & snps$plasmid=="plasmid2", "snpcount" ]))
      hpd = HPDinterval(as.mcmc(t$mutationRate4))
      
      dat = rbind(dat, data.frame(run=tmp[[2]], name="plasmid3", method="network",
                                  rate = median(t$mutationRate4), lower=hpd[1,"lower"], upper=hpd[1,"upper"], offset=1*runif(1),
                                  snp=snps[snps$run==tmp[[2]] & snps$plasmid=="plasmid3", "snpcount" ]))
      
      
      
      # t.r = true.rates[true.rates$run==tmp[[2]],]
      # for (j in seq(1,13)){
      #   hpd = HPDinterval(as.mcmc(t[,sprintf("popSize.%d",j)]))
      #   
      #   Ne = rbind(Ne, data.frame(run=tmp[[2]], x=j, true = t.r[1,sprintf("popSize.%d",j)],
      #                             median=median(t[,sprintf("popSize.%d",j)]), lower=hpd[1,"lower"], upper=hpd[1,"upper"], method="network"))
      # }
    }else{
      hpd = HPDinterval(as.mcmc(t$clockRate))
      
      dat = rbind(dat, data.frame(run=tmp[[2]], name=tmp[[3]], method="individual trees", 
                                  rate = median(t$clockRate), lower=hpd[1,"lower"], upper=hpd[1,"upper"], offset=-1*runif(1),
                                  snp=snps[snps$run==tmp[[2]] & snps$plasmid==tmp[[3]], "snpcount" ]))
      
      # if (grepl("core", tmp[[3]])){
      #   t.r = true.rates[true.rates$run==tmp[[2]],]
      #   
      #   for (j in seq(1,13)){
      #     hpd = HPDinterval(as.mcmc(t[,sprintf("popSize.%d",j)]))
      #     
      #     Ne = rbind(Ne, data.frame(run=tmp[[2]], x=j, true = t.r[1,sprintf("popSize.%d",j)],
      #                               median=median(t[,sprintf("popSize.%d",j)]), lower=hpd[1,"lower"], upper=hpd[1,"upper"], method="core only"))
      #   }
      #   
      # }
    }
  }
      
}

p = ggplot(dat[dat$name!="core1" & dat$snp>0,], aes(x=snp+offset*0.1, y=rate, color=method))+
  geom_point(aes(y=rate))+
  geom_pointrange(aes(ymin=lower,ymax=upper), alpha=0.5)+
  facet_wrap(.~method)+
  theme_minimal()+
  scale_y_log10()+
  scale_x_log10()+
  geom_hline(yintercept = 0.0005)+
  coord_cartesian(ylim=c(0.000001,0.01))+
  # facet_wrap(method~.)+
  # theme(axis.text.x=element_blank())+
  xlab("number of variable sites")+ylab("evolutionary rate") +
  scale_color_manual(values=c("#448C82", "#C35F35"))+
  theme(legend.position = "none")
plot(p)

  
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/clock_rate_plasmids.pdf", sep=""),width=9, height=4)

rate_shifts = c(seq(0,20,2),100,200)

p = ggplot(Ne, aes(x=rate_shifts[x], y=true))+
  geom_line(aes(y=median, color="inferred"))+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill="inferred"),alpha=0.5)+
  facet_wrap(.~run)+
  geom_line(aes(color="true"), linetype="dashed")+
  theme_minimal()+
  coord_cartesian(xlim=c(0,20))
  # scale_x_log10()
# geom_hline(yintercept = 0.0005)
plot(p)

ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Ne.pdf", sep=""),width=12, height=8)





