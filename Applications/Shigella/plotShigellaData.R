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
require(ggtext)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

names = c("chromosome", "pINV", "spA", "spB", "spC")
order = c("NC_007384", "NC_007385", "NC_009345", "NC_009346", "NC_009347")


colors = c("spA"="#0072B2",
            "entire spA"="#0072B2",
           "strAB + sul\n+ flanking"="#0072B2",
           "AMR genes\nonly"="#0072B2",
           'spB'="#009E73",
           'spC'="#D55E00",
           'pINV'="#E69F00")

# # prevents files hanging
# system("rm combined/*.*")
# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("*rep0.*", sep=""), full.names = TRUE)
for (i in seq(1, length(log))){
  print(log[i])
  if (grepl("trees", log[[i]]) || grepl("log", log[[i]])){
    in_command <- " -b 10 -resample 500000 -log"
    for (j in seq(0,2)){
      in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), log[i]), sep="")
    }

    out_command = gsub("rep0_", "", log[i])
    out_command = gsub("out", "combined", out_command)
#
    combined_command = out_command
    combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
    # combine the trees
# 
#     system(intern=T,ignore.stdout = TRUE, paste("/Applications/BEAST\\ 2.7.6/bin/logcombiner", in_command, combined_command, sep=" "))
    if (grepl("trees", log[[i]])){
      inc = gsub("-o ","",combined_command)
      if (grepl("network", log[[i]])){
        if (grepl("SonFlex", log[[i]])){
#           system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher PlasmidTransferCount -burnin 0 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt", inc, gsub(".trees",".txt",inc), sep=" "))
#           system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher PlasmidTreeMapper -burnin 0 -followSegment 1 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt", inc, gsub(".trees",".mapped.trees",inc), sep=" "))
#           system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/treeannotator -height mean", gsub(".trees",".mapped.trees",inc), gsub(".trees",".mapped.tree",inc), sep=" "))
        }
#         system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentNetworkSummarizer -positions mcc -followSegment 0", inc, gsub(".trees",".tree",inc), sep=" "))
      # system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher PlasmidLossRate -maxTipDistance 10 -burnin 0", inc, gsub(".trees",".loss.txt",inc)))
# 
      }else{
#         system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/treeannotator  -height mean", inc, gsub(".trees",".tree",inc), sep=" "))
      }
    }
  }
}
# 
# 
# system(intern=F, "/Applications/BEAST\\ 2.7.6/bin/applauncher PlasmidTreeMapper -burnin 0 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt combined/SonFlex_equal.network.trees combined/SonFlex_equal.mapped.trees")
# system(intern=T,"/Applications/BEAST\\ 2.7.6/bin/treeannotator -burnin 0 -height mean combined/SonFlex_equal.mapped.trees combined/SonFlex_equal.mapped.tree")



# Plots the rate at which each plasmid in the Shigella sonnei dataset is 
# transferred between bacterial lineages. The rate is per lineage and year.
names_spA = c(labels(colors)[[2]], labels(colors)[[3]], labels(colors)[[4]])
spAOnly = data.frame()
combinedFrame = data.frame()
for (son in c(1,2,3)){
  t <- read.table(paste("combined/Sonnei",son,"_equal.log", sep=""), header=TRUE, sep="\t")
  spAOnly = rbind(spAOnly, data.frame(y=t[,"plasmidTransferRate.2"], name=names_spA[[son]]))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"plasmidTransferRate.1"], x=names[2], fill="pINV",title=gsub("\n"," ",paste(names_spA[[son]]))))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"plasmidTransferRate.2"], x=names[3], fill=names_spA[[son]],title=gsub("\n"," ",paste(names_spA[[son]]))))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"plasmidTransferRate.3"], x=names[4], fill="spB",title=gsub("\n"," ",paste(names_spA[[son]]))))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"plasmidTransferRate.4"], x=names[5], fill="spC",title=gsub("\n"," ",paste(names_spA[[son]]))))
}
# reorder the levels of the titles in combinedFrame
combinedFrame$title = factor(combinedFrame$title, levels=unique(combinedFrame$title))
# plot all plots next to each other

p = ggplot(combinedFrame)+
  geom_violin(aes(x=x,y=y, fill=fill))+
  theme_minimal()+
  coord_cartesian(ylim=c(0,0.05))+
  scale_x_discrete(limits=names[2:5], labels=c(names[2], "spA", names[4], names[5]))+
  theme(legend.position = "none")+
  scale_fill_manual(values=colors)+
  ylab("Plasmid transfer rate in events per year")+
  xlab("")+
  facet_wrap(~title, ncol=2)
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_rate_equal.pdf", sep=""),width=8, height=5)


# plot spA's only
p = ggplot(spAOnly)+
  geom_violin(aes(x=name,y=y), fill=colors[1])+
  theme_minimal()+
  coord_cartesian(ylim=c(0,0.05))+
  scale_x_discrete(limits=c(names_spA[[1]], names_spA[[2]], names_spA[[3]]))+
  theme(legend.position = "none")+
  # scale_fill_manual(values=colors)+
  ylab("Plasmid transfer rate per year")+
  xlab("")
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_rate_spA_comp.pdf", sep=""),width=5, height=3)

spAOnly = data.frame()
combinedFrame = data.frame()
for (son in c(1,2,3)){
  t <- read.table(paste("combined/Sonnei",son,"_equal.log", sep=""), header=TRUE, sep="\t")
  spAOnly = rbind(spAOnly, data.frame(y=t[,"obsPlasmidTransferEvents.2"], name=names_spA[[son]]))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"obsPlasmidTransferEvents.1"], x=names[2], fill="pINV",title=gsub("\n"," ",paste(names_spA[[son]]))))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"obsPlasmidTransferEvents.2"], x=names[3], fill=names_spA[[son]],title=gsub("\n"," ",paste(names_spA[[son]]))))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"obsPlasmidTransferEvents.3"], x=names[4], fill="spB",title=gsub("\n"," ",paste(names_spA[[son]]))))
  combinedFrame = rbind(combinedFrame, data.frame(y=t[,"obsPlasmidTransferEvents.4"], x=names[5], fill="spC",title=gsub("\n"," ",paste(names_spA[[son]]))))
}

combinedFrame$title = factor(combinedFrame$title, levels=unique(combinedFrame$title))

p = ggplot(combinedFrame)+
  geom_violin(aes(x=x,y=y, fill=fill))+
  theme_minimal()+
  coord_cartesian(ylim=c(0,40))+
  scale_x_discrete(limits=names[2:5], labels=c("pInv", "spA", names[4], names[5]))+
  theme(legend.position = "none")+
  scale_fill_manual(values=colors)+
  ylab("Inferred number of\nplasmid transfer events")+
  xlab("")+
  facet_wrap(~title, ncol=2)

plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_events_equal.pdf", sep=""),width=8, height=5)

p = ggplot(spAOnly)+
  geom_violin(aes(x=name,y=y), fill=colors[1])+
  theme_minimal()+
  coord_cartesian(ylim=c(0,40))+
  scale_x_discrete(limits=c(names_spA[[1]], names_spA[[2]], names_spA[[3]]))+
  theme(legend.position = "none")+
  # scale_fill_manual(values=colors)+
  ylab("Inferred number of plasmid transfer events")+
  xlab("")
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_events_spA_comp.pdf", sep=""),width=5, height=3)

all_plots <- list()
combinedFrame = data.frame()
for (son in c(1,2,3)){
  t.loss <- read.table(paste("combined/Sonnei", son,"_equal.network.loss.txt", sep=""), header=TRUE, sep="\t")
  combinedFrame = rbind(combinedFrame, data.frame(nrevents=t.loss[,"nrevents"], length=t.loss[,"length"],
                                                  segment=t.loss[,"segment"], title=gsub("\n"," ",names_spA[[son]])))
  plot(p)
}
combinedFrame$title = factor(combinedFrame$title, levels=unique(combinedFrame$title))
p = ggplot(combinedFrame, aes(x = segment, y=nrevents/length, group=segment)) + 
  geom_violin(aes(fill=names[segment+1]))+
  theme_minimal() + 
  scale_y_log10() +
  scale_x_continuous(breaks=seq(1,length(order)-1), labels=c(names[2], names_spA[[1]], names[4], names[5]))+
  theme(legend.position = "none")+
  scale_fill_manual(values=colors)+
  ylab("rate of plasmid getting lost per year")+
  xlab("") +
  facet_wrap(~title, ncol=2)
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_loss_equal.pdf", sep=""),width=4.5, height=4)

# Initialize the list for storing the plots
all_plots <- list()
legend_plot <- NULL

for (son in c(1,2,3)) {
  t <- read.table(paste("combined/Sonnei", son, "_equal.log", sep=""), header=TRUE, sep="\t")
  seq_info <- read.table("shigella_lengthrate.tsv", header=TRUE, sep="\t")
  dat <- data.frame()
  prior_vals <- rlnorm(1000000, -9.6009, 2)
  
  for (i in which(grepl("mutationRate", labels(t)[[2]]))) {
    accession <- strsplit(labels(t)[[2]][i], split="\\.")[[1]][[3]]
    name <- names[which(order == accession)]
    snps <- as.numeric(seq_info[seq_info$plasmid == accession, "snps"][1])
    totlength <- as.numeric(seq_info[seq_info$plasmid == accession, "lengthforrate"][1])
    individual.table <- read.table(paste("./combined/Sonnei", son, "_equal_", accession, ".log", sep=""), header=TRUE, sep="\t")
    
    if (startsWith(accession, "NC_007")) {
      rates1 <- t[,"clockRate.c"] * t[, i] * snps / totlength
      rates2 <- individual.table[, "clockRate"] * snps / totlength
    } else {
      rates1 <- t[,"clockRate.c"] * t[, i]
      rates2 <- individual.table[, "clockRate"]
    }
    
    hpd1 <- HPDinterval(as.mcmc(rates1))
    hpd2 <- HPDinterval(as.mcmc(rates2))
    
    dat <- rbind(dat, data.frame(x = which(order == accession), name = name, evol.rate = median(rates1), lower = hpd1[1,"lower"], upper = hpd1[1,"upper"], method = "network", offset = -1))
    dat <- rbind(dat, data.frame(x = which(order == accession), name = name, evol.rate = median(rates2), lower = hpd2[1,"lower"], upper = hpd2[1,"upper"], method = "individual trees", offset = 1))
  }
  
  names_adapted <- c("chromosome", "pINV", "spA", "spB", "spC")
  names_adapted[[3]] <- names_spA[[son]]
  
  dat$label <- paste(signif(dat$evol.rate, digits = 2))
  
  expSup <- function(w, digits = 3) {
    base <- signif(w / 10^floor(log10(abs(w))), digits = digits)
    exponent <- floor(log10(abs(w)))
    bquote(.(base) %*% 10^.(exponent))
  }
  
  dat$label <- sapply(dat$evol.rate, function(x) expSup(x))
  
  p <- ggplot(dat, aes(x = x - offset * 0.05, y = evol.rate, color = method)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    geom_text(aes(x = x - offset * 0.2, y = evol.rate, label = label), angle = 90, color = "black", size = 3, parse = TRUE) +
    theme_minimal() + 
    scale_y_log10() +
    scale_x_continuous(breaks = seq(1, length(order)), labels = names_adapted) +
    ylab("clock rate per site and year") + xlab("") +
    scale_color_manual(values = c("#448C82", "#C35F35")) +
    coord_cartesian(ylim = c(0.0000001, 0.0001)) +
    theme(legend.position = "none")
  
  all_plots[[son]] <- p
  
  if (is.null(legend_plot)) {
    legend_plot <- get_legend(p + theme(legend.position = "top"))
  }
  
  ggsave(plot = p, paste("../../../Plasmids-Text/Figures/Sonnei", son, "_rate.pdf", sep = ""), width = 6, height = 4)
}

# Combine the individual plots and the legend into a single plot
combined_plot <- ggarrange(plotlist = c(all_plots, list(legend_plot)), ncol = 2, nrow = 2, heights = c(1, 1))
plot(combined_plot)

# Save and display the combined plot
ggsave(plot = combined_plot, filename = "../../../Plasmids-Text/Figures/combined_rate_plot.pdf", width = 9, height = 6)


t <- read.table("combined/SonFlex_equal.network.txt", header=TRUE, sep="\t")
tred = t[t$toheight<50 & t$fromheight<50,]


dat.counts = data.frame()
for (i in unique(tred$number)){
  dat.counts = rbind(dat.counts, data.frame(from="S. flexneri",to="into S. sonnei",plasmid="MDR plasmid",numbers=sum(tred$number==i & tred$from==1 & tred$to==0 & tred$segment==1), run=i))
  dat.counts = rbind(dat.counts, data.frame(from="S. sonnei",to="into S. flexneri",plasmid="MDR plasmid",numbers=sum(tred$number==i & tred$from==0 & tred$to==1 & tred$segment==1), run=i))
}

tred2 = t[t$toheight<50 & t$fromheight>50,]

for (i in unique(tred2$number)){
  dat.counts = rbind(dat.counts, data.frame(from="unknown", to="into S. sonnei",plasmid="MDR plasmid",numbers=sum(tred2$number==i & tred2$to==0 & tred2$segment==1), run=i))
  dat.counts = rbind(dat.counts, data.frame(from="unknown", to="into S. flexneri",plasmid="MDR plasmid",numbers=sum(tred2$number==i & tred2$to==1 & tred2$segment==1), run=i))
}


p = ggplot(data=dat.counts) +
  geom_histogram(aes(x=numbers, fill=from),position = "dodge", binwidth=0.5, color="black", size=0.2 )+
  theme_minimal()+
  ggtitle("Plasmids introduced into S. sonnei or S. flexneri") +
  xlab("number of events") +
  ylab("probability density")+
  facet_grid(to~.) +
  scale_fill_manual(name="from", values=c("#999933", "#44AA99", "#DDDDDD")) +
  scale_x_continuous(breaks=seq(0,10,1), limit=c(0,8)) +
  theme(        axis.text.y=element_blank(),  #remove y axis labels
                axis.ticks.y=element_blank(),  #remove y axis ticks
                legend.position="top"
  )
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/SonFlex_Crossspecies.pdf", sep=""),width=4.2, height=3.5)


networkfiles <- list.files(path="./combined/", pattern=paste("*.*network.trees", sep=""), full.names = TRUE)

nr_events=1
for (n in seq(2, length(networkfiles))){
  # if the network files contains SonFlex, add the clade files as input
  # if (grepl("SonFlex", networkfiles[[n]])){
  #   system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher LineagesThroughTime -conditionOnChromosome true -timepoints 0.1,20,0.1 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt", networkfiles[[n]], gsub(".trees",".ltt.txt",networkfiles[[n]]), sep=" "))
  # }else{
  #   system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher LineagesThroughTime -conditionOnChromosome true  -timepoints 0.1,20,0.1", networkfiles[[n]], gsub(".trees",".ltt.txt",networkfiles[[n]]), sep=" "))
  # }
  
  ## Plot the lineages through time for with and without plasmids
  t <- read.table(gsub(".trees",".ltt.txt",networkfiles[[n]]), header=TRUE, sep="\t")
  
  if (length(labels(t)[[2]])==7){
    colors = c("spA"="#0072B2",
               'spB'="#009E73",
               'spC'="#D55E00")
    
    names = c("chromosome", "pINV", "spA", "spB", "spC")
    use_segments = c(3,4,5)
  }else{
    names = c("chromosome", "pKSR100 both species","cinSonnei", "cInFlex","pKSR100 in S. sonnei", "pKSR100 in S. flexneri")
    use_segments = c(2,5,6)
    colors = c("pKSR100 both species"="#E69F00",
               'pKSR100 in S. sonnei'="#44AA99",
               'pKSR100 in S. flexneri'="#999933")
  } 

  # compute the Hpds using the different iterations for all segments. The segments
  # are denoted as "segmentprop_0" "segmentprop_1" "segmentprop_2" "segmentprop_3" "segmentprop_4"
  ltt = data.frame()
  for (segments in use_segments){
    for (time in unique(t$time)){
      values = t[t$time==time, labels(t)[[2]][segments+2]]
      # check if there are any NaN values, if so, skip this step
      if (any(is.na(values))){
        next
      }
      hpd = HPDinterval(as.mcmc(values))
      hpd.5 = HPDinterval(as.mcmc(values), prob=0.5)
      ltt = rbind(ltt, data.frame(time=time, 
                                  segment=names[segments], lower.5=hpd.5[1,"lower"], upper.5=hpd.5[1,"upper"], 
                                  lower.95=hpd[1,"lower"], upper.95=hpd[1,"upper"]))
    }
  }

  # subset ltt to only use segments in colors
  p = ggplot(ltt[ltt$segment %in% names[use_segments],])+
    geom_ribbon(aes(x=as.Date("2020-12-01")-365*time,
                    ymin=lower.5, ymax=upper.5,
                    fill=segment), alpha=1)+
    geom_ribbon(aes(x=as.Date("2020-12-01")-365*time,
                    ymin=lower.95, ymax=upper.95,
                    fill=segment), alpha=0.5)+
    geom_line(aes(x=as.Date("2020-12-01")-365*time, y=lower.5, color=segment))+
    geom_line(aes(x=as.Date("2020-12-01")-365*time, y=upper.5, color=segment))+
    scale_fill_manual(name="", values=colors)+
    scale_color_manual(name="", values=colors, guide="none")+
    ylab("proportion of lineages\nwith plasmid") +
    xlab('')+
    # scale_x_date()+
    # coord_cartesian(xlim=c(as.Date("2020-12-01") - 20*365,as.Date("2020-12-01")), ylim=c(0,1)) +
    theme_minimal() +
    theme(legend.position="top")+
    guides(fill=guide_legend(nrow=1)) +
    scale_x_date(limits = c(as.Date("2000-01-01"), as.Date("2020-01-01")), date_labels="%Y")  
  plot(p)
  

  if (length(labels(t)[[2]])==7){
    # add plasmid prevalence over time
    t_seq = read.table("sonnei_plasmid_info.tsv", header=T, sep="\t")
    t_seq$date = as.Date(t_seq$times)
    
    unique_times = seq(min(t_seq$date), max(t_seq$date),"1 month")
    
    moving_dat = data.frame()
    for (i in seq(3, length(unique_times)-3)){
      dates_inavg = seq(unique_times[i-2], unique_times[i+2],"1 month")
      values = t_seq[is.element(t_seq$date, dates_inavg),]
      moving_dat = rbind(moving_dat, data.frame(date=unique_times[i], mean=mean(values$NC_009345),plasmid="spA"))
      moving_dat = rbind(moving_dat, data.frame(date=unique_times[i], mean=mean(values$NC_009346),plasmid="spB"))
      moving_dat = rbind(moving_dat, data.frame(date=unique_times[i], mean=mean(values$NC_009347),plasmid="spC"))
    }
    
    p=p+
      # geom_bar(data=t_seq, aes(x=date), fill="#DDDDDD")+
      # geom_ma(ma_fun = SMA, n = 30, aes(y=NC_009345*40, fill="spA"))+
      geom_line(data=moving_dat, aes(x=date,y=mean, color=plasmid), linetype="dashed")+
      theme_minimal() +xlab("")+
      # scale_color_manual(name="", values=colors)+
      theme(legend.position="top")
    plot(p)
  }
  
  ggsave(plot=p,paste("../../../Plasmids-Text/Figures/", gsub(".network.trees",".ltt.pdf",strsplit(networkfiles[[n]], split='\\/')[[1]][[4]]), sep=""),
                                                              width=4, height=3)
}
