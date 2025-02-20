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

names = c("chromosome", "pINV", "spA", "spB", "spC")
order = c("NC_007384", "NC_007385", "NC_009345", "NC_009346", "NC_009347")


colors = c("spA"="#0072B2",
           'spB'="#009E73",
           'spC'="#D55E00",
           'pINV'="#E69F00")

# prevents files hanging
system("rm combined/*.*")
# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("*rep0.*", sep=""), full.names = TRUE)
for (i in seq(1, length(log))){
# for (i in seq(1, 9)){

  print(log[i])
  if (grepl("trees", log[[i]]) || grepl("log", log[[i]])){
    in_command <- " -b 10 -resample 100000 -log"
    for (j in seq(0,2)){
      in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), log[i]), sep="")
    }

    out_command = gsub("rep0_", "", log[i])
    out_command = gsub("out", "combined", out_command)

    combined_command = out_command
    combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
    # combine the trees

    system(intern=T,ignore.stdout = TRUE, ignore.stderr = TRUE, paste("/Applications/BEAST\\ 2.6.7/bin/logcombiner", in_command, combined_command, sep=" "))

    if (grepl("trees", log[[i]])){
      inc = gsub("-o ","",combined_command)
      if (grepl("network", log[[i]])){
        if (grepl("SonFlex", log[[i]])){
          system(intern=T, paste("java -jar ../../Software/TransferCounter.jar -burnin 0 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt", inc, gsub(".trees",".txt",inc), sep=" "))
          system(intern=T, paste("java -jar ../../Software/PlasmidTreeMapper.jar -burnin 0 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt", inc, gsub(".trees",".mapped.trees",inc), sep=" "))

          system(intern=T, paste("/Applications/BEAST\\ 2.6.7/bin/treeannotator -heights mean", gsub(".trees",".mapped.trees",inc), gsub(".trees",".mapped.tree",inc), sep=" "))
        }
        system(intern=T, paste("java -jar ../../Software/Summarizer.jar -positions mcc -followSegment 0", inc, gsub(".trees",".tree",inc), sep=" "))
        system(intern=T, paste("java -jar ../../Software/LossCounter.jar -burnin 0", inc, gsub(".trees",".loss.txt",inc), sep=" "))
      }else{
        system(intern=T, paste("/Applications/BEAST\\ 2.6.7/bin/treeannotator  -heights mean", inc, gsub(".trees",".tree",inc), sep=" "))
      }
    }else{
    }
  }
}






t <- read.table("combined/Sonnei_equal.log", header=TRUE, sep="\t")

p = ggplot(t)+
  geom_violin(aes(x=names[2],y=plasmidTransferRate.1, fill="pINV"))+
  geom_violin(aes(x=names[3],y=plasmidTransferRate.2, fill="spA"))+
  geom_violin(aes(x=names[4],y=plasmidTransferRate.3, fill="spB"))+
  geom_violin(aes(x=names[5],y=plasmidTransferRate.4, fill="spC"))+
  theme_minimal()+
  coord_cartesian(ylim=c(0,0.05))+
  scale_x_discrete(limits=names[2:5])+
  theme(legend.position = "none")+
  scale_fill_manual(values=colors)+
  ylab("Plasmid transfer rate per year")+
  xlab("")
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_rate_equal.pdf", sep=""),width=4.5, height=4)


p = ggplot(t)+
  geom_violin(aes(x=names[2],y=obsPlasmidTransferEvents.1, fill="pINV"))+
  geom_violin(aes(x=names[3],y=obsPlasmidTransferEvents.2, fill="spA"))+
  geom_violin(aes(x=names[4],y=obsPlasmidTransferEvents.3, fill="spB"))+
  geom_violin(aes(x=names[5],y=obsPlasmidTransferEvents.4, fill="spC"))+
  theme_minimal()+
  scale_x_discrete(limits=names[2:5])+
  coord_cartesian(ylim=c(0,60))+
  theme(legend.position = "none")+
  scale_fill_manual(values=colors)+
  ylab("Inferred number of plasmid transfer events")+
  xlab("")
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_events_equal.pdf", sep=""),width=4.5, height=4)








t.loss <- read.table("combined/Sonnei_equal.network.loss.txt", header=TRUE, sep="\t")

p = ggplot(data=t.loss, aes(x = segment, y=nrevents/length, group=segment)) + 
  geom_violin(aes(fill=names[segment+1]))+
  theme_minimal() + 
  # scale_y_log10() +
  scale_x_continuous(breaks=seq(1,length(order)-1), labels=names[2:length(names)])+
  theme(legend.position = "none")+
  scale_fill_manual(values=colors)+
  ylab("rate of plasmid getting lost per year")+
  xlab("") 
# scale_color_discrete_diverging()+coord_cartesian(ylim=c(0.0000001,0.0001))

plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Sonnei_plasmid_loss_equal.pdf", sep=""),width=4.5, height=4)




seq_info =  read.table("shigella_lengthrate.tsv", header=TRUE, sep="\t")
dat = data.frame()
prior_vals=rlnorm(1000000,-9.6009,2)
for (i in which(grepl("mutationRate", labels(t)[[2]]))){
  print(i)
  accession = strsplit(labels(t)[[2]][i], split="\\.")[[1]][[3]]
  
  snps = as.numeric(seq_info[seq_info$plasmid==accession, "snps"][1])
  totlength = as.numeric(seq_info[seq_info$plasmid==accession, "lengthforrate"][1])
  individual.table = read.table(paste("./combined/Sonnei_equal_",accession,".log",sep=""), header=TRUE, sep="\t")

  
  rates = t[, i]*snps/totlength
  hpd = HPDinterval(as.mcmc(rates))
  
  dat = rbind(dat, data.frame(x = which(order==accession),accession=accession, evol.rate = median(rates), lower=hpd[1,"lower"], upper=hpd[1,"upper"], method="network", offset=-1))
  
  
  rates = individual.table[, "clockRate"]*snps/totlength
  hpd = HPDinterval(as.mcmc(rates))
  dat = rbind(dat, data.frame(x = which(order==accession),accession=accession, evol.rate = median(rates), lower=hpd[1,"lower"], upper=hpd[1,"upper"], method="individual trees", offset=1))
  
  rates = prior_vals*snps/totlength
  hpd = HPDinterval(as.mcmc(rates))
  
  
  # dat = rbind(dat, data.frame(x = which(order==accession),accession=accession, evol.rate = mean(rates), lower=hpd[1,"lower"], upper=hpd[1,"upper"], method="prior", offset=1))
  
}



dat$label = gsub("e-","*10^-",paste(signif(dat$evol.rate, digits = 2), " (", signif(dat$lower, digits = 2),", ",signif(dat$upper, digits = 2), ")", sep=""))

# dat$name = factor(dat$name, levels=order)
require(ggtext)
p = ggplot(dat, aes(x = x-offset*0.05, y=evol.rate,color=method)) + 
  geom_pointrange(aes(ymin=lower,ymax=upper))+
  geom_richtext(aes(x=x-offset*0.2,y=evol.rate, label=label), angle = 90, color="black",fill = NA, label.color = NA, size=2)+
  theme_minimal() + 
  scale_y_log10() +
  scale_x_continuous(breaks=seq(1,length(order)), labels=names)+
  # facet_wrap(.~accession, scales="free")+
  ylab("rate of evolution pre site and year")+xlab("") +
  scale_color_discrete_diverging()+coord_cartesian(ylim=c(0.0000001,0.001))

plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Shigella_rate.pdf", sep=""),width=9, height=4)







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
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/Crossspecies.pdf", sep=""),width=4.2, height=3.5)





# add plasmid prevalence over time
t = read.table("sonnei_plasmid_info.tsv", header=T, sep="\t")
t$date = as.Date(t$times)

unique_times = seq(min(t$date), max(t$date),"1 month")

moving_dat = data.frame()
for (i in seq(2, length(unique_times)-1)){
  dates_inavg = seq(unique_times[i-1], unique_times[i+1],"1 month")
  values = t[is.element(t$date, dates_inavg),]
  moving_dat = rbind(moving_dat, data.frame(date=unique_times[i], mean=mean(values$NC_009345),plasmid="spA"))
  moving_dat = rbind(moving_dat, data.frame(date=unique_times[i], mean=mean(values$NC_009346),plasmid="spB"))
  moving_dat = rbind(moving_dat, data.frame(date=unique_times[i], mean=mean(values$NC_009347),plasmid="spC"))
}



p=ggplot(data=t, aes(x=date))+
  geom_bar(aes(x=date), fill="#DDDDDD")+
  # geom_ma(ma_fun = SMA, n = 30, aes(y=NC_009345*40, fill="spA"))+
  geom_line(data=moving_dat, aes(x=date,y=mean*40, color=plasmid))+
  scale_y_continuous(
    breaks=seq(0,40,10),
    limits=c(0,50),
    # Features of the first axis
    name = "sequenced cases per month",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~./40, name="3 month moving average\nof plasmid prevalence",breaks=seq(0,1,0.25))
  ) +
  scale_color_manual(values=colors)+
  theme_minimal() +xlab("")+
  theme(legend.position="none")
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/prevalence.pdf", sep=""),width=9, height=3)






networkfiles <- list.files(path="./combined/", pattern=paste("*.*network.trees", sep=""), full.names = TRUE)

nr_events=1
time_eval = seq(0.1,40,0.1)
for (n in seq(2, length(networkfiles))){
# for (n in c(2, 4,5)){
    system(intern=T, paste("java -jar ../../Software/LTT.jar -burnin 0", networkfiles[[n]], gsub(".trees",".ltt.txt",networkfiles[[n]]), sep=" "))
  
  ## Plot the lineages through time for with and without plasmids
  t <- read.table(gsub(".trees",".ltt.txt",networkfiles[[n]]), header=TRUE, sep="\t")
  
  times = data.frame()
  for (i in seq(1,max(t$segment))){
    dat = data.frame();
    for (j in seq(1,max(t$iteration),10)){
      print(j)
      vals_with = t[which(t$iteration==j & t$segments==i & t$withsegment=="with"),"ltt"]
      vals_without = t[which(t$iteration==j & t$segments==i & t$withsegment=="total"),"ltt"]
      vals = strsplit(vals_with, split="\\,")[[1]]
      vals_wo = strsplit(vals_without, split="\\,")[[1]]
      
      
      
      time_old = 0
      last_time = 0
      
      nr_lineages_old=0
      accumulator=0
      nCoal = 0
      
      nr_lineages_old_wo=0
      accumulator_wo=0
      nCoal_wo = 0
      
      
      
      for (k in seq(1,length(vals))){
        vals2 = strsplit(vals[k], split=":")[[1]]
        vals2_wo = strsplit(vals_wo[k], split=":")[[1]]
        if (as.numeric(vals2[2])<100){
          nr_lineages = as.numeric(vals2[1])
          nr_lineages_wo = as.numeric(vals2_wo[1])
          
          time = as.numeric(vals2[2])
          accumulator = accumulator + (time-time_old)*nr_lineages*(nr_lineages-1)/2
          accumulator_wo = accumulator_wo + (time-time_old)*nr_lineages_wo*(nr_lineages_wo-1)/2
          
          
          if (nr_lineages<nr_lineages_old){
            nCoal=nCoal+1
          }
          
          if (nr_lineages_wo<nr_lineages_old_wo){
            nCoal_wo=nCoal_wo+1
          }
          
          if (nCoal_wo==nr_events){
            
            
            if (length(which(time_eval>last_time & time_eval<time))>0){
              dat  = rbind(dat, data.frame(last_time=last_time, time=time, meant=mean(time, last_time),
                                           lineages = nr_lineages, lineages_tot = nr_lineages_wo,
                                           Ne=accumulator/nCoal, Ne_wo = accumulator_wo/nCoal_wo,
                                           iteration=j, segment=i))
            }
            nCoal=0
            accumulator=0
            nCoal_wo=0
            accumulator_wo=0
            last_time=time
          }
          
          
          
          nr_lineages_old=nr_lineages
          nr_lineages_old_wo=nr_lineages_wo
          time_old = time
        }
      }
    }
    for (tt in time_eval){
      vals = dat[which(dat$last_time<tt & dat$time>tt),]
      interval.95 = HPDinterval(as.mcmc(vals$lineages/vals$lineages_tot))
      interval.5 = HPDinterval(as.mcmc(vals$lineages/vals$lineages_tot), prob=0.5)
      times = rbind(times, data.frame(x=as.Date("2020-12-01")-tt*365, 
                                      lower.5=interval.5[1,"lower"],upper.5=interval.5[1,"upper"],
                                      lower.95=interval.95[1,"lower"],upper.95=interval.95[1,"upper"],
                                      segment=as.character(i)))
    }
  }
  
  if (max(t$segment)==1){
    curr_colors=c('1'='#117733')
  }else{
    curr_colors = c("2"="#0072B2",
               '3'="#009E73",
               '4'="#D55E00",
               '1'="#E69F00")
  }
  
  
  

  
  p = ggplot(times)+
    geom_ribbon(aes(x=x,ymin=lower.5, ymax=upper.5, fill=segment), alpha=1)+
    geom_ribbon(aes(x=x,ymin=lower.95, ymax=upper.95, fill=segment), alpha=0.5)+
    # geom_segment(aes(x=last_time_date, xend=time_date,
    #                  y=lineages/lineages_tot, 
    #                  yend=lineages/lineages_tot, group=iteration, color=segment), alpha=0.25) +
    scale_fill_manual(values=curr_colors)+
    ylab("proportion of ancestral lineages carrying plasmid") +
    xlab('')+
    scale_x_date()+
    coord_cartesian(xlim=c(as.Date("2020-12-01") - 20*365,as.Date("2020-12-01")), ylim=c(0,1)) +
    theme_minimal() +
    theme(legend.position = 'none')
  plot(p)
  
  ggsave(plot=p,paste("../../../Plasmids-Text/Figures/", gsub(".network.trees",".ltt.pdf",strsplit(networkfiles[[n]], split='\\/')[[1]][[4]]), sep=""),
                                                              width=9, height=4)
  
}
  
