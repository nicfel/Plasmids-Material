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
library(dplyr)

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

# function to compute the moving average for if a clade had a segment
movingavg_data_frame <- function(t, time_points, time_window, moving, name, name2){
  # Set up the moving average matrix for each time point and iteration
  moving_avg <- matrix(NA, nrow = length(time_points), ncol = max(t$iteration))
  for (i in 0:max(t$iteration)) {
    t_subset = t[t$iteration == i,]
    # Get all time points
    for (j in 1:length(time_points)) {
      # Get all values for the current time point
      moving_avg[j, i] <- mean(t_subset[t_subset$Time > time_points[j] & t_subset$Time < time_points[j] + time_window, "condition"], na.rm = TRUE)
    }
  }
  
  for (j in 1:length(time_points)) {
    # Get median and quantiles for the current time point
    median_val <- median(moving_avg[j, ], na.rm = TRUE)
    lower_quantile.1 <- quantile(moving_avg[j, ], 0.45, na.rm = TRUE)
    upper_quantile.1 <- quantile(moving_avg[j, ], 0.55, na.rm = TRUE)
    
    lower_quantile.5 <- quantile(moving_avg[j, ], 0.25, na.rm = TRUE)
    upper_quantile.5 <- quantile(moving_avg[j, ], 0.75, na.rm = TRUE)
    
    lower_quantile <- quantile(moving_avg[j, ], 0.025, na.rm = TRUE)
    upper_quantile <- quantile(moving_avg[j, ], 0.975, na.rm = TRUE)
    
    # Append results to the data frame
    moving <- rbind(moving, data.frame(time = time_points[j] + time_window/2, 
                                       median = median_val, 
                                       lower = lower_quantile, upper = upper_quantile,
                                       lower.5 = lower_quantile.5, upper.5 = upper_quantile.5,
                                       lower.1 = lower_quantile.1, upper.1 = upper_quantile.1, name=name,name2=name2))
  }
  return(moving)
}

# cpompute the average over the data
avg_data_frame = function(t, average, name, name2){
  vals=c()
  for (j in 0:max(t$iteration)){
    vals = c(vals, mean(t[t$iteration==j, "condition"]))
  }
  average = rbind(average, data.frame(mean=mean(vals),
                                        lower = quantile(vals, 0.025), upper = quantile(vals, 0.975),
                                        name=name,name2=name2))
}

sonnei_names = c("Sonnei1", "Sonnei2", "Sonnei3")

tau_vals =c(1)
# # compute the
# for (tau in tau_vals){
#   system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher SegmentLBI -burnin 0 -calculateAverageDifference true -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt -outTable combined/SonFlex_", tau, "_table.tsv -tau ",
#                          tau, " -windowSize 2.5 combined/SonFlex_equal.network.trees combined/SonFlex_" , tau, "_mapped.trees", sep=""))
# 
#   system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/treeannotator -height mean -burnin 0 combined/SonFlex_", tau, "_mapped.trees combined/SonFlex_", tau, "_mapped.tree", sep=""))
# 
#   for (son in sonnei_names){
#     system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher SegmentLBI -burnin 0 -calculateAverageDifference true -outTable combined/", son, "_", tau, "_table.tsv -tau ",
#                            tau, " -windowSize 2.5 combined/", son, "_equal.network.trees combined/", son, "_", tau, "_mapped.trees", sep=""))
#     system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/treeannotator -height mean -burnin 0 combined/", son, "_", tau, "_mapped.trees combined/", son, "_", tau, "_mapped.tree", sep=""))
#   }
# }
# 
# system(intern=F, "/Applications/BEAST\\ 2.7.6/bin/applauncher SegmentLBI -burnin 0 -calculateAverageDifference true -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt -outTable combined/SonFlex_equal.table.tsv combined/SonFlex_equal.network.trees combined/SonFlex_equal.mapped.trees")
# system(intern=T, "/Applications/BEAST\\ 2.7.6/bin/treeannotator -height mean -burnin 0 combined/SonFlex_equal.mapped.trees combined/SonFlex_equal.mapped.tree")
# system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ClusterSizeComparison -burnin 0 -cladeFileInput data/Sonnei_MonthYear_07082022.txt,data/ShigFlex_01092022_YearMonth.txt combined/SonFlex_equal.network.trees combined/SonFlex_comp.tsv", sep=""))
# 
# for (son in sonnei_names){
#   system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher SegmentLBI -burnin 0 -calculateAverageDifference true  -outTable combined/", son,"_equal.table.tsv combined/", son,"_equal.network.trees combined/", son,"_equal.mapped.trees", sep=""))
#   system(intern=T, paste("/Applications/BEAST\\ 2.7.6/bin/treeannotator -height mean -burnin 0 combined/", son,"_equal.mapped.trees combined/", son,"_equal.mapped.tree", sep=""))
#   system(intern=F, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ClusterSizeComparison -burnin 0 combined/", son,"_equal.network.trees combined/", son,"_comp.tsv", sep=""))
# }

# compute the 
moving_avg = data.frame()
mean_lbi = data.frame()
time_points <- seq(0, 20, 0.25)
# for (tau in tau_vals){
#   print(tau)
  # t <- read.table(paste("combined/SonFlex_", tau, "_table.tsv", sep=""), header=TRUE, sep="\t")
#   t$condition = t$seg1.LBI
#   moving_avg = movingavg_data_frame(t[t$Time<20,], time_points,5, moving_avg, "either species", paste("LBI with tau=", tau, sep=""))
#   mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, "either species", paste("LBI with tau=", tau, sep=""))
#   t$condition = t$seg1.state0.LBI
#   moving_avg = movingavg_data_frame(t[t$Time<20,], time_points,5, moving_avg, "S. sonnei", paste("LBI with tau=", tau, sep=""))
#   mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, "S. sonnei", paste("LBI with tau=", tau, sep=""))
#   t$condition = t$seg1.state1.LBI
#   moving_avg = movingavg_data_frame(t[t$Time<20,], time_points,5, moving_avg, "S. flexneri", paste("LBI with tau=", tau, sep=""))
#   mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, "S. flexneri", paste("LBI with tau=", tau, sep=""))
# }
  
t <- read.table(paste("combined/SonFlex_", tau_vals[[1]], "_table.tsv", sep=""), header=TRUE, sep="\t")
  
t$condition = t$seg1.CS
# remove any rows that have a NaN somewhere
t = t[complete.cases(t),]
moving_avg = movingavg_data_frame(t[t$Time<20,], time_points,5, moving_avg, "either species", "relative number of offsprings")
mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, "either species", "relative number of offsprings")
t$condition = t$seg1.state0.CS
moving_avg = movingavg_data_frame(t[t$Time<20,], time_points,5, moving_avg, "S. sonnei", "relative number of offsprings")
mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, "S. sonnei", "relative number of offsprings")
t$condition = t$seg1.state1.CS
moving_avg = movingavg_data_frame(t[t$Time<20,], time_points,5, moving_avg, "S. flexneri", "relative number of offsprings")
mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, "S. flexneri", "relative number of offsprings")



colors = c("either species"="#E69F00",
           'S. sonnei'="#44AA99",
           'S. flexneri'="#999933")


p = ggplot(moving_avg, aes(x=as.Date("2020-12-01")-365*time, y=median, color=name, fill=name))+
  geom_line()+
  # geom_point(aes(x=as.Date("2020-12-01")-365*time,y=3.15), color=NA)+
  geom_ribbon(aes(ymin=lower.5, ymax=upper.5), color=NA, alpha=1)+
  geom_ribbon(aes(ymin=lower, ymax=upper), color=NA, alpha=0.35)+
  theme_minimal()+
  scale_fill_manual(name="", values=colors)+
  scale_color_manual(name="", values=colors, guide="none")+
  ylab("average number of\noffpring on lineages with\npKSR100 over compared to without")+
  xlab("")+
  theme(legend.position="none")+
  # add line at y=1
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  guides(fill=guide_legend(nrow=1)) +
  scale_x_date(limits = c(as.Date("2000-01-01"), as.Date("2020-01-01")), date_labels="%Y")
# Plot the graph
plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/SonFlex_rel_lbi_moving.pdf", sep=""),
       width=5, height=3)

# make the x-axis of the plot into x values instead of factors, based on name2 and name
# name order is c("either species", "S. sonnei", "S. flexneri"), name2 order
# is c( paste("LBI with tau=", tau, sep=""), "relative number of offsprings").
mean_lbi$x = which(c("either species", "S. sonnei", "S. flexneri") %in% mean_lbi$name)
# for each value in mean_lbi, check which index in c(paste("LBI with tau=", tau, sep=""), "relative number of offsprings")
# has the same value
comparison = c("relative number of offsprings")
for (i in 1:nrow(mean_lbi)){
  mean_lbi$x[i] = mean_lbi$x[i]+ 0.075*which(comparison %in% mean_lbi$name2[i])
}


p = ggplot(mean_lbi, aes(x=x, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  theme_minimal()+
  ylab("Ratio with\nover without pKSR100")+
  xlab("")+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  # scale_color_manual(values=colors, name="")+
  scale_x_continuous(breaks=seq(1,3)+(length(comparison)+1)*0.075/2, 
                   labels=c("either species", "S. sonnei", "S. flexneri")) +
  theme(panel.grid.minor.x = element_blank()) +
  # scale_color_manual(values=c("#440154", "#3a528b", "#20908c", "#5ec961", "#fde724", "#ed8467"),
  #                    labels=c(paste("LBI: tau=",tau_vals, sep=""), "relative number\n of offsprings" ),
  #                    name="")+
  # scale_y_log10(breaks=c(1,2,4)) +
  theme(legend.position="top")

plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/SonFlex_rel_lbi.pdf", sep=""),
       width=4.5, height=3)

t <- read.table("combined/SonFlex_comp.tsv", header=TRUE, sep="\t")

t$sum = t$clusterSize1+t$clusterSize2

# estimate the weighted
mean_cluster_size = data.frame()
c=0
for (clade in c("any", 0, 1)){
  c2=-0.15-0.075/2
  for (years in c(20,40)){
    if (clade=="any"){
      subset = t[t$time < years,]
    } else {
      subset = t[t$clade==clade & t$time < years,]
    }
    
    for (offspring in c(1,4,9)){
       vals=c()
       for (j in 1:max(subset$iteration)){
        mean_vals = mean(subset[subset$iteration==j & subset$sum>offspring, "clusterSize1"]/(subset[subset$iteration==j & subset$sum>offspring, "sum"]))
        # mean_vals = mean(subset[subset$iteration==j & subset$sum>offspring, "clusterSize1"]>subset[subset$iteration==j & subset$sum>offspring, "clusterSize2"])
        vals = c(vals, mean_vals)
       }
       mean_cluster_size = rbind(mean_cluster_size, data.frame(x=c+c2,
                                                               mean=mean(vals),
                                                               lower = quantile(vals, 0.025), 
                                                               upper = quantile(vals, 0.975),
                                                               years = paste("last",years, "years"),
                                                               offspring = paste("at least",offspring+1),
                                                               clade=clade,
                                                               method=paste("more than ", offspring, "offsprings\nin the last", years, "years")))
       c2=c2+0.075
    }
  }
  c=c+1
}

mean_cluster_size$offspring = factor(mean_cluster_size$offspring, 
                                     levels=c("at least 2", "at least 5", "at least 10"))

mean_cluster_size$species = factor(mean_cluster_size$clade, levels=c("any", 0, 1), labels=c("either species", "S. sonnei", "S. flexneri"))

p = ggplot(mean_cluster_size, aes(x=x, y=mean, color=years, shape=offspring))+
  geom_point(aes(), size=3)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  theme_minimal()+
  ylab("Proportion of offsprings\n on child node with pKSR100")+
  xlab("Clade")+
  scale_color_manual(values=c("#440154", "#20908c"),name="time window")+
  scale_shape(name="total number of\noffspring")+
  # coord_cartesian(ylim=c(0.3,0.9))+
  geom_hline(yintercept=0.5, linetype="dashed", color="black")+
  xlab("")+
  theme(legend.position = "top")+
  guides(color=guide_legend(nrow=3), shape=guide_legend(nrow=3)) +
  scale_x_continuous(breaks=seq(0,2), labels=labels(colors))


plot(p)
ggsave(plot=p,paste("../../../Plasmids-Text/Figures/SonFlex_cluster_size_comp.pdf", sep=""),
       width=6, height=4.5)



fdas
# do the same, but for for Sonnei and spa spb and spc
segmentNames = c("pINV", "spA", "spB", "spC")

colors = c("without plasmid" = "#595959",
           "pINV" = "#CC79A7",
           "spA" = "#0072B2",
           "spB" = "#009E73",
           "spC" = "#D55E00")

for (son in sonnei_names){
  # compute the 
  moving_avg = data.frame()
  mean_lbi = data.frame()
  time_points <- seq(0, 20, 0.25)
  for (tau in tau_vals){
    print(tau)
    t <- read.table(paste("combined/", son, "_", tau, "_table.tsv", sep=""), header=TRUE, sep="\t")
    
    for (j in seq(1,length(segmentNames))){
      t$condition = t[[paste("seg",j,".LBI", sep="")]]
      # remove all rows that have a NaN somewhere
      t = t[complete.cases(t),]
      moving_avg = movingavg_data_frame(t[t$Time<25,], time_points,5, moving_avg, segmentNames[[j]], paste("LBI with tau=", tau, sep=""))
      mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, segmentNames[[j]], paste("LBI with tau=", tau, sep=""))
    }
  }
  for (j in seq(1,length(segmentNames))){
    t$condition = t[[paste("seg",j,".CS", sep="")]]
    moving_avg = movingavg_data_frame(t[t$Time<25,], time_points,5, moving_avg, segmentNames[[j]], "relative number of offsprings")
    mean_lbi = avg_data_frame(t[t$Time<max(time_points),], mean_lbi, segmentNames[[j]], "relative number of offsprings")
  }
  
  p = ggplot(moving_avg[moving_avg$name2=="relative number of offsprings", ], 
             aes(x=as.Date("2020-12-01")-365*time, y=median, color=name, fill=name))+
    geom_line()+
    geom_point(aes(x=as.Date("2020-12-01")-365*time,y=2.5), color=NA)+
    geom_ribbon(aes(ymin=lower.5, ymax=upper.5), color=NA, alpha=1)+
    geom_ribbon(aes(ymin=lower, ymax=upper), color=NA, alpha=0.35)+
    theme_minimal()+
    scale_fill_manual(name="", values=colors)+
    scale_color_manual(name="", values=colors, guide="none")+
    ylab("with over without pKSR100")+
    xlab("")+
    theme(legend.position="top")+
    facet_wrap(name2~., scales="free_y")+
    scale_y_log10()+
    # add line at y=1
    geom_hline(yintercept=1, linetype="dashed", color="black")+
    guides(fill=guide_legend(nrow=1)) +
    scale_x_date(limits = c(as.Date("2000-01-01"), as.Date("2020-01-01")), date_labels="%Y")
  plot(p)
  ggsave(plot=p,paste("../../../Plasmids-Text/Figures/", son, "_rel_lbi_moving.pdf", sep=""),
         width=8, height=4)
  
  
  # make the x-axis of the plot into x values instead of factors, based on name2 and name
  # name order is c("either species", "S. sonnei", "S. flexneri"), name2 order
  # is c( paste("LBI with tau=", tau, sep=""), "relative number of offsprings").
  mean_lbi$x = which(segmentNames %in% mean_lbi$name)
  # for each value in mean_lbi, check which index in c(paste("LBI with tau=", tau, sep=""), "relative number of offsprings")
  # has the same value
  comparison = c(paste("LBI with tau=", tau_vals, sep=""), "relative number of offsprings")
  for (i in 1:nrow(mean_lbi)){
    mean_lbi$x[i] = mean_lbi$x[i]+ 0.075*which(comparison %in% mean_lbi$name2[i])
  }
  
  p = ggplot(mean_lbi, aes(x=x, y=mean, color=name2))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
    theme_minimal()+
    ylab("Ratio with\nover without pKSR100")+
    xlab("")+
    geom_hline(yintercept=1, linetype="dashed", color="black")+
    # scale_color_manual(values=colors, name="")+
    scale_x_continuous(breaks=seq(1,4)+(length(comparison)+1)*0.075/2, 
                       labels=segmentNames[1:4]) +
    theme(panel.grid.minor.x = element_blank()) +
    scale_color_manual(values=c("#440154", "#3a528b", "#20908c", "#5ec961", "#fde724", "#ed8467"),
                       labels=c(paste("LBI: tau=",tau_vals, sep=""), "relative number\n of offsprings" ),
                       name="")+
    scale_y_log10(breaks=c(1,2,4)) +
    theme(legend.position="top")
  plot(p)
  
  ggsave(plot=p,paste("../../../Plasmids-Text/Figures/", son, "_rel_lbi.pdf", sep=""),width=4.5, height=3)
  
  jk
  
  t <- read.table(paste("combined/", son, "_comp.tsv", sep=""), header=TRUE, sep="\t")
  
  t$sum = t$clusterSize1+t$clusterSize2
  
  # estimat the weighted
  mean_cluster_size = data.frame()
  c=0
  for (segment in seq(1,4)){
    c2=-0.15-0.075/2
    for (years in c(20,40)){
      subset = t[t$segment==segment & t$time < years,]
      for (offspring in c(1,4,9)){
        vals=c()
        for (j in 1:max(subset$iteration)){
          mean_vals = mean(subset[subset$iteration==j & subset$sum>offspring, "clusterSize1"]/(subset[subset$iteration==j & subset$sum>offspring, "sum"]))
          # mean_vals = mean(subset[subset$iteration==j & subset$sum>offspring, "clusterSize1"]>subset[subset$iteration==j & subset$sum>offspring, "clusterSize2"])
          vals = c(vals, mean_vals)
        }
        mean_cluster_size = rbind(mean_cluster_size, data.frame(x=c+c2,
                                                                mean=mean(vals),
                                                                lower = quantile(vals, 0.025), 
                                                                upper = quantile(vals, 0.975),
                                                                years = paste("last",years, "years"),
                                                                offspring = paste("at least",offspring+1),
                                                                segment=segmentNames[segment],
                                                                method=paste("more than ", offspring, "offsprings\nin the last", years, "years")))
        c2=c2+0.075
      }
    }
    c=c+1
  }
  
  mean_cluster_size$offspring = factor(mean_cluster_size$offspring, 
                                       levels=c("at least 2", "at least 5", "at least 10"))
  
  p = ggplot(mean_cluster_size, aes(x=x, y=mean, color=years, shape=offspring))+
    geom_point(aes(), size=3)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
    theme_minimal()+
    ylab("Proportion of offsprings\n on child node with pKSR100")+
    xlab("Clade")+
    scale_color_manual(values=c("#440154", "#20908c"),name="time window")+
    scale_shape(name="total number of\noffspring")+
    # coord_cartesian(ylim=c(0.3,0.9))+
    geom_hline(yintercept=0.5, linetype="dashed", color="black")+
    xlab("")+
    theme(legend.position = "top")+
    guides(color=guide_legend(nrow=3), shape=guide_legend(nrow=3)) +
    scale_x_continuous(breaks=seq(0,3)+0.2, labels=labels(colors[2:5]))
  
  
  plot(p)
  ggsave(plot=p,paste("../../../Plasmids-Text/Figures/", son, "_cluster_size_comp.pdf", sep=""),
         width=7, height=4.5)
}
