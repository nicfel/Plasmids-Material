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
library(ggpattern)  # For pattern fills

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

names = c("chromosome", "pINV", "spA", "spB", "spC")


colors = c("spA"="#0072B2",
            "entire spA"="#0072B2",
           "strAB + sul\n+ flanking"="#0072B2",
           "AMR genes\nonly"="#0072B2",
           'spB'="#009E73",
           'spC'="#D55E00",
           'pINV'="#E69F00")


log_dta = read.table("./supplementalxmls/dta.log", header=TRUE, sep="\t")
# remove first 10%
log_dta = log_dta[seq(ceiling(0.1*nrow(log_dta)), nrow(log_dta)), ]

# for all plasmids pInV, spA, spB, spC, get traitClockRate.pInv and traitClockRate.pInv_rand on the x axis and their values on the y-axis
# Extract trait clock rates for all plasmids
plasmids <- c("pInv", "spA", "spB", "spC")

plot_data <- data.frame()

c = 1
for (plasmid in plasmids) {
  for (isrand in c(F, T)){
      # Get columns for this plasmid's trait clock rates
    if (!isrand){
      trait_col <- paste0("traitClockRate.", plasmid)
    }else{
      trait_col <- paste0("traitClockRate.", plasmid, "_rand")
    }
    
    plot_data <- rbind(plot_data, data.frame(x=c, y=log_dta[,trait_col], plasmid=plasmid, isrand=isrand))
    c=c+1
  }
}

# Fix color mapping - use uppercase pINV for color mapping
plot_data$plasmid_color <- ifelse(plot_data$plasmid == "pInv", "pINV", plot_data$plasmid)

# Create custom labels for x-axis - only show plasmid names at appropriate positions
x_labels <- rep("", max(plot_data$x))
plasmid_positions <- aggregate(x ~ plasmid, plot_data, mean)
for(i in 1:nrow(plasmid_positions)) {
  x_labels[round(plasmid_positions$x[i])] <- plasmid_positions$plasmid[i]
}

# Create the plot with pattern fill
p1 = ggplot(plot_data) +
  geom_violin_pattern(aes(x=x, y=y, 
                          fill=plasmid_color, 
                          pattern=ifelse(isrand, "stripe", "none")), 
                      color="black", 
                      size=0.3,
                      pattern_density=0.3,
                      pattern_spacing=0.02,
                      pattern_angle=45,
                      pattern_color="white",
                      pattern_size=0.3) +
  scale_fill_manual(values=colors) +
  scale_pattern_manual(values=c("none"="none", "stripe"="stripe"), 
                       name="", 
                       labels=c("none"="Observed Presence/Absence", 
                                "stripe"="Randomized Presence/Absence")) +
  scale_x_continuous(breaks=seq(1.5, by = 2, length.out = length(plasmids)), 
                     labels=plasmids) +
  scale_y_log10() +
  theme_minimal() +
  labs(
    x = "Plasmid",
    y = "Rate of change pre year",
    title = "Rate of change in plasmid presence absence"
  ) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = "none")

print(p1)

ggsave(plot=p1,paste("../../../Plasmids-Text/Figures/Randomized_rates.pdf", sep=""),
       width=6, height=4)

c=1
plot_data <- data.frame()

for (plasmid in plasmids) {
  for (gain in c("gain", "loss")){
    # Get columns for this plasmid's trait clock rates
    trait_col <- paste0("traitClockRate.", plasmid)
    if (grepl("gain", gain)){
      # relativeGeoRates.pInv.1
      rel_rate_col <- paste0("relativeGeoRates.", plasmid, ".1")
    }else{
      # relativeGeoRates.pInv.2
      rel_rate_col <- paste0("relativeGeoRates.", plasmid, ".2")
    }
    
    plot_data <- rbind(plot_data, data.frame(x=c, y=log_dta[,trait_col]*log_dta[,rel_rate_col], 
                                              plasmid=plasmid, isrand=F, gain=gain))
    c=c+1
  }
}

# Fix color mapping - use uppercase pINV for color mapping
plot_data$plasmid_color <- ifelse(plot_data$plasmid == "pInv", "pINV", plot_data$plasmid)
# Create custom labels for x-axis - only show plasmid names at appropriate positions
x_labels <- rep("", max(plot_data$x))
plasmid_positions <- aggregate(x ~ plasmid, plot_data, mean)
for(i in 1:nrow(plasmid_positions)) {
  x_labels[round(plasmid_positions$x[i])] <- plasmid_positions$plasmid[i]
}

plot_data$isgain = ifelse(plot_data$gain == "gain", T, F)

# Create the plot with pattern fill for gain/loss
p2 = ggplot(plot_data[!is.na(plot_data$gain), ]) +
  geom_violin_pattern(aes(x=x, y=y, 
                          fill=plasmid_color, 
                          pattern=ifelse(isgain, "stripe", "none")), 
                      color="black", 
                      size=0.3,
                      pattern_density=0.3,
                      pattern_spacing=0.02,
                      pattern_angle=45,
                      pattern_color="white",
                      pattern_size=0.3) +
  scale_fill_manual(values=colors) +
  scale_pattern_manual(values=c("none"="none", "stripe"="stripe"), 
                       name="", 
                       labels=c("none"="gain", 
                                "stripe"="loss")) +
  scale_x_continuous(breaks=seq(1.5, by = 2, length.out = length(plasmids)), 
                     labels=plasmids) +
  scale_y_log10() +
  theme_minimal() +
  labs(
    x = "Plasmid",
    y = "Rate of change",
    title = "Posterior Estimates of Plasmids being gained or lost"
  ) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = "none")
print(p2)

library(dplyr)

# Per-violin summaries for annotation
anno <- plot_data %>%
  filter(!is.na(gain)) %>%
  group_by(x, plasmid_color, gain) %>%
  summarize(
    med  = median(y, na.rm = TRUE),
    ymax = max(y, na.rm = TRUE),
    .groups = "drop"
  )

# Add median labels above each violin
p2 <- p2 +
  geom_text(
    data = anno,
    aes(x = x, y = ymax * 1.15, label = signif(med, 2)),
    vjust = 0, size = 3
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(y = max(anno$ymax, na.rm = TRUE) * 1.2) +
  theme(plot.margin = margin(t = 12, r = 5, b = 5, l = 5))


p2 <- p2 + stat_summary(fun = median, aes(x = x, y = y), geom = "point", size = 1.4, color = "black")
print(p2)

ggsave(plot=p2,paste("../../../Plasmids-Text/Figures/DTA_gain_loss.pdf", sep=""),
       width=6, height=4)



log_com = read.table("./supplementalxmls/plasmidComp.log", header=TRUE, sep="\t")
# remove first 10%
log_com = log_com[seq(ceiling(0.2*nrow(log_com)), nrow(log_com)), ]
# loop over the row 
reassort = data.frame()
for (i in 1:nrow(log_com)) {
  len = log_com[i, "Tree.spA.treeLength"]
  # multiply by exprnd(0.01)
  valrnd = len*rexp(1, 1/0.01)
  # take a poisson random number with mean valrnd
  val = rpois(1, valrnd)
  reassort = rbind(reassort, data.frame(x=val, value="prior"))
  reassort = rbind(reassort, data.frame(x=log_com[i, "networkCwR.alltrees.reassortmentNodeCount"], value="posterior"))
}

# Create the plot for reassortment
p3 = ggplot(reassort, aes(x=x, fill=value, color=value)) +
  geom_histogram(binwidth=1, alpha=0.5, position="identity") +
  scale_fill_manual(values=c("prior"="#D55E00", "posterior"="#0072B2")) +
  scale_color_manual(values=c("prior"="#D55E00", "posterior"="#0072B2"), guide="none") +
  labs(
    x = "Reassortment Node Count",
    y = "Frequency",
    title = "Reassortment Node Count Distribution"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  guides(fill=guide_legend(title=""))
print(p3)



log1 = read.table("./supplementalxmls/RelaxedClock.log", header=TRUE, sep="\t")
log2 = read.table("./supplementalxmls/RelaxedClockStrict.log", header=TRUE, sep="\t")
log1 = log1[seq(ceiling(0.3*nrow(log1)), nrow(log1)), ]
log2 = log2[seq(ceiling(0.3*nrow(log2)), nrow(log2)), ]

data = data.frame(y=log1[, "ORCRatesStat.coefficientOfVariation"], value="true data")
data = rbind(data, data.frame(y=log2[, "ORCRatesStat.coefficientOfVariation"], value="data simulated using a strict clock model"))

p4 = ggplot(data, aes(x=y, fill=value, color=value)) +
  geom_histogram(binwidth=0.01, alpha=0.5, position="identity") +
  scale_fill_manual(values=c("true data"="#0072B2", "data simulated using a strict clock model"="#D55E00")) +
  scale_color_manual(values=c("true data"="#0072B2", "data simulated using a strict clock model"="#D55E00"), guide="none") +
  labs(
    x = "Coefficient of Variation",
    y = "Frequency",
    title = "Coefficient of Variation for\ntrue and simulated data"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  guides(fill=guide_legend(title=""))
print(p4)

ggsave(plot=p4,paste("../../../Plasmids-Text/Figures/relaxedClock.pdf", sep=""),
       width=6, height=4)

