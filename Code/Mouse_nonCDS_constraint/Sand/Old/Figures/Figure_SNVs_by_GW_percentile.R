rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)

### FUNCTIONS

gobbler <- function(promoter){
  percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                        "51-75", "76-90", "91-95", "96-98", "99", "100")
  mean_fold_change <- rep(NA, 12)
  mean_fold_change[1] <- mean(promoter$fold_change[promoter$percentile %in% c(1)])
  mean_fold_change[2] <- mean(promoter$fold_change[promoter$percentile %in% c(2)])
  mean_fold_change[3] <- mean(promoter$fold_change[promoter$percentile %in% c(3:5)])
  mean_fold_change[4] <- mean(promoter$fold_change[promoter$percentile %in% c(6:10)])
  mean_fold_change[5] <- mean(promoter$fold_change[promoter$percentile %in% c(11:25)])
  mean_fold_change[6] <- mean(promoter$fold_change[promoter$percentile %in% c(26:50)])
  mean_fold_change[7] <- mean(promoter$fold_change[promoter$percentile %in% c(51:75)])
  mean_fold_change[8] <- mean(promoter$fold_change[promoter$percentile %in% c(76:90)])
  mean_fold_change[9] <- mean(promoter$fold_change[promoter$percentile %in% c(91:95)])
  mean_fold_change[10] <- mean(promoter$fold_change[promoter$percentile %in% c(96:98)])
  mean_fold_change[11] <- mean(promoter$fold_change[promoter$percentile %in% c(99)])
  mean_fold_change[12] <- mean(promoter$fold_change[promoter$percentile %in% c(100)])
  out <- data.frame(percentile_group = percentile_group,
                    rank = c(1:12),
                    mean_fold_change = mean_fold_change,
                    annotation = rep(promoter$category[1], 12))
  return(out)
}

gobbler2 <- function(promoter){
  percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                        "51-75", "76-90", "91-95", "96-98", "99", "100")
  mean_fraction <- rep(NA, 12)
  mean_fraction[1] <- mean(promoter$fraction[promoter$percentile %in% c(1)])
  mean_fraction[2] <- mean(promoter$fraction[promoter$percentile %in% c(2)])
  mean_fraction[3] <- mean(promoter$fraction[promoter$percentile %in% c(3:5)])
  mean_fraction[4] <- mean(promoter$fraction[promoter$percentile %in% c(6:10)])
  mean_fraction[5] <- mean(promoter$fraction[promoter$percentile %in% c(11:25)])
  mean_fraction[6] <- mean(promoter$fraction[promoter$percentile %in% c(26:50)])
  mean_fraction[7] <- mean(promoter$fraction[promoter$percentile %in% c(51:75)])
  mean_fraction[8] <- mean(promoter$fraction[promoter$percentile %in% c(76:90)])
  mean_fraction[9] <- mean(promoter$fraction[promoter$percentile %in% c(91:95)])
  mean_fraction[10] <- mean(promoter$fraction[promoter$percentile %in% c(96:98)])
  mean_fraction[11] <- mean(promoter$fraction[promoter$percentile %in% c(99)])
  mean_fraction[12] <- mean(promoter$fraction[promoter$percentile %in% c(100)])
  out <- data.frame(percentile_group = percentile_group,
                    rank = c(1:12),
                    mean_fraction = mean_fraction,
                    annotation = rep(promoter$annotation[1], 12))
  return(out)
}

# FUNCTION to calculate fold change vs the 100th percentile
fold.change <- function(sub){
  sub$fold_change <- sub$fraction/sub$fraction[sub$percentile == 100]
  return(sub)
}

### SET ARGS

### IMPORT
m_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/ClinVar_SNV_by_OE_percentile_750_50_MGP_allMUSMUS.csv")

### FORMAT

# calculate fold change vs the 100th percentile
m_fc <- ddply(m_dt, "category", fold.change)
# group percentiles with custom function
m_fc <- ddply(m_fc, "category", gobbler)
# set levels order
levels(m_fc$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                                   "51-75", "76-90", "91-95", "96-98", "99", "100")


### PLOT FOLD CHANGE

m_p_fc <- ggplot(m_fc, aes(x=rank, y=mean_fold_change, color=category)) +
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Pathogenic SNV sites\n(fold change versus 100th percentile)") +
  # scale_color_brewer(palette="Set3") +
  scale_y_continuous(
    trans = log2_trans(),
    # breaks = trans_breaks("log2", function(x) 2^x)) +
    breaks = c(0.5, 1, 2, 5, 10, 20),
    limits = c(0.5, 35)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
m_p_fc
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_SNVs_by_GW_percentile_line.jpg", plot = m_p_fc, height = 5, width = 7)


p <- ggplot(m_dt, aes(x=percentile, y=fraction, color=category)) +
  geom_point(alpha = 1/3) +
  # geom_line() 
  geom_smooth()
p
