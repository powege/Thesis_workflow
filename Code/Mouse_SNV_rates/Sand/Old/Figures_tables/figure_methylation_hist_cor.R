rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

### IMPORT 
meth <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz")

### FORMAT

NPC_coverage <- nrow(meth[NPC_coverage != 0]) / nrow(meth)
ES_coverage <- nrow(meth[ES_coverage != 0]) / nrow(meth)
total_coverage <- nrow(meth[ES_coverage != 0 | NPC_coverage != 0]) / nrow(meth)

NPC_coverage <- nrow(meth[NPC_coverage > 5]) / nrow(meth)
ES_coverage <- nrow(meth[ES_coverage > 5]) / nrow(meth)
total_coverage <- nrow(meth[coverage > 5]) / nrow(meth)

cor.test(meth$NPC_percent_methylated, meth$ES_percent_methylated)
cor.test(meth$NPC_percent_methylated[meth$NPC_coverage != 0 & meth$ES_coverage != 0], 
         meth$ES_percent_methylated[meth$ES_coverage != 0 & meth$NPC_coverage != 0])
cor.test(meth$NPC_percent_methylated[meth$NPC_coverage > 5 & meth$ES_coverage > 5], 
         meth$ES_percent_methylated[meth$ES_coverage > 5 & meth$NPC_coverage > 5])

hist(meth$NPC_coverage[meth$NPC_coverage < 100])
hist(meth$ES_coverage[meth$ES_coverage < 100])
hist(meth$coverage[meth$coverage < 100])

hist(meth$NPC_percent_methylated[meth$NPC_coverage != 0])
hist(meth$ES_percent_methylated[meth$ES_coverage != 0])

hist(meth$NPC_percent_methylated[meth$NPC_coverage > 5])
hist(meth$ES_percent_methylated[meth$ES_coverage > 5])
hist(meth$percent_methylated[meth$coverage > 5])


dt_hist <- data.table(percent_methylated = meth$percent_methylated[meth$coverage > 5])
fig_hist <- ggplot(dt_hist, aes(x=percent_methylated)) + 
  geom_histogram(color="black", fill="white")
fig_hist

dt_bar <- data.table(percent_methylated = meth$percent_methylated[meth$coverage > 5])
dt_bar <- as.data.table(table(round(dt_bar$percent_methylated/5)*5))
dt_bar$percentage <- (dt_bar$N/sum(dt_bar$N))*100
dt_bar$V1 <- factor(dt_bar$V1, levels = dt_bar$V1)
fig_bar <- ggplot(data=dt_bar, aes(x=V1, y=percentage)) + 
  geom_bar(stat="identity", color="black", fill="black") +
  geom_rect(aes(xmin = 0.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "transparent", color = "red", size = 1) +
  geom_rect(aes(xmin = 12.5, xmax = 21.5, ymin = -Inf, ymax = Inf),
            fill = "transparent", color = "blue", size = 1) +
  xlab("percentage methylated") +
  ylab("percentage of CGs") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
fig_bar
  

 





