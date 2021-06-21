rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(tidyr)

### FUNCTION for plotting lm equation on ggplot
linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 2),
            yy = format(abs(coef(k)[2]), digits = 2),
            r2 = format(summary(k)$r.squared, digits = 2));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)   
  }
  as.character(as.expression(eq));               
}

### FUNCTION for plotting r2 and p.val on ggplot
# k <- mod_all
r2p = function(k) {
  eq <- substitute(italic(r)^2~"="~rvalue*","~italic(p)~"="~pvalue, 
                   list(rvalue = format(summary(k)$r.squared, digits = 2), 
                        pvalue = format(anova(k)$'Pr(>F)'[1], digits = 2)))
    as.character(as.expression(eq))            
}

### FUNCTION for plotting Spearman's correlation on ggplot
# k <- cor.test(df_md_all$x, df_md_all$Category, method = "spearman")
rs_corr_eqn <- function(k, digits) {
  output <- substitute(italic(rho)~"="~rhovalue*","~italic(p)~"="~pvalue, 
                   list(rhovalue = unname(format(k$estimate, digits = digits)), 
                        pvalue = unname(format(k$p.value, digits = digits))))
  as.character(as.expression(output))  
  }

### IMPORT DATA
dt_top_level <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_OE_IMPC_pleiotropy_top_level.csv")
dt_procedure <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_OE_IMPC_pleiotropy_procedure.csv")
dt_parameter <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_OE_IMPC_pleiotropy_parameter.csv")


### FORMAT

dt_top_level_avg <- aggregate(list(mean_hit_rate=dt_top_level$hit_rate), by=list(score_percentile=dt_top_level$score_percentile), FUN=median)
dt_procedure_avg <- aggregate(list(mean_hit_rate=dt_procedure$hit_rate), by=list(score_percentile=dt_procedure$score_percentile), FUN=median)
dt_parameter_avg <- aggregate(list(mean_hit_rate=dt_parameter$hit_rate), by=list(score_percentile=dt_parameter$score_percentile), FUN=median)

# lm
summary(lm(dt_top_level_avg$mean_hit_rate~dt_top_level_avg$score_percentile))
summary(lm(dt_procedure_avg$mean_hit_rate~dt_procedure_avg$score_percentile))
summary(lm(dt_parameter_avg$mean_hit_rate~dt_parameter_avg$score_percentile))
cor.test(dt_top_level_avg$mean_hit_rate, dt_top_level_avg$score_percentile, method = "spearman")
cor.test(dt_procedure_avg$mean_hit_rate, dt_procedure_avg$score_percentile, method = "spearman")
cor.test(dt_parameter_avg$mean_hit_rate, dt_parameter_avg$score_percentile, method = "spearman")

equation_top_level <- linear(lm(dt_top_level_avg$mean_hit_rate~dt_top_level_avg$score_percentile))
r2p_top_level <- r2p(lm(dt_top_level_avg$mean_hit_rate~dt_top_level_avg$score_percentile))
rho_top_level <- rs_corr_eqn(cor.test(dt_top_level_avg$mean_hit_rate, dt_top_level_avg$score_percentile, method = "spearman"), 2)

equation_procedure <- linear(lm(dt_procedure_avg$mean_hit_rate~dt_procedure_avg$score_percentile))
r2p_procedure <- r2p(lm(dt_procedure_avg$mean_hit_rate~dt_procedure_avg$score_percentile))
rho_procedure <- rs_corr_eqn(cor.test(dt_procedure_avg$mean_hit_rate, dt_procedure_avg$score_percentile, method = "spearman"), 2)


### PLOT SCATTER
arbok <- function(dt_top_level_avg, rho_top_level, title, colour){
plot <- ggplot(dt_top_level_avg, aes(x = score_percentile, y = mean_hit_rate)) +
  geom_point(alpha = 1/2, colour = colour) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "blue", size = 0.6, fullrange = T) +
  annotate("text", x = 37, y = max(dt_top_level_avg$mean_hit_rate), label = rho_top_level, colour="black", size = 5, parse=TRUE) +
  xlab("NOER percentile rank") +
  ylab('pleiotropy (median hit-rate)') +
  ggtitle(title) +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))
return(plot)
}

p_top_level_scatter <- arbok(dt_top_level_avg, rho_top_level, "(D) Top-level MP terms", "#009E73")
p_procedure_scatter <- arbok(dt_procedure_avg, rho_procedure, "(C) Procedures", "#CC79A7")


### PLOT HIST

ekans <- function(dt_top_level, title, fill){
dt_bar <- data.table(hit_rate = (dt_top_level$hit_rate)*100)
dt_bar <- as.data.table(table(round(dt_bar$hit_rate/10)*10))
dt_bar$V1 <- as.numeric(dt_bar$V1)/100
dt_bar$percentage <- (dt_bar$N/sum(dt_bar$N))*100
dt_bar$V1 <- factor(dt_bar$V1, levels = dt_bar$V1)
fig_bar <- ggplot(data=dt_bar, aes(x=V1, y=percentage)) + 
  geom_bar(stat="identity", color="black", fill=fill) +
  ggtitle(title) +
  xlab("pleiotropy (hit-rate bins)") +
  ylab("percentage of KOs") +
  theme_classic() +
  theme(plot.title =  element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        # axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
return(fig_bar)
}
p_top_level_hist <- ekans(dt_top_level, "(B) Top-level MP terms", "#009E73")
p_procedure_hist <- ekans(dt_procedure, "(A) Procedures", "#CC79A7")

### EXPORT 
pout <- grid.arrange(p_procedure_hist, p_top_level_hist, p_procedure_scatter, p_top_level_scatter, nrow = 2, ncol = 2, widths = c(1, 1))
ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/figure_main_IMPC_pleiotropy_hist_scatter.jpg", plot = pout, height = 9, width = 9)


#######
