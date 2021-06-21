rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plyr)

### SET VARS
in.gwas.null <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_NULL_by_chr.csv"
out.gwas.null.figure <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Figures/Figure_GWAS_DDS_null_sample_size_SD.jpeg"

### IMPORT
gwas_null <- fread(in.gwas.null)

out_list_A <- list()
for(j in 1:length(unique(gwas_null$annotation))){
  sub <- gwas_null[annotation == unique(gwas_null$annotation)[j]]
  out_list_B <- list()
  for (i in 1:length(unique(sub$chromosome))){
    sub_sub <- sub[chromosome == unique(sub$chromosome)[i]]
    sub_sub$iteration <- 1:nrow(sub_sub)
    out_list_B[[i]] <- sub_sub
  }
  out_list_A[[j]] <- do.call("rbind", out_list_B)
}
gwas_null <- do.call("rbind", out_list_A)
gwas_null <- gwas_null[iteration <= 400]

out_list_A <- list()
for(j in 1:length(unique(gwas_null$annotation))){
  sub <- gwas_null[annotation == unique(gwas_null$annotation)[j]]
  out_list_B <- list()
  for (i in 1:length(unique(sub$iteration))){
    sub_sub <- sub[iteration == unique(sub$iteration)[i]]
    out_list_B[[i]] <- c(annotation = sub_sub$annotation[1],
                         iteration = sub_sub$iteration[1],
                         apply(sub_sub[,c("n_total", "n_aligned", "n_conserved")], 2, sum))
  }
  out_list_A[[j]] <- do.call("rbind", out_list_B)
}
gwas_null <- as.data.table(do.call("rbind", out_list_A))

gwas_null[, c(2:5)] <- lapply(gwas_null[, c(2:5)], as.numeric)
gwas_null$percentage_aligned_null <- (gwas_null$n_aligned / gwas_null$n_total)*100
gwas_null$percentage_conserved_null <- (gwas_null$n_conserved / gwas_null$n_total)*100

# sub <- gwas_null[annotation == "Exon - CDS"]
catapilar <- function(sub){
sub$percentage_aligned_iteration_mean <- NA
sub$percentage_aligned_iteration_SD <- NA
sub$percentage_conserved_iteration_mean <- NA
sub$percentage_conserved_iteration_SD <- NA
for(i in 1:nrow(sub)){
  sub$percentage_aligned_iteration_mean[i] <- mean(sub$percentage_aligned_null[1:i])
  sub$percentage_aligned_iteration_SD[i] <- sd(sub$percentage_aligned_null[1:i])
  sub$percentage_conserved_iteration_mean[i] <- mean(sub$percentage_conserved_null[1:i])
  sub$percentage_conserved_iteration_SD[i] <- sd(sub$percentage_conserved_null[1:i])
}
return(sub)
}
gwas_null <- ddply(gwas_null, "annotation", catapilar)

# sub <- gwas_null[gwas_null$annotation == "Exon - CDS",]
beetle <- function(sub){
  sub$percentage_aligned_iteration_mean_r <- sub$percentage_aligned_iteration_mean / sub$percentage_aligned_iteration_mean[nrow(sub)]
  sub$percentage_aligned_iteration_SD_r <- sub$percentage_aligned_iteration_SD / sub$percentage_aligned_iteration_SD[nrow(sub)]
  sub$percentage_conserved_iteration_mean_r <- sub$percentage_conserved_iteration_mean / sub$percentage_conserved_iteration_mean[nrow(sub)]
  sub$percentage_conserved_iteration_SD_r <- sub$percentage_conserved_iteration_SD / sub$percentage_conserved_iteration_SD[nrow(sub)]
  return(sub)
}
gwas_null <- ddply(gwas_null, "annotation", beetle)



### PLOT

plot_function <- function(dt, y, xlab, ylab, title){
  
out <- ggplot(dt, aes(x=iteration, y=y, color=annotation)) +
  geom_line() +
  xlab(xlab) +
  ylab(ylab) +
  ggtitle(title) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )

return(out)
}

plot_function(dt = gwas_null, 
              y = gwas_null$percentage_aligned_iteration_mean_r, 
              xlab = "Null sample size (n)", 
              ylab = "Null sample mean",
              title = "Alignment")
plot_function(dt = gwas_null, 
              y = gwas_null$percentage_aligned_iteration_SD_r, 
              xlab = "Null sample size (n)", 
              ylab = "Null sample standard deviation",
              title = "Alignment")
plot_function(dt = gwas_null, 
              y = gwas_null$percentage_conserved_iteration_mean_r, 
              xlab = "Null sample size (n)", 
              ylab = "Null sample mean",
              title = "Conservation")
plot_function(dt = gwas_null, 
              y = gwas_null$percentage_conserved_iteration_SD_r, 
              xlab = "Null sample size (n)", 
              ylab = "Null sample standard deviation",
              title = "Conservation")

#####


# ggplot(gwas_null, aes(x=iteration, y=percentage_aligned_iteration_mean, color=annotation)) +
#   geom_line() +
#   facet_grid(annotation ~ ., scales = "free")
#   
# ggplot(gwas_null, aes(x=iteration, y=percentage_aligned_iteration_SD, color=annotation)) +
#   geom_line() +
#   facet_grid(annotation ~ ., scales = "free")
# 
# ggplot(gwas_null, aes(x=iteration, y=percentage_conserved_iteration_mean, color=annotation)) +
#   geom_line() +
#   facet_grid(annotation ~ ., scales = "free")
# 
# ggplot(gwas_null, aes(x=iteration, y=percentage_conserved_iteration_SD, color=annotation)) +
#   geom_line() +
#   facet_grid(annotation ~ ., scales = "free")



