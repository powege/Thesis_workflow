rm(list = ls())
graphics.off()

library(data.table)

### SET ARGS
ann.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
cs.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/mmus_GRC38_constraint_by_window_1000_100_WMallSP.csv.gz"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/mmus_annotation_by_constraint_percentile_750_50_WMallSP.csv.gz"
functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"

### FUNCTIONS
source(functions.file)

### IMPORT
ann <- fread(paste0("gunzip -cq ", ann.infile))
cs <- fread(paste0("gunzip -cq ", cs.infile))

### FORMAT
annotation <- unique(ann$category)
out_list <- list()
for (i in 1:100){
  
  percentile_sub <- cs[residual_percentile == i]
  
  n_bp_overlap <- rep(NA, length(annotation))
  for(j in 1:length(annotation)){
    
    ann_sub <- ann[category == annotation[j]]
    ann_sub <- collapse.overlap(ann_sub[,c("chromosome", "start", "end")])
    overlap_sub <- bed.intersect(ann_sub, 
                                 percentile_sub[,c("chromosome", "start", "end")])
    n_bp_overlap[j] <- sum((overlap_sub$end +1) - overlap_sub$start)
  }
  
  out_list[[i]] <- data.table(percentile = rep(i, length(annotation)),
                              annotation = annotation,
                              n_bp_percentile = rep(sum((percentile_sub$end +1) - percentile_sub$start)),
                              n_bp_overlap = n_bp_overlap)
  print(i)
}

output <- do.call("rbind", out_list)
output$fraction <- output$n_bp_overlap / output$n_bp_percentile

### EXPORT
fwrite(output, out.file)

#####
x <- output[annotation == "Exon - CDS"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Exon - 5'UTR"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Exon - 3'UTR"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Exon - other"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Promoter"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Enhancer - proximal"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Enhancer - distal"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "TAD boundry"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Miscellaneous"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "CTCF binding"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Intron - proximal"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Intron - distal"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]

x <- output[annotation == "Unannotated"]
plot(x$percentile, x$fraction)
x$fraction[1] / x$fraction[100]
