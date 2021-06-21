### SCRIPT that QCs annotation model sequences

rm(list = ls())
graphics.off()

library(data.table)

### FUNCTION
alakazam <- function(mask_file,
                     CG_file,
                     SNV_file,
                     length_min,
                     length_max,
                     sm_max,
                     Nm_max,
                     cm_max,
                     interest){
  
### IMPORT

mask_dt <- fread(mask_file)
CG_dt <- fread(CG_file)
SNV_dt <- fread(SNV_file)

### FORMAT 

# colnames
colnames(mask_dt) <- c("chromosome", "start", "end", "annotation", "strand", "ID", "n_sm", "n_Nm", "n_cm")
colnames(CG_dt) <- c("chromosome", "start", "end", "annotation", "strand", "ID", "n_CG")
colnames(SNV_dt) <- c("chromosome", "start", "end", "annotation", "strand", "ID", "n_SNV")

# subset
mask_dt <- mask_dt[annotation %in% interest]
CG_dt <- CG_dt[annotation %in% interest]
SNV_dt <- SNV_dt[annotation %in% interest]

# merge
dt <- SNV_dt[mask_dt, on = c("chromosome", "start", "end", "annotation", "strand", "ID")]
dt <- dt[CG_dt, on = c("chromosome", "start", "end", "annotation", "strand", "ID")]

# annotation length
dt$length <- (dt$end + 1) - dt$start

# fractions
dt$f_sm <- dt$n_sm / dt$length
dt$f_Nm <- dt$n_Nm / dt$length
dt$f_cm <- dt$n_cm / dt$length

### QC

removed <- data.frame() # set dt for removed 

# filter by length
hist(dt$length, breaks = 100)
rm.id <- which(dt$length < length_min | dt$length > length_max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by N mask fraction
hist(dt$f_Nm)
rm.id <- which(dt$f_Nm > Nm_max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by soft mask fraction
hist(dt$f_sm)
rm.id <- which(dt$f_sm > sm_max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by coverage mask fraction
hist(dt$f_cm)
rm.id <- which(dt$f_cm > cm_max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

dt <- dt[,c("chromosome", "start", "end", "annotation", "strand", "ID", "length", "n_SNV", "n_sm", "n_Nm", "n_cm", "n_CG")]
out <- list(dt, nrow(removed))
return(out)
}

### RUN

M_QCed <- alakazam(
  mask_file = "~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_nMASK.csv",
  CG_file = "~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_nCG.csv",
  SNV_file = "~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_nSNV.csv",
  length_min = 400,
  length_max = 10000,
  sm_max = 0.9,
  Nm_max = 0.05,
  cm_max = 0.5,
  interest = c("Exon - UTR", # annotations of interest
    "Exon - other",
    "Promoter",
    "Enhancer - proximal",
    "Enhancer - distal",
    "CTCF binding")
)
table(M_QCed[[1]]$annotation)

H_QCed <- alakazam(
  mask_file = "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_nMASK.csv",
  CG_file = "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_nCG.csv",
  SNV_file = "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_nSNV.csv",
  length_min = 400,
  length_max = 10000,
  sm_max = 0.9,
  Nm_max = 0.05,
  cm_max = 0.5,
  interest = c("Exon - UTR", # annotations of interest
    "Exon - other",
    "Promoter",
    "Enhancer - proximal",
    "Enhancer - distal",
    "CTCF binding")
)
table(H_QCed[[1]]$annotation)

### EXPORT
fwrite(M_QCed[[1]], "~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_varQCed.csv")
fwrite(H_QCed[[1]], "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_varQCed.csv")


#####

# H_dt <- H_QCed[[1]]
# table(H_dt$annotation)
# 
# H_promoter <- H_dt[annotation == "Promoter"]
# hist(H_promoter$length)
# hist(H_promoter$n_CG)
# 
# mod <- lm(n_SNV ~ length + n_CG, data = H_promoter)
# summary(mod)
# plot(mod)

