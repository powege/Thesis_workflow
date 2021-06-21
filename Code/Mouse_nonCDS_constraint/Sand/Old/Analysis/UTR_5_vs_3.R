rm(list = ls())
graphics.off()

library(data.table)
library(reshape2)


cs <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_CS.csv")
utr <- fread("~/Dropbox/PhD/Data/Target_gene/Human_UTR_tranID.csv")
can <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv")
colnames(utr) <- c("chromosome", "ensembl_transcript_id", "transcript_biotype", "prime", "annotation", "start", "end")
can <- unique(can$ensembl_transcript_id)
dt <- cs[utr, on = c("chromosome", "start", "end", "annotation")]
dt <- dt[complete.cases(dt),]
dt <- unique(dt)
dt <- dt[ensembl_transcript_id %in% can]
dt <- dt[order(ensembl_transcript_id),]

hist(dt$residual[dt$prime == "five_prime_utr"])
hist(dt$residual[dt$prime == "three_prime_utr"])
t.test(dt$residual[dt$prime == "five_prime_utr"], dt$residual[dt$prime == "three_prime_utr"])

dupes <- dt[ensembl_transcript_id %in% dt$ensembl_transcript_id[duplicated(dt$ensembl_transcript_id)]]
dupes <- dupes[order(ensembl_transcript_id),]

dupes_wide <- dupes[,c("ensembl_transcript_id", "prime", "residual")]
dupes_wide <- reshape(dupes_wide, idvar = "ensembl_transcript_id", 
                       timevar = "prime", direction = "wide")
hist(dupes_wide$residual.five_prime_utr)
hist(dupes_wide$residual.three_prime_utr)
t.test(dupes_wide$residual.five_prime_utr, dupes_wide$residual.three_prime_utr, paired = T)

