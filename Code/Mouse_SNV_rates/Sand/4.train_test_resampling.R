rm(list = ls())
graphics.off()

library(data.table)
library(stringi)
library(plyr)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### SET VARS
k11.p.snv.specific.file <- args[1]
out.specific.AC.file <- args[2]
out.specific.CCG.file <- args[3]

# k11.p.snv.specific.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz"

### IMPORT
p_snv_specific <- fread(paste0("gunzip -cq ", k11.p.snv.specific.file))

### FORMAT

## Split into train and test datasets with 80:20 split. 

# subset A from
p.snv.A <- p_snv_specific[k1_from == "A"]
p.snv.A <- p.snv.A[,c("to", "k11_from", "k11_from_N", "k11_to_N")]
p.snv.A.wide <- dcast(p.snv.A, k11_from + k11_from_N ~ to, value.var="k11_to_N")
p.snv.A.wide$X <- p.snv.A.wide$k11_from_N - (p.snv.A.wide$C + p.snv.A.wide$G + p.snv.A.wide$T)
tmp_list <- split(p.snv.A.wide, rep(1:10, length.out = nrow(p.snv.A.wide), each = ceiling(nrow(p.snv.A.wide)/10))) # split as otherwise vector memory reached
for(i in 1:length(tmp_list)){
  tmp <- tmp_list[[i]]
  tmp_list[[i]] <- dcast(
    tmp[,k11_from:=k11_from
                 ][, .(draw=sample(c(rep("C",C),
                                     rep("G",G),
                                     rep("T",T),
                                     rep("X",X)),
                                   size=round(.8*k11_from_N),
                                   replace = F)  )
                   , by=k11_from
                   ][,.N,
                     by=.(k11_from,draw)],
    k11_from~draw,value.var="N")
  print(i)
}
dt_80_A <- do.call("rbind", tmp_list)
dt_80_A[is.na(dt_80_A)] <- 0
dt_80_A$k11_from_N <- dt_80_A$C + dt_80_A$G + dt_80_A$T + dt_80_A$X
dt_80_A$X <- NULL
p.snv.A.wide$X <- NULL
ind <- which(p.snv.A.wide$k11_from %in% dt_80_A$k11_from)
dt_80_A <- rbind(p.snv.A.wide[-ind,], dt_80_A)
dt_80_A <- dt_80_A[order(k11_from),]
p.snv.A.wide <- p.snv.A.wide[order(k11_from),]
dt_20_A <- data.table(
  k11_from = dt_80_A$k11_from,
  k11_from_N = p.snv.A.wide$k11_from_N - dt_80_A$k11_from_N,
  C = p.snv.A.wide$C - dt_80_A$C,
  G = p.snv.A.wide$G - dt_80_A$G,
  T = p.snv.A.wide$T - dt_80_A$T
)
dt_20_A$k11_from_N <- p.snv.A.wide$k11_from_N - dt_80_A$k11_from_N
dt_20_A$k11_from <- dt_80_A$k11_from
dt_80_A_long <- melt(dt_80_A,
                   id.vars=c("k11_from", "k11_from_N"),
                   measure.vars=c("C", "G", "T" ),
                   variable.name="to",
                   value.name="k11_to_N"
)
dt_80_A_long <- dt_80_A_long[order(k11_from),]
dt_20_A_long <- melt(dt_20_A,
                   id.vars=c("k11_from", "k11_from_N"),
                   measure.vars=c("C", "G", "T" ),
                   variable.name="to",
                   value.name="k11_to_N"
)
dt_20_A_long <- dt_20_A_long[order(k11_from),]


# subset C from
p.snv.C <- p_snv_specific[k1_from == "C"]
p.snv.C <- p.snv.C[,c("to", "k11_from", "k11_from_N", "k11_to_N")]
p.snv.C.wide <- dcast(p.snv.C, k11_from + k11_from_N ~ to, value.var="k11_to_N")
p.snv.C.wide$X <- p.snv.C.wide$k11_from_N - (p.snv.C.wide$A + p.snv.C.wide$G + p.snv.C.wide$T)
tmp_list <- split(p.snv.C.wide, rep(1:10, length.out = nrow(p.snv.C.wide), each = ceiling(nrow(p.snv.C.wide)/10))) # split as otherwise vector memory reached
for(i in 1:length(tmp_list)){
  tmp <- tmp_list[[i]]
  tmp_list[[i]] <- dcast(
    tmp[,k11_from:=k11_from
        ][, .(draw=sample(c(rep("A",A),
                            rep("G",G),
                            rep("T",T),
                            rep("X",X)),
                          size=round(.8*k11_from_N),
                          replace = F)  )
          , by=k11_from
          ][,.N,
            by=.(k11_from,draw)],
    k11_from~draw,value.var="N")
  print(i)
}
dt_80_C <- do.call("rbind", tmp_list)
dt_80_C[is.na(dt_80_C)] <- 0
dt_80_C$k11_from_N <- dt_80_C$A + dt_80_C$G + dt_80_C$T + dt_80_C$X
dt_80_C$X <- NULL
p.snv.C.wide$X <- NULL
ind <- which(p.snv.C.wide$k11_from %in% dt_80_C$k11_from)
dt_80_C <- rbind(p.snv.C.wide[-ind,], dt_80_C)
dt_80_C <- dt_80_C[order(k11_from),]
p.snv.C.wide <- p.snv.C.wide[order(k11_from),]
dt_20_C <- data.table(
  k11_from = dt_80_C$k11_from,
  k11_from_N = p.snv.C.wide$k11_from_N - dt_80_C$k11_from_N,
  A = p.snv.C.wide$A - dt_80_C$A,
  G = p.snv.C.wide$G - dt_80_C$G,
  T = p.snv.C.wide$T - dt_80_C$T
)
dt_20_C$k11_from_N <- p.snv.C.wide$k11_from_N - dt_80_C$k11_from_N
dt_20_C$k11_from <- dt_80_C$k11_from
dt_80_C_long <- melt(dt_80_C,
                     id.vars=c("k11_from", "k11_from_N"),
                     measure.vars=c("A", "G", "T" ),
                     variable.name="to",
                     value.name="k11_to_N"
)
dt_80_C_long <- dt_80_C_long[order(k11_from),]
dt_20_C_long <- melt(dt_20_C,
                     id.vars=c("k11_from", "k11_from_N"),
                     measure.vars=c("A", "G", "T" ),
                     variable.name="to",
                     value.name="k11_to_N"
)
dt_20_C_long <- dt_20_C_long[order(k11_from),]

dt_80 <- rbind(dt_80_A_long, dt_80_C_long)
dt_20 <- rbind(dt_20_A_long, dt_20_C_long)
rm(dt_80_A_long, dt_80_C_long, dt_20_A_long, dt_20_C_long, dt_80_A, dt_80_C, dt_20_A, dt_20_C,
   p.snv.A, p.snv.A.wide, p.snv.C, p.snv.C.wide, tmp_list, tmp, p_snv_specific)

## test mu rate
test_specific <- dt_20
test_specific$k1_from <- stri_sub(test_specific$k11_from, 6, -6)
test_specific$k11_mu_rate <- test_specific$k11_to_N / test_specific$k11_from_N
test_specific$k11_mu_rate[test_specific$k11_from_N == 0] <- 0
rm(dt_20)

## train mu rates
k11_specific <- dt_80 # 11mer

k9_specific <- k11_specific # 9mer
colnames(k9_specific) <- gsub("k11", "k9", colnames(k11_specific)) # colnames
k9_specific$k9_from <- stri_sub(k9_specific$k9_from, 2, -2)
k9_specific <- setDT(k9_specific)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, to)] # sum across kmers

k7_specific <- k11_specific # 7mer
colnames(k7_specific) <- gsub("k11", "k7", colnames(k11_specific)) # colnames
k7_specific$k7_from <- stri_sub(k7_specific$k7_from, 3, -3)
k7_specific <- setDT(k7_specific)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, to)] # sum across kmers

k5_specific <- k11_specific # 5mer
colnames(k5_specific) <- gsub("k11", "k5", colnames(k11_specific)) # colnames
k5_specific$k5_from <- stri_sub(k5_specific$k5_from, 4, -4)
k5_specific <- setDT(k5_specific)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, to)] # sum across kmers

k3_specific <- k11_specific # 3mer
colnames(k3_specific) <- gsub("k11", "k3", colnames(k11_specific))
k3_specific$k3_from <- stri_sub(k3_specific$k3_from, 5, -5)
k3_specific <- setDT(k3_specific)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, to)] # sum across kmers

k1_specific <- k11_specific # 1mer
colnames(k1_specific) <- gsub("k11", "k1", colnames(k11_specific))
k1_specific$k1_from <- stri_sub(k1_specific$k1_from, 6, -6)
k1_specific <- setDT(k1_specific)[, .(k1_from_N = sum(k1_from_N), k1_to_N = sum(k1_to_N)), by = .(k1_from, to)] # sum across kmers

k1CG_1 <- k3_specific[which(stri_sub(k3_specific$k3_from, 2, -2) == "A" | stri_sub(k3_specific$k3_from, 2, -2) == "C")] # k1 CG
k1CG_1$k1CG_from <- NA
k1CG_1$k1CG_from[ unique(c( grep("CG", k1CG_1$k3_from))) ] <- "C (CG)"
k1CG_1$k1CG_from[is.na(k1CG_1$k1CG_from)] <- "C (nonCG)"
k1CG_1$k1CG_from[which(stri_sub(k1CG_1$k3_from, 2, -2) == "A")] <- "A"
if (length(unique(k1CG_1$k3_from_N)) == 32){
  CG_sub <- k1CG_1[k1CG_from == "C (CG)"]
  nonCG_sub <- k1CG_1[k1CG_from == "C (nonCG)"]
  A_sub <- k1CG_1[k1CG_from == "A"]
  k1CG_2 <- data.table(k1CG_from = c(rep("C (CG)", 3), rep("C (nonCG)", 3), rep("A", 3)),
                       to = c(rep(c("A", "G", "T"), 2), rep(c("C", "G", "T"), 1)),
                       k1CG_from_N = c(rep(sum(unique(CG_sub$k3_from_N)), 3),
                                       rep(sum(unique(nonCG_sub$k3_from_N)), 3),
                                       rep(sum(unique(A_sub$k3_from_N)), 3)),
                       k1CG_to_N =c(sum(CG_sub$k3_to_N[CG_sub$to == "A"]),
                                    sum(CG_sub$k3_to_N[CG_sub$to == "G"]),
                                    sum(CG_sub$k3_to_N[CG_sub$to == "T"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$to == "A"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$to == "G"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$to == "T"]),
                                    sum(A_sub$k3_to_N[A_sub$to == "C"]),
                                    sum(A_sub$k3_to_N[A_sub$to == "G"]),
                                    sum(A_sub$k3_to_N[A_sub$to == "T"]))
  )
  rm(CG_sub, nonCG_sub, A_sub)
}
k1CG_1 <- k1CG_1[,c("k3_from", "to", "k1CG_from")]
k1CG_specific <- k1CG_2[k1CG_1, on = c("k1CG_from", "to")]
rm(k1CG_1, k1CG_2)

# kmer rates of change
k11_specific$k11_mu_rate <- k11_specific$k11_to_N / k11_specific$k11_from_N
k9_specific$k9_mu_rate <- k9_specific$k9_to_N / k9_specific$k9_from_N
k7_specific$k7_mu_rate <- k7_specific$k7_to_N / k7_specific$k7_from_N
k5_specific$k5_mu_rate <- k5_specific$k5_to_N / k5_specific$k5_from_N
k3_specific$k3_mu_rate <- k3_specific$k3_to_N / k3_specific$k3_from_N
k1_specific$k1_mu_rate <- k1_specific$k1_to_N / k1_specific$k1_from_N
k1CG_specific$k1CG_mu_rate <- k1CG_specific$k1CG_to_N / k1CG_specific$k1CG_from_N

k11_specific$k11_mu_rate[k11_specific$k11_from_N == 0] <- 0
k9_specific$k9_mu_rate[k9_specific$k9_from_N == 0] <- 0
k7_specific$k7_mu_rate[k7_specific$k7_from_N == 0] <- 0
k5_specific$k5_mu_rate[k5_specific$k5_from_N == 0] <- 0
k3_specific$k3_mu_rate[k3_specific$k3_from_N == 0] <- 0
k1_specific$k1_mu_rate[k1_specific$k1_from_N == 0] <- 0
k1CG_specific$k1CG_mu_rate[k1CG_specific$k1CG_from_N == 0] <- 0

# merge
train_specific <- k11_specific
train_specific$k9_from <- stri_sub(train_specific$k11_from, 2, -2)
train_specific$k7_from <- stri_sub(train_specific$k11_from, 3, -3)
train_specific$k5_from <- stri_sub(train_specific$k11_from, 4, -4)
train_specific$k3_from <- stri_sub(train_specific$k11_from, 5, -5)
train_specific$k1_from <- stri_sub(train_specific$k11_from, 6, -6)
train_specific <- k9_specific[train_specific, on = c("k9_from", "to")]
train_specific <- k7_specific[train_specific, on = c("k7_from", "to")]
train_specific <- k5_specific[train_specific, on = c("k5_from", "to")]
train_specific <- k3_specific[train_specific, on = c("k3_from", "to")]
train_specific <- k1_specific[train_specific, on = c("k1_from", "to")]
train_specific <- k1CG_specific[train_specific, on = c("k3_from", "to")]
rm(k1_specific, k1CG_specific, k3_specific, k5_specific, k7_specific, k9_specific, k11_specific, dt_80)

## merge test and train
names(test_specific)[names(test_specific) == 'k11_mu_rate'] <- 'k11_mu_rate_obs'
test_specific <- test_specific[,c("k11_from", "to", "k11_mu_rate_obs")]
train_specific <- train_specific[,c("k11_from",
                        "k1_from",
                        "to",
                        "k1CG_mu_rate",
                        "k1_mu_rate",
                        "k3_mu_rate",
                        "k5_mu_rate",
                        "k7_mu_rate",
                        "k9_mu_rate",
                        "k11_mu_rate")]
k11_dt <- test_specific[train_specific, on = c("k11_from", "to")]

## calculate absolute error between train and test for each k11
k11_dt$k1_AE <- abs(k11_dt$k1_mu_rate - k11_dt$k11_mu_rate_obs)
k11_dt$k1CG_AE <- abs(k11_dt$k1CG_mu_rate - k11_dt$k11_mu_rate_obs)
k11_dt$k3_AE <- abs(k11_dt$k3_mu_rate - k11_dt$k11_mu_rate_obs)
k11_dt$k5_AE <- abs(k11_dt$k5_mu_rate - k11_dt$k11_mu_rate_obs)
k11_dt$k7_AE <- abs(k11_dt$k7_mu_rate - k11_dt$k11_mu_rate_obs)
k11_dt$k9_AE <- abs(k11_dt$k9_mu_rate - k11_dt$k11_mu_rate_obs)
k11_dt$k11_AE <- abs(k11_dt$k11_mu_rate - k11_dt$k11_mu_rate_obs)

## calculate the total absolute error for each mutation type 
k11_dt_A <- k11_dt[k1_from == "A"]
out_A <- list()
for (j in 1:length(unique(k11_dt_A$to))){
  out_A[[j]] <- data.table(
    k1_from = k11_dt_A$k1_from[1],
    to = unique(k11_dt_A$to)[j],
    k1_AE = sum(k11_dt_A$k1_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]]),
    k1CG_AE = sum(k11_dt_A$k1CG_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]]),
    k3_AE = sum(k11_dt_A$k3_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]]),
    k5_AE = sum(k11_dt_A$k5_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]]),
    k7_AE = sum(k11_dt_A$k7_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]]),
    k9_AE = sum(k11_dt_A$k9_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]]),
    k11_AE = sum(k11_dt_A$k11_AE[k11_dt_A$to == unique(k11_dt_A$to)[j]])
  )
}

k11_dt_C <- k11_dt[k1_from == "C"]
out_C <- list()
for (j in 1:length(unique(k11_dt_C$to))){
  out_C[[j]] <- data.table(
    k1_from = k11_dt_C$k1_from[1],
    to = unique(k11_dt_C$to)[j],
    k1_AE = sum(k11_dt_C$k1_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]]),
    k1CG_AE = sum(k11_dt_C$k1CG_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]]),
    k3_AE = sum(k11_dt_C$k3_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]]),
    k5_AE = sum(k11_dt_C$k5_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]]),
    k7_AE = sum(k11_dt_C$k7_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]]),
    k9_AE = sum(k11_dt_C$k9_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]]),
    k11_AE = sum(k11_dt_C$k11_AE[k11_dt_C$to == unique(k11_dt_C$to)[j]])
  )
}
out_A_C <- do.call("rbind", c(out_A, out_C))

k11_dt_C_noCG <- k11_dt_C[which(stri_sub(k11_dt_C$k11_from, 6, -5) != "CG"),]
out_C_noCG <- list()
for (j in 1:length(unique(k11_dt_C_noCG$to))){
  out_C_noCG[[j]] <- data.table(
    k1_from = k11_dt_C_noCG$k1_from[1],
    to = unique(k11_dt_C_noCG$to)[j],
    k1_AE = sum(k11_dt_C_noCG$k1_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]]),
    k1CG_AE = sum(k11_dt_C_noCG$k1CG_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]]),
    k3_AE = sum(k11_dt_C_noCG$k3_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]]),
    k5_AE = sum(k11_dt_C_noCG$k5_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]]),
    k7_AE = sum(k11_dt_C_noCG$k7_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]]),
    k9_AE = sum(k11_dt_C_noCG$k9_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]]),
    k11_AE = sum(k11_dt_C_noCG$k11_AE[k11_dt_C_noCG$to == unique(k11_dt_C_noCG$to)[j]])
  )
}
out_C_noCG <- do.call("rbind", out_C_noCG)
out_C_noCG$CG_status <- "nonCG"

k11_dt_C_CG <- k11_dt_C[which(stri_sub(k11_dt_C$k11_from, 6, -5) == "CG"),]
out_C_CG <- list()
for (j in 1:length(unique(k11_dt_C_CG$to))){
  out_C_CG[[j]] <- data.table(
    k1_from = k11_dt_C_CG$k1_from[1],
    to = unique(k11_dt_C_CG$to)[j],
    k1_AE = sum(k11_dt_C_CG$k1_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]]),
    k1CG_AE = sum(k11_dt_C_CG$k1CG_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]]),
    k3_AE = sum(k11_dt_C_CG$k3_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]]),
    k5_AE = sum(k11_dt_C_CG$k5_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]]),
    k7_AE = sum(k11_dt_C_CG$k7_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]]),
    k9_AE = sum(k11_dt_C_CG$k9_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]]),
    k11_AE = sum(k11_dt_C_CG$k11_AE[k11_dt_C_CG$to == unique(k11_dt_C_CG$to)[j]])
  )
}
out_C_CG <- do.call("rbind", out_C_CG)
out_C_CG$CG_status <- "CG"
out_C_CG <- rbind(out_C_CG, out_C_noCG)

### EXPORT
fwrite(out_A_C, out.specific.AC.file, append = T)
fwrite(out_C_CG, out.specific.CCG.file, append = T)





