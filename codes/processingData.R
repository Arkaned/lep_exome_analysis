rm(list = ls())

setwd("~/Documents/lab/butterfly_evol/")

library(phytools)
library(stringi)

## loading frequency aa as percentage for orthologs and non-orthologs
freq_all <- read.csv("data/aa_perce.csv")
freq_ort <- read.csv("data/aa_perce_ortho.csv")
freq_non <- read.csv("data/aa_perce_non-ortho.csv")

## tirar stop codon and unknown aminoacids
freq_all <- freq_all[!(freq_all$AA.FREQ %in% c("X", "*")), ]
freq_ort <- freq_ort[!(freq_ort$AA %in% c("X", "*")), ]

freq <- read.csv("data/Insect_DBR_AAfreq_clean.csv", row.names = 1)
freq$species <- stri_replace_all_fixed(freq$species, " ", "_")

## changing the names
for (i in 1:ncol(freq_all)) {
  spp_i <- colnames(freq_all)[i]
  x <- strsplit(spp_i, "[.]")[[1]]
  
  colnames(freq_all)[i] <- x[1]
}
freq_all <- freq_all[, -37]

for (i in 1:ncol(freq_ort)) {
  spp_i <- colnames(freq_ort)[i]
  x <- strsplit(spp_i, "[.]")[[1]]
  
  colnames(freq_ort)[i] <- x[1]
}

for (i in 1:ncol(freq_non)) {
  spp_i <- colnames(freq_non)[i]
  x <- strsplit(spp_i, "[.]")[[1]]
  
  colnames(freq_non)[i] <- x[1]
}

## standardizing the data
freq_all_stand <- freq_all[, -1]
rownames(freq_all_stand) <- freq_all[, 1]
for (i in 21:nrow(freq_all_stand)) {
  for (j in 1:ncol(freq_all_stand)) {
    num_i <- freq_all_stand[i, j]
    spp_i <- colnames(freq_all_stand)[j]
    aa_i <- rownames(freq_all_stand)[i]
    
    div_i <- freq$numb_codons[freq$species == spp_i & freq$aa_letter == aa_i]
    
    if (length(div_i) == 1) {
      freq_all_stand[i, j] <- num_i/div_i
    } 
    
    if (length(div_i) == 0) {
      freq_all_stand[i, j] <- NA
    }
    
  }
}

freq_ort_stand <- freq_ort[, -1]
rownames(freq_ort_stand) <- freq_ort[, 1]
for (i in 21:nrow(freq_ort_stand)) {
  for (j in 1:ncol(freq_ort_stand)) {
    num_i <- freq_ort_stand[i, j]
    spp_i <- colnames(freq_ort_stand)[j]
    aa_i <- rownames(freq_ort_stand)[i]
    
    div_i <- freq$numb_codons[freq$species == spp_i & freq$aa_letter == aa_i]
    
    if (length(div_i) == 1) {
      freq_ort_stand[i, j] <- num_i/div_i
    } 
    
    if (length(div_i) == 0) {
      freq_ort_stand[i, j] <- NA
    }
  
    }
}

freq_non_stand <- freq_non[, -1]
rownames(freq_non_stand) <- freq_non[, 1]
for (i in 21:nrow(freq_non_stand)) {
  for (j in 1:ncol(freq_non_stand)) {
    num_i <- freq_non_stand[i, j]
    spp_i <- colnames(freq_non_stand)[j]
    aa_i <- rownames(freq_non_stand)[i]
    
    div_i <- freq$numb_codons[freq$species == spp_i & freq$aa_letter == aa_i]
    
    if (length(div_i) == 1) {
      freq_non_stand[i, j] <- num_i/div_i
    } 
    
    if (length(div_i) == 0) {
      freq_non_stand[i, j] <- NA
    }
    
  }
}

write.csv(freq_all_stand, "data/freq_all_stand.csv")
write.csv(freq_ort_stand, "data/freq_ort_stand.csv")
write.csv(freq_non_stand, "data/freq_non_stand.csv")
