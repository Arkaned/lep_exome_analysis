rm(list = ls())

setwd("~/Documents/lab/butterfly_evol/")

library(phytools)
library(caper)
library(geodata)
library(stringi)
library(dplyr)
library(tidyr)
library(CoordinateCleaner)
library(rgbif)
library(lmerTest)

## loading phylogeny
tr <- read.tree("data/bs.tre")

## frequency data
freq_all <- read.csv("data/freq_all_stand.csv", row.names = 1)
freq_ort <- read.csv("data/freq_ort_stand.csv", row.names = 1)
freq_non <- read.csv("data/freq_non_stand.csv", row.names = 1)

pca_all <- prcomp(t(freq_all))
sum_all <- summary(pca_all)
pcdata_all <- as.data.frame(pca_all$x[, 1:3])
pcloa_all <- pca_all$rotation[, 1:3]
write.csv(pcloa_all, "tables/pca_loadings_all.csv")

pca_ort <- prcomp(t(freq_ort))
sum_ort <- summary(pca_ort)
pcdata_ort <- as.data.frame(pca_ort$x[, 1:3])
pcloa_ort <- pca_ort$rotation[, 1:3]
write.csv(pcloa_ort, "tables/pca_loadings_ort.csv")

pca_non <- prcomp(t(freq_non))
sum_non <- summary(pca_non)
pcdata_non <- as.data.frame(pca_non$x[, 1:3])
pcloa_non <- pca_non$rotation[, 1:3]
write.csv(pcloa_non, "tables/pca_loadings_non.csv")

## ecological data
diet <- read.csv("data/Insect_DBR_AAfreq_clean.csv", row.names = 1)
diet$species <- stri_replace_all_fixed(diet$species, " ", "_")

mod <- lmer(sp_freq ~ numb_codons + (1|species), data = diet) 
summary(mod)
anova(mod)

## occurrence data
dat <- read.csv("data/Insect_clean_DBR_Exome.csv", row.names = 1)
dat$species <- stri_replace_all_fixed(dat$species, " ", "_")

dat$species[dat$species == "Aglais_io"] <- "Nymphalis_io"
dat$species[dat$species == "Mycalesis_anynana"] <- "Bicyclus_anynana"

# query current environmental data
clim_current <- worldclim_global("bio", res = 2.5, path = tempdir())
names(clim_current) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
                         "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", 
                         "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
clim_current <- raster::stack(clim_current)

## extracting data for each occurrence
extract_env <- function(spp, data, clim, diet) {
  
  dat_sub <- data[data$species == spp, ]
  
  values_env <- terra::extract(clim, dat_sub[, c("decimalLongitude", "decimalLatitude")])
  values_env <- as.data.frame(values_env)
  
  temp_niche_breadth <- (max(values_env$bio5, na.rm = T)) - 
    (min(values_env$bio6, na.rm = T))
  names(temp_niche_breadth) <- "temp_niche_breadth"
  prec_niche_breadth <- max(values_env$bio16, na.rm = T) - 
    min(values_env$bio17, na.rm = T)
  names(prec_niche_breadth) <- "prec_niche_breadth"
  max_temp <- max(values_env$bio5, na.rm = T)
  min_temp <- min(values_env$bio6, na.rm = T)
  prec_max <- max(values_env$bio16, na.rm = T)
  prec_min <- min(values_env$bio17, na.rm = T)
  
  diet <- as.numeric(names(table(diet[diet$species == stri_replace_all_fixed(spp, " ", "_"), "n"])))
  
  res <- data.frame(Species = spp,
                    TNB = temp_niche_breadth,
                    PNB = prec_niche_breadth,
                    Temp_max = max_temp,
                    Temp_min = min_temp,
                    Prec_max = prec_max,
                    Prec_min = prec_min,
                    DNB = diet)
  
  return(res)
}

spp <- unique(dat$species)

dat_clim_all <- dat_clim_ort <- dat_clim_non <- as.data.frame(matrix(nrow = length(spp), ncol = 11))
dat_clim_all[, 1] <- dat_clim_ort[, 1] <- dat_clim_non[, 1] <- spp
colnames(dat_clim_all) <- colnames(dat_clim_ort) <- colnames(dat_clim_non) <- 
  c("Spp", "TNB", "PNB", "Temp_max", "Temp_min", "Prec_max", "Prec_min", "DNB", 
    "PC1", "PC2", "PC3")
for (i in 1:length(spp)) {
  dat_clim_all[i, 2:8] <- extract_env(spp[i], dat, clim_current, diet)[, -1]
  dat_clim_all[i, 9:11] <- pcdata_all[rownames(pcdata_all) == spp[i], ]
  
  dat_clim_ort[i, 2:8] <- extract_env(spp[i], dat, clim_current, diet)[, -1]
  dat_clim_ort[i, 9:11] <- pcdata_ort[rownames(pcdata_ort) == spp[i], ]
  
  dat_clim_non[i, 2:8] <- extract_env(spp[i], dat, clim_current, diet)[, -1]
  dat_clim_non[i, 9:11] <- pcdata_non[rownames(pcdata_non) == spp[i], ]
}

write.csv(dat_clim_ort, "data/dat_clim_ort.csv")

stats_pgls <- function(x) {
  
  df <- as.data.frame(matrix(ncol = 10, nrow = length(x)))
  colnames(df) <- c("intercept", 
                    "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                    "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                    "R_squared")
  for (i in 1:length(x)) {
    df[i, 1] <- summary(x[[i]])$coefficients[1, 1]
    df[i, 2] <- summary(x[[i]])$coefficients[2, 1]
    df[i, 3] <- summary(x[[i]])$coefficients[2, 2]
    df[i, 4] <- summary(x[[i]])$coefficients[2, 3]
    df[i, 5] <- summary(x[[i]])$coefficients[2, 4]
    df[i, 6] <- summary(x[[i]])$coefficients[3, 1]
    df[i, 7] <- summary(x[[i]])$coefficients[3, 2]
    df[i, 8] <- summary(x[[i]])$coefficients[3, 3]
    df[i, 9] <- summary(x[[i]])$coefficients[3, 4]
    df[i, 10] <- summary(x[[i]])$r.squared
  }
  
  return(df)
  
}

pgls_pca_all <- as.data.frame(matrix(ncol = 10, nrow = 3))
colnames(pgls_pca_all) <- c("intercept", 
                            "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                            "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                            "R_squared")
rownames(pgls_pca_all) <- c("PC1", "PC2", "PC3")
for (i in 1:3) {
  dat_i <- dat_clim_all[, c(1, 4, 8)]
  dat_i <- cbind(dat_i, 
                 PC = dat_clim_all[, colnames(dat_clim_all) == rownames(pgls_pca_all)[i]])
  pgls_i <- list()
  for (j in 1:length(tr)) {
    comp_dat_j <- comparative.data(data = dat_i, phy = tr[[j]],
                                   names.col = "Spp", vcv.dim = 2)
    pgls_i[[j]] <- pgls(PC ~ log(DNB) + log(Temp_max), data = comp_dat_j)
  }
  
  pgls_stats_i <- stats_pgls(pgls_i)
  
  pgls_pca_all[i, ] <- rbind(
    paste0(apply(pgls_stats_i, 2, function(x) round(median(x), 3)), " (", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"))
}
write.csv(pgls_pca_all, "tables/pgls_pca_all.csv")

pgls_pca_ort <- as.data.frame(matrix(ncol = 10, nrow = 3))
colnames(pgls_pca_ort) <- c("intercept", 
                        "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                        "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                        "R_squared")
rownames(pgls_pca_ort) <- c("PC1", "PC2", "PC3")
for (i in 1:3) {
  dat_i <- dat_clim_ort[, c(1, 4, 8)]
  dat_i <- cbind(dat_i, 
                 PC = dat_clim_ort[, colnames(dat_clim_ort) == rownames(pgls_pca_ort)[i]])
  pgls_i <- list()
  for (j in 1:length(tr)) {
    comp_dat_j <- comparative.data(data = dat_i, phy = tr[[j]],
                                   names.col = "Spp", vcv.dim = 2)
    pgls_i[[j]] <- pgls(PC ~ log(DNB) + log(Temp_max), data = comp_dat_j)
    summary(pgls_i[[905]])$coefficients
  }
  
  pgls_stats_i <- stats_pgls(pgls_i)
  
  pgls_pca_ort[i, ] <- rbind(
    paste0(apply(pgls_stats_i, 2, function(x) round(median(x), 3)), " (", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"))
}
write.csv(pgls_pca_ort, "tables/pgls_pca_ort.csv")

pgls_pca_non <- as.data.frame(matrix(ncol = 10, nrow = 3))
colnames(pgls_pca_non) <- c("intercept", 
                            "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                            "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                            "R_squared")
rownames(pgls_pca_non) <- c("PC1", "PC2", "PC3")
for (i in 1:3) {
  dat_i <- dat_clim_non[, c(1, 4, 8)]
  dat_i <- cbind(dat_i, 
                 PC = dat_clim_non[, colnames(dat_clim_non) == rownames(pgls_pca_non)[i]])
  pgls_i <- list()
  for (j in 1:length(tr)) {
    comp_dat_j <- comparative.data(data = dat_i, phy = tr[[j]],
                                   names.col = "Spp", vcv.dim = 2)
    pgls_i[[j]] <- pgls(PC ~ log(DNB) + log(Temp_max), data = comp_dat_j)
  }
  
  pgls_stats_i <- stats_pgls(pgls_i)
  
  pgls_pca_non[i, ] <- rbind(
    paste0(apply(pgls_stats_i, 2, function(x) round(median(x), 3)), " (", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"))
}
write.csv(pgls_pca_non, "tables/pgls_pca_non.csv")

#plot(fitted(pgls_i[[1]]), residuals(pgls_i[[1]]))
#abline(h = 0, col = "red")

#qqnorm(residuals(pgls_i[[1]]), main = "Q-Q Plot")
#qqline(residuals(pgls_i[[1]]), col="red")

#hist(residuals(pgls_i[[1]]))

#plot(log(comp_dat_j$data$Temp_max), residuals(pgls_i[[1]]), 
#     xlab = "Temperatura", 
#     ylab = "Resíduos", 
#     main = "Resíduos vs. Temperatura")
#abline(h = 0, col = "red")
#plot(log(comp_dat_j$data$DNB), residuals(pgls_i[[1]]), 
#     xlab = "Temperatura", 
#     ylab = "Resíduos", 
#     main = "Resíduos vs. Temperatura")
#abline(h = 0, col = "red")

spp <- unique(dat$species)

dat_clim_all_sep <- dat_clim_ort_sep <- dat_clim_non_sep <- list()
for (i in 1:length(spp)) {
  dat_clim_all_sep[[i]] <- as.data.frame(matrix(nrow = 20, ncol = 10))
  dat_clim_all_sep[[i]][1:20, 1:8] <- extract_env(spp[i], dat, clim_current,
                                                  diet)
  dat_clim_all_sep[[i]][, 9] <- freq_all[, spp[i]]
  dat_clim_all_sep[[i]][, 10] <- rownames(freq_all)
  
  dat_clim_ort_sep[[i]] <- as.data.frame(matrix(nrow = 20, ncol = 10))
  dat_clim_ort_sep[[i]][1:20, 1:8] <- extract_env(spp[i], dat, clim_current,
                                                  diet)
  dat_clim_ort_sep[[i]][, 9] <- freq_ort[, spp[i]]
  dat_clim_ort_sep[[i]][, 10] <- rownames(freq_ort)
  
  dat_clim_non_sep[[i]] <- as.data.frame(matrix(nrow = 20, ncol = 10))
  dat_clim_non_sep[[i]][1:20, 1:8] <- extract_env(spp[i], dat, clim_current,
                                                  diet)
  dat_clim_non_sep[[i]][, 9] <- freq_non[, spp[i]]
  dat_clim_non_sep[[i]][, 10] <- rownames(freq_non)
}

dat_clim_all_sep <- do.call(rbind, dat_clim_all_sep)
dat_clim_ort_sep <- do.call(rbind, dat_clim_ort_sep)
dat_clim_non_sep <- do.call(rbind, dat_clim_non_sep)

colnames(dat_clim_all_sep)[1:8] <- colnames(extract_env(spp[1], dat, clim_current,
                                                        diet))
colnames(dat_clim_all_sep)[9:10] <- c("prop", "aa")
colnames(dat_clim_ort_sep)[1:8] <- colnames(extract_env(spp[1], dat, clim_current,
                                                    diet))
colnames(dat_clim_ort_sep)[9:10] <- c("prop", "aa")
colnames(dat_clim_non_sep)[1:8] <- colnames(extract_env(spp[1], dat, clim_current,
                                                    diet))
colnames(dat_clim_non_sep)[9:10] <- c("prop", "aa")

stats_pgls_sep <- function(x, n) {
  
  df <- as.data.frame(matrix(ncol = 10, nrow = length(x)))
  colnames(df) <- c("intercept", 
                    "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                    "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                    "R_squared")
  for (i in 1:length(x)) {
    df[i, 1] <- summary(x[[i]])$coefficients[1, 1]
    df[i, 2] <- summary(x[[i]])$coefficients[2, 1]
    df[i, 3] <- summary(x[[i]])$coefficients[2, 2]
    df[i, 4] <- summary(x[[i]])$coefficients[2, 3]
    df[i, 5] <- p.adjust(summary(x[[i]])$coefficients[2, 4], method = "BH", n = n)
    df[i, 6] <- summary(x[[i]])$coefficients[3, 1]
    df[i, 7] <- summary(x[[i]])$coefficients[3, 2]
    df[i, 8] <- summary(x[[i]])$coefficients[3, 3]
    df[i, 9] <- p.adjust(summary(x[[i]])$coefficients[3, 4], method = "BH", n = n)
    df[i, 10] <- summary(x[[i]])$r.squared
  }
  
  return(df)
  
}

aa <- unique(dat_clim_ort_sep$aa)

pgls_all_sep <- as.data.frame(matrix(ncol = 10, nrow = length(aa)))
colnames(pgls_all_sep) <- c("intercept", 
                            "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                            "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                            "R_squared")
rownames(pgls_all_sep) <- aa
for (i in 1:length(aa)) {
  dat_i <- dat_clim_all_sep[dat_clim_all_sep$aa == aa[i], ]
  pgls_i <- list()
  for (j in 1:length(tr)) {
    comp_dat_j <- comparative.data(data = dat_i, phy = tr[[j]],
                                   names.col = "Species", vcv.dim = 2)
    pgls_i[[j]] <- pgls(scale(prop) ~ log(DNB) + log(Temp_max), data = comp_dat_j)
  }
  
  pgls_stats_i <- stats_pgls_sep(pgls_i, n = length(aa)*3)
  
  pgls_all_sep[i, ] <- rbind(
    paste0(apply(pgls_stats_i, 2, function(x) round(median(x), 3)), " (", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"))
}
write.csv(pgls_all_sep, "tables/pgls_all_sep.csv")

pgls_ort_sep <- as.data.frame(matrix(ncol = 10, nrow = length(aa)))
colnames(pgls_ort_sep) <- c("intercept", 
                        "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                        "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                        "R_squared")
rownames(pgls_ort_sep) <- aa
for (i in 1:length(aa)) {
  dat_i <- dat_clim_ort_sep[dat_clim_ort_sep$aa == aa[i], ]
  pgls_i <- list()
  for (j in 1:length(tr)) {
    comp_dat_j <- comparative.data(data = dat_i, phy = tr[[j]],
                                   names.col = "Species", vcv.dim = 2)
    pgls_i[[j]] <- pgls(scale(prop) ~ log(DNB) + log(Temp_max), data = comp_dat_j)
  }
  
  pgls_stats_i <- stats_pgls_sep(pgls_i, n = length(aa)*3)
  
  pgls_ort_sep[i, ] <- rbind(
    paste0(apply(pgls_stats_i, 2, function(x) round(median(x), 3)), " (", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"))
}
write.csv(pgls_ort_sep, "tables/pgls_ort_sep.csv")

pgls_non_sep <- as.data.frame(matrix(ncol = 10, nrow = length(aa)))
colnames(pgls_non_sep) <- c("intercept", 
                            "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                            "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                            "R_squared")
rownames(pgls_non_sep) <- aa
for (i in 1:length(aa)) {
  dat_i <- dat_clim_non_sep[dat_clim_non_sep$aa == aa[i], ]
  pgls_i <- list()
  for (j in 1:length(tr)) {
    comp_dat_j <- comparative.data(data = dat_i, phy = tr[[j]],
                                   names.col = "Species", vcv.dim = 2)
    pgls_i[[j]] <- pgls(scale(prop) ~ log(DNB) + log(Temp_max), data = comp_dat_j)
  }
  
  pgls_stats_i <- stats_pgls_sep(pgls_i, n = length(aa)*3)
  
  pgls_non_sep[i, ] <- rbind(
    paste0(apply(pgls_stats_i, 2, function(x) round(median(x), 3)), " (", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
           apply(pgls_stats_i, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"))
}
write.csv(pgls_non_sep, "tables/pgls_non_sep.csv")
