## 01. Import Library #########
library(multtest)
library(devtools)
library('MASS')
library(gplots)
library(compiler)
library(glue)
library(GAPIT3)
library(qqman)
library(rrBLUP)
args = commandArgs(trailingOnly=TRUE)

today <- "date"; print(today)

## 02. Set working directory ######
dir_0 <- "working directory"; dir_0
dir.create(path = dir_0)
setwd(dir_0)

## 03. copy files ####
# source code
source(file = "20210103_GAPIT_source_code_version_3_png.r")

# Phenotypic data
file.copy(from = "phenotypic data", to = getwd())

# Genotypic data
file.copy(from = glue("genotypic data"), to = getwd())

## 04. Sort phenotypic data ######
rda <- load("phenotypic data"); rda
pheno <- pheno
head(pheno); dim(pheno)
myY <- pheno[, c(1, 4)]; head(myY); colnames(myY); dim(myY)
myY <- pheno[, c(1, as.numeric(args[3]))]; head(myY); colnames(myY); dim(myY)
myY_rm_NA <- myY[!is.na(myY[, 2]), ]; dim(myY); dim(myY_rm_NA)
myY <- myY_rm_NA; head(myY)
colnames(myY)[1] <- "Taxa"
head(myY)

## 05. Sort Genotypinc data #####
list_geno <- list.files(path = getwd(), pattern = "genotype"); list_geno
selected_rda <- list_geno[1]; selected_rda
selected_rda <- list_geno[as.numeric(args[1])]; selected_rda
print(glue("selected rda: {selected_rda}"))
rda <- load(file = selected_rda); rda

geno_temp_2 <- geno_CC
geno_temp_2[1:10, 1:10]; dim(geno_temp_2)

myGM <- geno_temp_2[, c(1:3)]; head(myGM); dim(myGM)
sapply(myGM, class)
myGM[, 1] <- as.character(myGM[, 1])
myGM[, 2] <- as.character(myGM[, 2])
sapply(myGM, class)

colnames(myGM) <- c("SNP", "Chromosome", "Position")
sapply(myGM, class)

# myGD
geno_temp_desig_rownames <- geno_temp_2; geno_temp_desig_rownames[1:10, 1:10]
rownames(geno_temp_desig_rownames) <- geno_temp_desig_rownames$ID; geno_temp_desig_rownames[1:10, 1:10]
geno_temp_t <- t(geno_temp_desig_rownames[, -c(1:3)]); geno_temp_t[1:10, 1:10]
myGD_temp <- data.frame(rownames(geno_temp_t), geno_temp_t, stringsAsFactors = F); myGD_temp[1:10, 1:10]
rownames(myGD_temp) <- c(1:nrow(myGD_temp)); colnames(myGD_temp)[1] <- "taxa"; myGD_temp[1:5, 1:5]
colnames(myGD_temp) <- gsub(colnames(myGD_temp), pattern = "X", replacement = ""); myGD_temp[1:5, 1:5]
myGD_temp_2 <- matrix(NA, nrow(myGD_temp), ncol(myGD_temp)); myGD_temp_2[1:10, 1:10]
myGD_temp_2[myGD_temp == -1] <- 0; myGD_temp_2[1:10, 1:10]
myGD_temp_2[myGD_temp == 0] <- 1; myGD_temp_2[1:10, 1:10]
myGD_temp_2[myGD_temp == 1] <- 2; myGD_temp_2[1:10, 1:10]
myGD_temp_2 <- as.data.frame(myGD_temp_2); myGD_temp_2[1:10, 1:10]
myGD_temp_2[, 1] <- myGD_temp[, 1]; myGD_temp_2[1:10, 1:10]
colnames(myGD_temp_2) <- colnames(myGD_temp); myGD_temp_2[1:10, 1:10]
myGD <- myGD_temp_2

## 06. leave intersect accessions ####
intersectedName <- intersect(myGD[, 1], myY[, 1]); intersectedName; length(intersectedName)
myGD_final <- myGD[myGD[, 1] %in% intersectedName, ]; myGD_final[1:10, 1:10]; dim(myGD_final)
myY_final <- myY[myY[, 1] %in% intersectedName, ]; head(myY_final); dim(myY_final)

## 07. Set save name ####
phenoName <- "phenoName"
genoName <- "genoName"

## 08. Kinship matrix ####
print("# Kinship matrix 만들기 ---------------")
g <- geno_temp_2[, -c(1:3)]; g[1:10, 1:10]
n.core <- 30
K <- A.mat(t(g), shrink = T, n.core = n.core)

print("PCA -----------------")
pca.res <- prcomp(K, scale = TRUE)
print("PC_number ------------------------")
PCA.X <- prcomp(K)
eigenvalues <- PCA.X$sdev^2
evp=eigenvalues/sum(eigenvalues)
nout=min(10,length(evp))
xout=1:nout

for(i in c(1:length(eigenvalues))){
  if(i == 1){
    temp <- 0
  }
  temp <- temp + eigenvalues[i]
  if((temp/sum(eigenvalues)*100) > 95){
    break
  }
}
PC_number <- i; print(paste0("PC_number: ", i))

## 9. models ####
model = c("FarmCPU", "MLMM", "GLM", "MLM", "CMLM", "ECMLM", "BLINK") 
model_sel <- model[as.numeric(args[2])]; model_sel


## 10. Run GAPIT 3 #### 
save_dir <- "save directory"
dir.create(save_dir)
setwd(save_dir)

set.seed(961022)
myGAPIT <- GAPIT(Y = myY_final,
                 GD = myGD_final,
                 GM = myGM,
                 PCA.total = PC_number, ### 얼마나 많은 PC 개수를 고려할 것인지를 지정
                 # cutOff=0.05,
                 model = model_sel
)

save(myGAPIT, file = glue("{today}_10_myGAPIT_{model_sel}_result_geno-{genoName}_pheno-{phenoName}.rda"))
save.image(glue("{today}_17_GAPIT_environment_{model_sel}_geno-{genoName}_pheno-{phenoName}.RData"))

## 11. Sort GAPIT result ####
head(myGAPIT$GWAS)
amfw <- matrix(NA, nrow = nrow(myGAPIT$GWAS), ncol = 4)
amfw <- as.data.frame(amfw)
colnames(amfw) <- c("SNP", "CHR", "BP", "P")
amfw$SNP <- myGAPIT$GWAS$SNP; head(amfw)
amfw$CHR <- as.numeric(myGAPIT$GWAS$Chromosome); head(amfw)
amfw$BP <- as.numeric(myGAPIT$GWAS$Position); head(amfw)
amfw$P <- myGAPIT$GWAS$P.value; head(amfw)

## 12. Set figure size ####
width <- 10
height <- 6.6666666

## 13. Set Bonferroni correction ####
n.markers <- nrow(amfw)
bonf <- -log10(0.05/ n.markers)

if(sum(is.infinite(-log10(amfw$P)))){
  ylim_max_temp <- (-log10(amfw$P)); ylim_max_temp
  ylim_max_temp_2 <- ylim_max_temp[which(ylim_max_temp != Inf)]; ylim_max_temp_2
  
  ylim_max <- floor(max(ylim_max_temp_2)+1); ylim_max
}else{
  ylim_max <- floor(max(-log10(amfw$P))+1); ylim_max
}

png(filename = glue("{today}_11_GAPIT_GWAS_BonfferoniCorrection_{model_sel}_geno-{genoName}_pheno-{phenoName}.png"), width = width, height = height, units = "in", res = 300)
try(manhattan(amfw, 
              main = glue("{model_sel}_{colnames(myY_final)[2]}"), 
              suggestiveline = FALSE, 
              genomewideline = FALSE, 
              ylim = c(0, ylim_max)) + 
      abline(h = bonf, col = 'blue', lty = 2), silent = T)
dev.off()

## 14. False Discovery Rate (FDR, BH method) ####
bh.res <- p.adjust(amfw$P, method = "BH")
sig <- amfw[bh.res < 0.05,]
bottom.bh <- min(-log10(sig$P))

png(filename = glue("{today}_12_GAPIT_GWAS_FDR-BH_{model_sel}_geno-{genoName}_pheno-{phenoName}.png"), width = width, height = height, units = "in", res = 300)
try(manhattan(amfw, 
              main = glue("{model_sel}_{colnames(myY_final)[2]}"), 
              suggestiveline = FALSE, 
              genomewideline = FALSE, 
              ylim = c(0, ylim_max)) + 
      abline(h= bottom.bh, col = "forestgreen", lty = 2), silent = T)
dev.off()

## 15. Export significant SNPs ####
head(amfw)
`-log10(P)` <- -log10(x = amfw$P); head(`-log10(P)`)
`amfw_-log10(P)` <- cbind(amfw, `-log10(P)`); head(`amfw_-log10(P)`)

amfw_sig <- `amfw_-log10(P)`[`amfw_-log10(P)`$`-log10(P)` > bonf, ]; amfw_sig; dim(amfw_sig)
name_amfw_sig <- glue("{today}_14_GAPIT_sigSNPs_{model_sel}_geno-{genoName}_pheno-{phenoName}.tsv");name_amfw_sig
write.table(x = amfw_sig, file = name_amfw_sig, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

head(myY_final); nrow(myY_final)
head(amfw_sig); dim(amfw_sig)
colnames(myGD_final) <- gsub(pattern = "^X", replacement = "", x = colnames(myGD_final))

## 16. ANOVA significant SNPs #####
amfw_sig_anova <- matrix(data = NA, nrow = 1, ncol = (ncol(amfw_sig)+1)); amfw_sig_anova
if (nrow(amfw_sig) != 0){
  if (sum(myGD_final[, 1] != myY_final[, 1]) == 0){
    colnames(myGD_final) <- gsub(pattern = ".HRSCAF", replacement = ";HRSCAF", x = colnames(myGD_final))
    pheno_geno_df <- cbind(myY_final, myGD_final[, amfw_sig[, 1]]);head(pheno_geno_df); dim(pheno_geno_df)
    colnames(pheno_geno_df) <- c(colnames(myY_final), amfw_sig[, 1]); head(pheno_geno_df)
  }else{
    temp_myY <- myY_final
    rownames(temp_myY) <- temp_myY[, 1]
    ordered_myY <- data.frame(myGD_final[, 1], temp_myY[myGD_final[, 1], ][,2]); head(ordered_myY)
    colnames(ordered_myY) <- colnames(myY_final)
    myY_final <- ordered_myY
    
    colnames(myGD_final) <- gsub(pattern = ".HRSCAF", replacement = ";HRSCAF", x = colnames(myGD_final))
    pheno_geno_df <- cbind(myY_final, myGD_final[, amfw_sig[, 1]]);head(pheno_geno_df); dim(pheno_geno_df)
    colnames(pheno_geno_df) <- c(colnames(myY_final), amfw_sig[, 1]); head(pheno_geno_df)
  }
  
  marker <- amfw_sig[, 1][1]; marker
  for(marker in amfw_sig[, 1]){
    print(marker)
    aov_sum <- summary(aov(pheno_geno_df[, 2] ~ pheno_geno_df[, marker], data = pheno_geno_df)); aov_sum
    
    if(aov_sum[[1]]$`Pr(>F)`[1] <= 0.05){
      print(glue("Significant: {marker}"))
      amfw_sig_anova <- rbind(amfw_sig_anova, unlist(c(amfw_sig[which(amfw_sig[, 1] == marker), ], aov_sum[[1]]$`Pr(>F)`[1]))); amfw_sig_anova
    }else{
      print(glue("Non-significant: {marker}"))
    }
  }
  
}

name_amfw_sig_anova <- glue("{today}_15_GAPIT_sigSNPs_anova_{model_sel}_geno-{genoName}_pheno-{phenoName}.tsv");name_amfw_sig_anova
write.table(x = amfw_sig_anova, file = name_amfw_sig_anova, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

png(filename = glue("{today}_16_GAPIT_GWAS_Bonf_sig_anova_{model_sel}_geno-{genoName}_pheno-{phenoName}.png"), width = width, height = height, units = "in", res = 300)
try(manhattan(amfw, 
              main = glue("{model_sel}_{colnames(myY_final)[2]}"), 
              suggestiveline = FALSE, 
              genomewideline = FALSE, 
              ylim = c(0, ylim_max),
              highlight = amfw_sig_anova[, 1]) + 
      abline(h = bonf, col = 'blue', lty = 2), silent = T)
dev.off()

## 17. environment 저장하기 #######
save.image(glue("{today}_17_GAPIT_environment_{model_sel}_geno-{genoName}_pheno-{phenoName}.RData"))
