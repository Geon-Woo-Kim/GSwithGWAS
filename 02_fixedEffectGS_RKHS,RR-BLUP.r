# library ----
library(vcfR)
library(regress)
library(synbreed)
library(tidyverse)
library(glue)
library(rrBLUP)
library(glmnet)
library(VIGoR)
library(randomForest)
args = commandArgs(trailingOnly=TRUE)

# function -----
# Design matrix
make.design.mat <- function(row.id, col.id){
  Z <- matrix(0, nrow=length(row.id), ncol=length(col.id))
  rownames(Z) <- row.id
  colnames(Z) <- col.id	
  Z[cbind(1:nrow(Z),match(row.id,col.id))] <- 1
  return(Z)
} # make.design.mat()

# directory ------
dir0 <- "/home/rstudio/Public/20_KGW/12_GSwithGWAS/20220107_01_modelingAndPrediction_CCannuumComplex_HN_capsaicinoidBLUP_LD0.30AndGWAS/"
setwd(dir0); getwd()

# training set genotypic data
rda <- load("training and test set genomic data"); rda
g <- t(geno_CC[,-(1:3)])
g[1:10,1:5]
snp_names <- geno_CC[, 1]; head(snp_names); length(snp_names); ncol(g)
class(snp_names)
snp_names <- as.character(snp_names); class(snp_names); head(snp_names)

# test population genotype -------
# geno_HN <- geno_HL
geno_HN[1:10, 1:10]
g_pre <- t(geno_HN[, -(1:3)]); g_pre[1:10, 1:10]
rownames(g_pre)
nameorder <- order(rownames(g_pre)); nameorder
g_pre_order <- g_pre[nameorder, ]
rownames(g_pre_order)
dim(g_pre)

# training set phenotypic data ------
file.copy(from = "/home/rstudio/Public/20_KGW/08_Phenotype_data/coreCollection/20201012_phenotype_coreCollection_capsacinoids_mean_blup_16to18_rmOutlier/20201012_phenotype_CCannuumComplex284_capsaicinoids_mean_BLUP_rmOutlier.rda", to = getwd())
rda <- load("20201012_phenotype_CCannuumComplex284_capsaicinoids_mean_BLUP_rmOutlier.rda"); rda
head(pheno); colnames(pheno)
which(is.na(pheno$totalCapBLUP_ug_g))
pheno <- pheno[which(!is.na(pheno$totalCapBLUP_ug_g)), ]
dim(pheno)
select_col <- c(4)
i <- 4
for(i in select_col){
  if(i == select_col[1]){
    list_y <- list()
  }
  list_y[[which(i == select_col)]] <- pheno[, i]
  names(list_y[[which(i == select_col)]]) <- pheno[,1]
  list_y[[which(i == select_col)]] <- list_y[[which(i == select_col)]][!is.na(list_y[[which(i == select_col)]])]
  if(i == select_col[length(select_col)]){
    names(list_y) <- colnames(pheno[select_col])
  }
}
names(list_y); length(list_y[[1]])

# testing set phenotypic data for calculating accuracy of model -------
file.copy(from = "/home/rstudio/Public/20_KGW/08_Phenotype_data/HanaEliteLine/20200930_phenotype_hanaEliteLine_capsaicinoids_average_final.rda", to = getwd())
rda <- load("20200930_phenotype_hanaEliteLine_capsaicinoids_average_final.rda"); rda
head(phenotypeAverage)
phenotypeAverage$Name
pheno_temp <- phenotypeAverage
colnames(pheno_temp)[1] <- "Name"
head(pheno_temp)
pheno_temp$Name[1:9] <- paste("HNAL-0", substr(pheno_temp$Name[1:9], 4,4), sep = ""); pheno_temp$Name[1:9]
pheno_temp$Name[10:length(pheno_temp$Name)] <- paste("HNAL-", substr(pheno_temp$Name[10:length(pheno_temp$Name)], 4,5), sep = "")
pheno_temp$Name
order(pheno_temp$Name)
pheno_temp <- pheno_temp[order(pheno_temp$Name), ]; pheno_temp
pheno <- data.frame(pheno_temp)
for(i in 2:ncol(pheno)){
  pheno[, i] <- as.numeric(pheno[, i])
}
head(pheno)
head(pheno)
ncol(pheno)
pheno_trait <- pheno[, c(1, 4)]; dim(pheno_trait)
pheno_trait_temp <- pheno_trait[which(!is.na(pheno_trait[, 2])), ]; dim(pheno_trait_temp)
pheno_trait <- pheno_trait_temp
head(pheno_trait); dim(pheno_trait)

names(list_y)
i <- 1
y <- list_y[[i]]; head(y)
intersect_sample <- intersect(x = names(y), y = rownames(g)); length(intersect_sample)
nrow(g)
g <- g[intersect_sample, ]; nrow(g)
length(y)
sum(rownames(g) != names(y))

# naming -----
today <- "date"
trainingPop <- glue("name of training population"); trainingPop
testPop <- glue("name of testing population"); testPop
traitName <- "trait"; traitName

# make design matrix relating observations to lines in the training set
Z <- make.design.mat(row.id=names(y), col.id=rownames(g))
length(y)
dim(Z)
dim(g)

# format of genotypic and phenotypic data for GS modeling ------
Zg <- Z %*% g
hist(y, breaks = 20)
hist(log10(y), breaks = 20)
dev.off()
temp.y <- y

NAs <- rep(NA)
Models <- list(gblupRR = NAs,
               RKHS = NAs
               )
NAs <- rep(NA, nrow(Zg))
GEBVs <- data.frame(
  IDs=rownames(Zg),
  Yobs=NAs,
  gblupRR=NAs,
  RKHS=NAs
)
GEBVs$Yobs <- temp.y
head(GEBVs)

# fixed effect ------
csv <- read.table(file = "fixed effect list", header = T, sep = "\t", stringsAsFactors = F); head(csv)

snp_num <- which(snp_names %in% csv[, 1]); snp_num
sum(snp_names[snp_num] %in% csv[, 1]) == sum(csv[, 1] %in% snp_names[snp_num]) # significant 한 마커가 모두 선택되었는지 확인하기.


save_dir <- "save directory"
dir.create(path = save_dir)

for(j in 1:length(snp_num)){ # 1개부터 significant marker 개수 만큼 몇 개의 marker를 fixed effect로 지정할 것인가?
# for(j in 1){ # fixed effect로 하나의 marker만 지정하려 함.
  combn_SNPs <- combn(x = snp_num, m = j); combn_SNPs
  for(k in 1:ncol(combn_SNPs)){
    print(glue("선택한 SNP 개수 (j): {j}/{length(snp_num)}, 순서 (k): {k}/{ncol(combn_SNPs)}"))
    fixed_SNP <- combn_SNPs[, k]; print(glue("SNP num: {fixed_SNP}"))
    selected_SNP <- snp_names[fixed_SNP]; print(glue("SNP name: {selected_SNP}"))
    fixed_SNP_vec <- paste0(fixed_SNP, collapse = ", "); fixed_SNP_vec; print(glue("SNP num: {fixed_SNP_vec}")) 
    selected_SNP_vec <- paste0(selected_SNP, collapse = ", "); selected_SNP_vec; print(glue("SNP name: {selected_SNP_vec}"))
    
    # model construction -------
    # GBLUP-RR (gblupRR)
    print(glue("gblupRR training start"))
    res <- kinship.BLUP(y=temp.y, G.train=Zg, G.pred = g_pre, K.method="RR", X = Zg[, fixed_SNP])
    GEBVs$gblupRR <- res$g.train + as.numeric(res$beta)
    fixed_SNP_num <- j
    repeat_num <- k
    
    temp_matrix <- matrix(data = NA, nrow = 1, ncol = length(GEBVs$gblupRR)+4)
    temp_matrix[1, 1] <- fixed_SNP_num  
    temp_matrix[1, 2] <- repeat_num
    temp_matrix[1, 3] <- fixed_SNP_vec
    temp_matrix[1, 4] <- selected_SNP_vec
    temp_matrix[1, 5:ncol(temp_matrix)] <- GEBVs$gblupRR
    
    temp_df <- data.frame(temp_matrix, stringsAsFactors = F); colnames(temp_df) <- c("fixed_SNP_num", "repeat_num", "SNP_num_order", "SNP_name", names(temp.y)); temp_df[, 1:10]
    
    if(j == 1 && k == 1){
      write.table(x = temp_df, file = glue("{save_dir}/{today}_01_GEBVs_RR_BLUP_{testPop}_{traitName}.tsv"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    }else{
      write.table(x = temp_df, file = glue("{save_dir}/{today}_01_GEBVs_RR_BLUP_{testPop}_{traitName}.tsv"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
    } # if
    
    cor(x = res$g.pred[pheno_trait[, 1]], y = pheno_trait[, 2])
    temp_cor_df <- data.frame(fixed_SNP_num, repeat_num, "SNP_num_order" = fixed_SNP_vec, "SNP_name" = selected_SNP_vec, "accuracy_RR_BLUP" = cor(x = res$g.pred[pheno_trait[, 1]], y = pheno_trait[, 2]))
    if(j == 1 && k == 1){
      write.table(x = temp_cor_df, file = glue("{save_dir}/{today}_02_accuracy_RR_BLUP_{testPop}_{traitName}.tsv"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    }else{
      write.table(x = temp_cor_df, file = glue("{save_dir}/{today}_02_accuracy_RR_BLUP_{testPop}_{traitName}.tsv"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
    } # if
    
    # GBLUP-GAUSS (RKHS)
    print(glue("RKHS training start"))
    res <- kinship.BLUP(y=temp.y, G.train=Zg, G.pred = g_pre, K.method="GAUSS", X = Zg[, fixed_SNP])
    GEBVs$RKHS <- res$g.train + as.numeric(res$beta)
    # Models$RKHS <- res # RKHS
    fixed_SNP_num <- j
    repeat_num <- k
    
    temp_matrix <- matrix(data = NA, nrow = 1, ncol = length(GEBVs$RKHS)+4)
    temp_matrix[1, 1] <- fixed_SNP_num
    temp_matrix[1, 2] <- repeat_num
    temp_matrix[1, 3] <- fixed_SNP_vec
    temp_matrix[1, 4] <- selected_SNP_vec
    temp_matrix[1, 5:ncol(temp_matrix)] <- GEBVs$RKHS

    temp_df <- data.frame(temp_matrix, stringsAsFactors = F); colnames(temp_df) <- c("fixed_SNP_num", "repeat_num", "SNP_num_order", "SNP_name", names(temp.y)); temp_df[, 1:10]

    
    if(j == 1 && k == 1){
      write.table(x = temp_df, file = glue("{save_dir}/{today}_01_GEBVs_RKHS_{testPop}_{traitName}.tsv"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    }else{
      write.table(x = temp_df, file = glue("{save_dir}/{today}_01_GEBVs_RKHS_{testPop}_{traitName}.tsv"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
    }

    
    cor(x = res$g.pred[pheno_trait[, 1]], y = pheno_trait[, 2])
    temp_cor_df <- data.frame(fixed_SNP_num, repeat_num, "SNP_num_order" = fixed_SNP_vec, "SNP_name" = selected_SNP_vec, "accuracy_RKHS" = cor(x = res$g.pred[pheno_trait[, 1]], y = pheno_trait[, 2]))
    if(j == 1 && k == 1){
      write.table(x = temp_cor_df, file = glue("{save_dir}/{today}_02_accuracy_RKHS_{testPop}_{traitName}.tsv"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    }else{
      write.table(x = temp_cor_df, file = glue("{save_dir}/{today}_02_accuracy_RKHS_{testPop}_{traitName}.tsv"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
    }
  }#for K
} # j

save.image(file = glue("saveImage.RData"))
