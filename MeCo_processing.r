# Hi 👋, this is our code for computing the MeCo refined scores. 
# Please feel free to give us feedbacks! 

-------------------------------------------------------------------------------------------------------------------------------------------

library(readxl)
library(edgeR)

BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)

soft <- read_excel("3D_NF_so.xlsx")
stiff <- read_excel("3D_NF_st.xlsx")

# convert Ensembl_gene_ID in gene symbols
geneID_soft <- ensembldb::select(EnsDb.Hsapiens.v86, keys= soft$Ensembl_gene_ID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneID_stiff <- ensembldb::select(EnsDb.Hsapiens.v86, keys= stiff$Ensembl_gene_ID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

`%notin%` <- Negate(`%in%`)

# remove all the Ensembl_gene_IDs for which there isn't the corresponding gene symbol
soft <- soft[-c(which(soft$Ensembl_gene_ID %notin% geneID_soft$GENEID)),]
stiff <- stiff[-c(which(stiff$Ensembl_gene_ID %notin% geneID_stiff$GENEID)),]

# change Ensembl_gene_ID column with the correspoding gene symbols and rename the column
soft$Ensembl_gene_ID <- geneID_soft$SYMBOL
stiff$Ensembl_gene_ID <- geneID_stiff$SYMBOL

colnames(soft) <- c('Gene_symbol', colnames(soft)[-1])
colnames(stiff) <- c('Gene_symbol', colnames(stiff)[-1])

# Remove all the genes we are not interested in: ribosomial protein genes, mitochondrial genes and others
soft <- soft[-(which(startsWith(soft$Gene_symbol, 'RP'))),]
stiff <- stiff[-(which(startsWith(stiff$Gene_symbol, 'RP'))),]
soft <- soft[-(which(startsWith(soft$Gene_symbol, 'MRPL'))),]
stiff <- stiff[-(which(startsWith(stiff$Gene_symbol, 'MRPL'))),]

soft <- soft[-(which(startsWith(soft$Gene_symbol, 'LOC'))),]
soft <- soft[-(which(startsWith(soft$Gene_symbol, 'LINC'))),]
soft <- soft[-(which(startsWith(soft$Gene_symbol, 'MIR'))),]
soft <- soft[-(which(startsWith(soft$Gene_symbol, 'SNORD'))),]

stiff <- stiff[-(which(startsWith(stiff$Gene_symbol, 'LOC'))),]
stiff <- stiff[-(which(startsWith(stiff$Gene_symbol, 'LINC'))),]
stiff <- stiff[-(which(startsWith(stiff$Gene_symbol, 'MIR'))),]
stiff <- stiff[-(which(startsWith(stiff$Gene_symbol, 'SNORD'))),]

# DE analysis
d <- data.frame(soft[,c(5, 6, 7)], stiff[, c(5, 6, 7)])
countData <- as.matrix(d)
rownames(countData) <- soft$Gene_symbol

y <- DGEList(counts=countData)
y

rownames(y$samples) <- as.factor(c('soft1', 'soft2', 'soft3', 'stiff1', 'stiff2', 'stiff3'))

group <- as.factor(c('soft', 'soft', 'soft', 'stiff', 'stiff', 'stiff'))
y$samples$group <- group
y

y$samples$dex <- as.factor(c(rep(control, 3), rep(treatment, 3)))
y

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs, keep.lib.sizes=FALSE]
dim(y)
# after filtering out those genes with very few counts in the different sample, 15027 genes remain

y <- calcNormFactors(y, method = "TMM")
y

design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

qlf <- glmQLFTest(fit, contrast = c(-1,1))

summary(decideTests(qlf, p.value=0.01, lfc=1))
# 397 up in stiff, 338 up in soft, 14292 don't change

table <- topTags(qlf, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

table <- as.data.frame(table)

table <- table[which((table$logFC > 1 | table$logFC < -1) & table$FDR <0.01),]

table <- cbind(table, upSTIFF = "", upSOFT = "")
for (i in 1:nrow(table)) {
  if (table[i,]$logFC > 1) {table[i,]$upSTIFF <- rownames(table)[i]} else
  {table[i,]$upSOFT <- rownames(table)[i]}
}
head(table)

table <- table[-(which(table$logCPM<0)),]


# Load prototypical information about the patients and their gene expression profile
pheno_data <- TCGA.KIRC.GDC_phenotype
geno_data <- TCGA.KIRC.htseq_fpkm

# Define the sequenced patients from the genomic data set
sequenced_patients <- colnames(geno_data)
sequenced_patients_phenoinfo <- intersect(gsub('-','.',pheno_data$submitter_id.samples),sequenced_patients)

# The geno_data and pheno_data datasets share information about 607 patients
# We store the phenotype of only the shared patients
pheno_data <-pheno_data[which(gsub('-','.',pheno_data$submitter_id.samples) %in% sequenced_patients_phenoinfo),]

# Also we want to keep only the data of Primary Tumor samples. 
pheno_data <- pheno_data[-(which(pheno_data$sample_type.samples != 'Primary Tumor')),]

which(pheno_data$tumor_stage.diagnoses=='not reported') 
submitter_not_reported <- pheno_data$submitter_id.samples[which(pheno_data$tumor_stage.diagnoses=='not reported')]
pheno_data <- pheno_data[-c(115, 423, 504),] # Those are the index of the patients that we want to discard. 
geno_data[,'TCGA.B4.5838.01A'] <- NULL 
geno_data[,'TCGA.MM.A563.01A'] <- NULL 
geno_data[,'TCGA.BP.4798.01A'] <- NULL 

# Updating all after polishing geno and pheno datasets 
patients_phenoinfo <- intersect(gsub('-','.',pheno_data$submitter_id.samples),sequenced_patients)

# find the patients for which we have the expression profile but not prototypical information and remove them from 
# the geno dataset
only_pheno <- colnames(geno_data)[-1][which(colnames(geno_data)[-1] %notin% patients_phenoinfo)]

for (i in only_pheno) {
  geno_data[,i] <- NULL
}

# convert the gene Ensembl_IDs in geno_data in the corresponding gene symbols and remove all the genes we are not interested in
geno_data$Ensembl_ID <- gsub("\\..*","",geno_data$Ensembl_ID)
geneID_geno <- ensembldb::select(EnsDb.Hsapiens.v86, keys= geno_data$Ensembl_ID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geno_data <- geno_data[-c(which(geno_data$Ensembl_ID %notin% geneID_geno$GENEID)),]
geno_data$Ensembl_ID <- geneID_geno$SYMBOL

colnames(geno_data) <- c('Gene_symbol', colnames(geno_data)[-1])

geno_data <- geno_data[-(which(startsWith(geno_data$Gene_symbol, 'RP'))),]
geno_data <- geno_data[-(which(startsWith(geno_data$Gene_symbol, 'MRPL'))),]
geno_data <- geno_data[-(which(startsWith(geno_data$Gene_symbol, 'LOC'))),]
geno_data <- geno_data[-(which(startsWith(geno_data$Gene_symbol, 'LINC'))),]
geno_data <- geno_data[-(which(startsWith(geno_data$Gene_symbol, 'MIR'))),]
geno_data <- geno_data[-(which(startsWith(geno_data$Gene_symbol, 'SNORD'))),]

# extract the lists of genes significantly up or down regulated
geno_data_down <- table$upSOFT[table$upSOFT != ""]
geno_data_up <- table$upSTIFF[table$upSTIFF != ""]

# compute the log2 transformation of all the expression values 
col1 <- geno_data[,1]
geno_data2 <- log2(geno_data[,-1]+1)
geno_data2 <- cbind(col1, geno_data2)

# compute the MeCo score for each patient as the difference between his mean expression value of the 
# up regulated genes and his mean expression value of the down regulated ones
MecoScore <- data.frame(patients = gsub('-','.',pheno_data$submitter_id.samples))
MecoScore$MeSc <- NA

for (i in MecoScore$patients){
  mean_up <- mean(geno_data2[which(geno_data2[,1] %in% geno_data_up),i])
  mean_down <- mean(geno_data2[which(geno_data2[,1] %in% geno_data_down),i])
  MecoScore[which(MecoScore$patients == i),]$MeSc <- mean_up - mean_down
}

pheno_data$submitter_id.samples <- gsub('-','.',pheno_data$submitter_id.samples)


# Refining the MeCo
# Retrieving the Up and Down regulated genes and saving them in a excel file.
library(writexl)

write_xlsx(data.frame(geno_data_up), 'Up regulated genes in Stiff.xlsx')
write_xlsx(data.frame(geno_data_down), 'Down regulated genes in Stiff.xlsx')

# Working on Up regulated genes in STIFF
# Loading the GO_list obtained through Metascape
GOStiff <- read.csv('GO_AllLists_Stiff.csv')

# Filter out annotations with a LogP-value higher than -10
GOStiff <- GOStiff[GOStiff$LogP < -10, ]

library(tidyr)
# We need to split the Hit column into sub_columns, one for each gene. 
columns <- paste0("M_",seq(1, 50))
GOStiff <- separate(GOStiff, Hits, columns)

# At this point we can define refined MeCos based on common GOParental therms: 
# Having an idea of unique GO parent therms
GO_parentals <- (GOStiff$PARENT_GO)

genes_col <-  grep('M_', colnames(GOStiff))

for (GOparental in GO_parentals){
    GO_df <- GOStiff[GOStiff$PARENT_GO == GOparental,] 
    gene_list <- NA
    for (row in seq(1, nrow(GO_df))){
        for (gene in genes_col){
          ifelse(GO_df[row, gene] %in% gene_list, next, gene_list <- c(gene_list, (GO_df[row, gene]))) 
        }
    } 
    assign(paste0("GLT_",GOparental), gene_list)
}

# Working on Down regulated genes in STIFF (so Up regulated genes in SOFT)

#Loading the GO_list obtained through Metascape
GOSoft <- read.csv('GO_AllList_Soft.csv')

# Filter out annotations with a LogP-value higher than -10
GOSoft <- GOSoft[GOSoft$LogP < -10, ]

# We need to split the Hit column into sub_columns, one for each gene. 
columns <- paste0("M_",seq(1, 50))
GOSoft <- separate(GOSoft, Hits, columns)

# At this point we can define refined MeCos based on common GOParental therms: 
# Having an idea of unique GO parent therms
GO_parentals <- (GOSoft$PARENT_GO)

genes_col <-  grep('M_', colnames(GOSoft))

for (GOparental in GO_parentals){
  GO_df <- GOSoft[GOSoft$PARENT_GO == GOparental,] 
  gene_list <- NA
  for (row in seq(1, nrow(GO_df))){
    for (gene in genes_col){
      ifelse(GO_df[row, gene] %in% gene_list, next, gene_list <- c(gene_list, (GO_df[row, gene]))) 
    }
  } 
  assign(paste0("GLO_",GOparental), gene_list)
}

# Computing the MeCo refined Stimulus

MecoScore$MR_Stimulus <- NA

for (i in MecoScore$patients){
  mean_up <- mean(geno_data2[which(geno_data2[,1] %in% `GLT_19_GO:0050896 response to stimulus`),i])
  mean_down <- mean(geno_data2[which(geno_data2[,1] %in% `GLO_19_GO:0050896 response to stimulus`),i])
  MecoScore[which(MecoScore$patients == i),]$MR_Stimulus <- mean_up - mean_down
}


# Computing the MeCo refined Regulation

MecoScore$MR_Regulation <- NA

for (i in MecoScore$patients){
  mean_up <- mean(geno_data2[which(geno_data2[,1] %in% `GLT_19_GO:0050789 regulation of biological process`),i])
  mean_down <- mean(geno_data2[which(geno_data2[,1] %in% `GLO_19_GO:0050789 regulation of biological process`),i])
  MecoScore[which(MecoScore$patients == i),]$MR_Regulation <- mean_up - mean_down
}


# Computing the MeCo refined Developmental

MecoScore$MR_Developmental <- NA

for (i in MecoScore$patients){
  mean_up <- mean(geno_data2[which(geno_data2[,1] %in% `GLT_19_GO:0032502 developmental process`),i])
  mean_down <- mean(geno_data2[which(geno_data2[,1] %in% `GLO_19_GO:0032502 developmental process`),i])
  MecoScore[which(MecoScore$patients == i),]$MR_Developmental <- mean_up - mean_down
}



# DATASET DEFINITION

# set up a dataset with the variables of interest:
surv_data <- data.frame(pheno_data$submitter_id.samples) 
#age
surv_data$age <- pheno_data$age_at_initial_pathologic_diagnosis
#time to event
surv_data$time <- NA
time_to_death <- pheno_data$days_to_death.demographic
time_to_last_FU <- pheno_data$days_to_last_follow_up.diagnoses

for (i in 1:length(pheno_data$submitter_id.samples)){
  if (is.na(time_to_death[i]))
    surv_data$time[i] <- time_to_last_FU[i]
  else
    surv_data$time[i] <- time_to_death[i]
}

#status
surv_data$status <- as.numeric(as.factor(pheno_data$vital_status.demographic)) # 1 = patient alive, 2 = patient dead
#metastasis
surv_data$metastasis <- NA
meta <- pheno_data$additional_surgery_metastatic_procedure

for (i in 1:length(pheno_data$submitter_id.samples)){
  if (meta[i]=="")
    surv_data$metastasis[i] <- 'NO'
  else
    surv_data$metastasis[i] <- 'YES'
}
surv_data$metastasis <- as.factor(surv_data$metastasis)
#sex
surv_data$gender <- as.factor(pheno_data$gender.demographic)
#tumor stage
surv_data$stage <- as.factor(pheno_data$tumor_stage.diagnoses)
# laterality
surv_data$laterality <- as.factor(pheno_data$laterality)
# MeCo scores
surv_data$MeCo <- MecoScore$MeSc
surv_data$MeCo_st <- MecoScore$MR_Stimulus
surv_data$MeCo_reg <- MecoScore$MR_Regulation
surv_data$MeCo_dev <- MecoScore$MR_Developmental

# check the presence of NA values
which(is.na(surv_data))

# remove time=0
surv_data <- surv_data[-which(surv_data$time==0),]
