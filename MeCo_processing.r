# Hi 👋, this is our code for computing the MeCo refined scores. 
# Please feel free to give us feedbacks! 

-------------------------------------------------------------------------------------------------------------------------------------------
         Differential Gene Expression 
-------------------------------------------------------------------------------------------------------------------------------------------

# Set up the data 💾
    # First of all we need to load the full expression matrix for CAFs (Cancer Associated Fibroblasts) 
    # and NF (Normal Fibroblasts) cultured on 3D_soft and 3D_stiff substrates.
  
        install.packages('readxl')
        library(readxl)

        soft <- read_excel("3D_NF_so.xlsx")
        stiff <- read_excel("3D_NF_st.xlsx")

    # Combine the two expression matrices (soft & stiff) to have a single one with only the expression levels of CAFs 
    # We want to keep the ENSEMBL id for each gene for later use in DESeq2, so we keep the first column. 

        countData <- data.frame(soft[,c(1, 5, 6, 7)], stiff[, c(5, 6, 7)])

    # DESeq2 Analysis 🧬
        # Metadata and preparation of the data
          # Create a matrix for the metadata associated to each sample
          # This will enable DESeq2 to recognize which sample is associated which condition, for simplicity reported here as 'control' and 'treatment'.
          # This way we can compare the gene expression of CAF grown on 3D_soft substrate (control) 
          # with the gene expression of CAF grown on 3D_stiff (treatment).

            control <- 'control'
            treatment <- 'treatment'
            id <- c('soft1', 'soft2', 'soft3', 'stiff1', 'stiff2', 'stiff3')
            dex <- c(rep(control, 3), rep(treatment, 3))

            metadata <-  data.frame(id=as.factor(id), dex=as.factor(dex))
            metadata
        
        # Performing DGE and taking care of FDR(False Discovery Rate) 🔬
          # We use DESeq2 to obtain differential expressed genes. We used a threshold of p.adjusted <= 0.01 and log2FC >= |±3|
          # Opting for such a high FC is needed because we have only 3 replicates for each condition. While the p.adjusted
          # is used for selecting only highly significant differentially expressed genes. 
          # Discard all the lines that have at least two cells with 9 or less counts, this will help us with the FDR.
          # For FDR we used the Benjamini-Yosef methods.
            
            BiocManager::install('DESeq2')
            library(DESeq2)
            
            dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metadata, 
                                design= ~ dex, tidy = TRUE)
            keep <- rowSums(counts(dds) >= 10) >= 2
            dds <- dds[keep,]
            dds <- DESeq(dds)
            res <- results(dds,alpha = 0.05, pAdjustMethod = 'BH')
            res 

          # We store the results in a data frame and label the differentially expressed genes with UP or DOWN, 
          # if the gene is up regulated or down regulated, respectively.
            
            de <- data.frame(res)
            de$diffexpressed <- "NO"
            de$diffexpressed[de$log2FoldChange > 2 & de$padj < 0.05] <- "UP" 
            de$diffexpressed[de$log2FoldChange < -2 & de$padj < 0.05] <- "DOWN"

    # Vulcano Plot 🌋
        # To visualize the results we use a VulcanPlot 
            if (!require("pacman")) install.packages("pacman")
                pacman::p_load(here,  
                        tidyverse, 
                        janitor,     # Cleaning column names  
                        scales,      # Transform axis scales   
                        ggrepel)     # Optimise plot label separation  

        # Polishing the data for later use…            
            dati <- de %>%
                    mutate(gene_type = case_when(log2FoldChange >= 2 &  padj <= 0.05 ~ "up", 
                                                log2FoldChange <= -2 &  padj <= 0.05 ~ "down", 
                                                TRUE ~ "ns"))   
            count(dati, gene_type)


            cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
            sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
            alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

        # Using ggplot2 for the VulcanPLot
            v <- ggplot(data = dati, aes(x = log2FoldChange,
                    y = -log10(padj),
                    fill = gene_type,    
                    size = gene_type,
                    alpha = gene_type)) + 
                geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
                            colour = "black") + 
                geom_hline(yintercept = -log10(0.05),
                            linetype = "dashed") + 
                
                geom_vline(xintercept = c(2, -2),
                            linetype = "dashed") +
                scale_fill_manual(values = cols) + # Modify point colour
                scale_size_manual(values = sizes) + # Modify point size
                scale_alpha_manual(values = alphas) + # Modify point transparency
                scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                                    limits = c(-10, 10))+
                ylim(0,125)
            v

-------------------------------------------------------------------------------------------------------------------------------------------
        Working on TGCA-KIRC Project data
-------------------------------------------------------------------------------------------------------------------------------------------

# Polishing and checking the data 🔎
    # Now we want to prepare the TGCA-KIRC data for later use in computing the MeCo refined scores

    # Load the TCGA-KIRC.htseq_fpkm dataset. For each patient we have the expression level of almost all the genes.
    # Also load the phenotype dataset of the TCGA-KIRC project. 

        geno_data <- read.delim("TCGA-KIRC.htseq_fpkm.tsv")     #importing manually the dataset may be more convenient due to its size 
        pheno_data <- read.delim("TCGA-KIRC.GDC_phenotype.tsv")
        
    # Define the sequenced patients from the genomic data set
        sequenced_patients <- colnames(geno_data)
        all <- intersect(gsub('-','.',pheno_data$submitter_id.samples),sequenced_patients)
        all[seq(1,10)]

    # The geno_data and pheno_data datasets share information about 607 patients
    # We store the phenotype of only the shared patients
        pheno_data <-pheno_data[which(gsub('-','.',pheno_data$submitter_id.samples) %in% all),]

    # Also we want to keep only the data of Primary Tumor samples. 
        pheno_data <- pheno_data[-(which(pheno_data$sample_type.samples != 'Primary Tumor')),]

        # Double checking
            which(pheno_data$sample_type.samples != 'Primary Tumor')

    # For some patients may not be reported the stage. 
    # We must discard them for the final dataset. 
        which(pheno_data$tumor_stage.diagnoses=='not reported') 
        submitter_not_reported <- pheno_data$submitter_id.samples[which(pheno_data$tumor_stage.diagnoses=='not reported')]
        pheno_data <- pheno_data[-c(115, 423, 504),] # Those are the index of the patients that we want to discard. 
        geno_data[,'TCGA.B4.5838.01A'] <- NULL 
        geno_data[,'TCGA.MM.A563.01A'] <- NULL 
        geno_data[,'TCGA.BP.4798.01A'] <- NULL 

    # Updating all after polishing geno and pheno datasets 
        all <- intersect(gsub('-','.',pheno_data$submitter_id.samples),sequenced_patients)

    # We must remove the version number from the Ensembl.id for later use 
        Ensembl_id <- geno_data[,1]
        Ensembl_id <- gsub("\\.[0-9]*$", "", Ensembl_id)

    # Finally update the geno_data with only the patient id of primary tumoral samples, with a defined stage and keeping the Ensembl id
        geno_data <- data.frame(Ensembl_id, geno_data[,which(colnames(geno_data) %in% all)])
        sequenced_patients <- colnames(geno_data)[-1]

-------------------------------------------------------------------------------------------------------------------------------------------
        Calculating the MeCo Score 
-------------------------------------------------------------------------------------------------------------------------------------------

# We first propose a general MeCo score 📊
    # This way of computing the MeCo does not take in consideration the impact of each gene on a specific pathway. 
    # So it is quite general and later will be refined using Pathway Analysis. 
    # It is useful to get a general understanding of the mechanical conditioning of the cell due to the substrate stiffness. 

# UP and DOWN regulated genes 🔺 🔻
      # First of all we must create a list of UP and DOWN regulated genes using the results of DGE. 
    
        UP <- rownames(de[which(de$diffexpressed == 'UP'),])
        DOWN <- rownames(de[which(de$diffexpressed == 'DOWN'),])

      # Than we use the lists to obtain two distinct expression matrix for genes UP and DOWN regulated. 
      # Both matrix are linked by the patient id. 
    
        geno_data_down <- geno_data[which(geno_data$Ensembl_id %in% DOWN), ]
        geno_data_up <- geno_data[which(geno_data$Ensembl_id %in% UP), ]

# Finally calculating the genral MeCo score for each patient ⭐️
        
        MecoScore <- data.frame(sequenced_patients) 
        MecoScore$MeSc <- NA

        for (i in 1:nrow(MecoScore)){
            MecoScore$MeSc[i] <- mean(as.numeric(geno_data_up[,i+1])) - mean(as.numeric(geno_data_down[,i+1]))
        }

        summary(MecoScore)

# Computing the general MeCo score for each stage 💻
        
        # Check presence of stage i patients in the genomic dataset
            # takes the first column containing the submitter id sample
                samples_stage1 <- pheno_data[which(pheno_data$tumor_stage.diagnoses=='stage i'),1] 
            # replace - with . as patients are reported differently in the two datasets
                samples_stage1 <- gsub('-','.', samples_stage1)
            # get the list of patients for which we have stage 1 and gene expression data
                sequenced_stage1 <- intersect(sequenced_patients,samples_stage1)

        # Check presence of stage ii patients in the genomic dataset
            samples_stage2 <- pheno_data[which(pheno_data$tumor_stage.diagnoses=='stage ii'),1]
            samples_stage2 <- gsub('-','.', samples_stage2)
            sequenced_stage2 <- intersect(sequenced_patients,samples_stage2)

        # Check presence of stage iii patients in the genomic dataset
            samples_stage3 <- pheno_data[which(pheno_data$tumor_stage.diagnoses=='stage iii'),1]
            samples_stage3 <- gsub('-','.', samples_stage3)
            sequenced_stage3 <- intersect(sequenced_patients,samples_stage3)

        # Check presence of stage iv patients in the genomic dataset
            samples_stage4 <- pheno_data[which(pheno_data$tumor_stage.diagnoses=='stage iv'),1]
            samples_stage4 <- gsub('-','.', samples_stage4)
            sequenced_stage4 <- intersect(sequenced_patients,samples_stage4)

        # General MeCo for: 
        # stage i
            mean(MecoScore[which(MecoScore$sequenced_patients %in% sequenced_stage1), 'MeSc'])

        # stage ii
            mean(MecoScore[which(MecoScore$sequenced_patients %in% sequenced_stage2), 'MeSc'])

        # stage iii
            mean(MecoScore[which(MecoScore$sequenced_patients %in% sequenced_stage3), 'MeSc'])

        # stage iv
            mean(MecoScore[which(MecoScore$sequenced_patients %in% sequenced_stage4), 'MeSc'])

-------------------------------------------------------------------------------------------------------------------------------------------
        Refining the Meco
-------------------------------------------------------------------------------------------------------------------------------------------

# We performed pathway analysis with Metascape. 
# Given a list of genes, Metascape analyzes in which pathway each gene is involved and returns the most influenced pathways. 
    # Metascape 📈
        # We need an excel file with all the differentialy expressed genes
              
            install.packages("writexl")
            library(writexl)

            DifferentialExpressedGenes <- c(DOWN,UP) 
            write_xlsx(data.frame(DifferentialExpressedGenes), 'PathwayAnalysis.xlsx')

        # At this point we must use Metascape (https://metascape.org/gp/index.html)
        # And save the results in a excel file for later use. 
    
    # Manually curated Pathways 🛣
        # Since Metascape gives us a huge number of pathways in which the differentialy expressed genes are involved,
        # we need to manually curate the principal pathways in which we are interested. 
        # We chose to use the 5 most representative pathways in our sample of CAF and 
        # that are strictly involved in tumoral progression. 

        # So we focused our attention on: 
        # • ECM (ExtraCellular Matrix)
        # • Proliferation
        # • Chemotaxis
        # • Inflammation
        # • Antitumoral mechanism 

        # Get all the distinct pathways 
            metascape_result <- read_excel("metascape_result.xlsx") # PLease, insert the name of your file! 
            distinct_pathways <- colnames(metascape_result)[-seq(1, 15)]
            head(distinct_pathways)

        # For each pathway get the list of genes.
        # We store this information in a List of pathways made of Lists of genes (List of Lists, aka LoL) 

            sum(as.integer(metascape_result$`R-HSA-1474244 Extracellular matrix organizat`))

            LoL <- list()

            for(i in seq(1, length(distinct_pathways))){
            LoL[i] <- metascape_result[which(metascape_result[,(15 + i)] == '1.0'), 'DifferentialExpressedGenes']
            }

            names(LoL) <- distinct_pathways

        # manually curate the main pathways (yeah! more fun!)
            ECM_path <- unique(c(LoL$`M5884 NABA CORE MATRISOME`,
                                LoL$`R-HSA-1474244 Extracellular matrix organizat`, 
                                LoL$`GO:0030198 extracellular matrix organizat`, 
                                LoL$`M5882 NABA PROTEOGLYCANS`, 
                                LoL$`R-HSA-1474228 Degradation of the extracellul`, 
                                LoL$`R-HSA-2129379 Molecules associated with elas`))
                    
            Proliferation_path <- unique(c(LoL$`GO:0001501 skeletal system development`, 
                                        LoL$`GO:0061061 muscle structure development`))

            Chemotaxis_path <- LoL$`GO:0006935 chemotaxis`

            Inflammatory_path <- LoL$`GO:0006954 inflammatory response`

            Antitumoral_path <- LoL$`GO:0008285 negative regulation of cell po`
            
        # For each pathway subdivided the genes in UP and DOWN regulated according to the
        # Differential Gene Expression analysis on CAFs 
            ECM_path_Down <- intersect(ECM_path, DOWN)
            ECM_path_Up <- intersect(ECM_path, UP)

            Proliferation_path_Down <- intersect(Proliferation_path, DOWN)
            Proliferation_path_Up <- intersect(Proliferation_path, UP)

            Chemotaxis_path_Down <- intersect(Chemotaxis_path, DOWN)
            Chemotaxis_path_Up <- intersect(Chemotaxis_path, UP)

            Inflammatory_path_Down <- intersect(Inflammatory_path, DOWN)
            Inflammatory_path_Up <- intersect(Inflammatory_path, UP)

            Antitumoral_path_Down <- intersect(Antitumoral_path, DOWN)
            Antitumoral_path_Up <- intersect(Antitumoral_path, UP)

    # Computing the MeCo refined scores ⭐️⭐️⭐️
        # We call it refined because it is pathway based, 
        # so it is more easily interpretable and also more meaningful for assessing the role of 
        # mechanotransduction on tumoral cells. 
        
        # ECM MeCo refined
            geno_data_ECMdown <- geno_data[which(geno_data$Ensembl_id %in% ECM_path_Down), ]
            geno_data_ECMup <- geno_data[which(geno_data$Ensembl_id %in% ECM_path_Up), ]

            MecoScoreECM <- data.frame(sequenced_patients) 
            MecoScoreECM$MeSc <- NA

            for (i in 1:nrow(MecoScoreECM)){
            MecoScoreECM$MeSc[i] <- mean(as.numeric(geno_data_ECMup[,i+1])) - mean(as.numeric(geno_data_ECMdown[,i+1]))
            }

            summary(MecoScoreECM)
    
        # Proliferation MeCo refined
            geno_data_down <- geno_data[which(geno_data$Ensembl_id %in% Proliferation_path_Down), ]
            geno_data_up <- geno_data[which(geno_data$Ensembl_id %in% Proliferation_path_Up), ]

            MecoScorePro <- data.frame(sequenced_patients) 
            MecoScorePro$MeSc <- NA

            for (i in 1:nrow(MecoScorePro)){
            MecoScorePro$MeSc[i] <- mean(as.numeric(geno_data_up[,i+1])) - mean(as.numeric(geno_data_down[,i+1]))
            }

            summary(MecoScorePro)

        # Chemotaxis MeCo refined
            geno_data_down <- geno_data[which(geno_data$Ensembl_id %in% Chemotaxis_path_Down), ]
            geno_data_up <- geno_data[which(geno_data$Ensembl_id %in% Chemotaxis_path_Up), ]

            MecoScoreCh <- data.frame(sequenced_patients) 
            MecoScoreCh$MeSc <- NA

            for (i in 1:nrow(MecoScoreCh)){
            MecoScoreCh$MeSc[i] <- mean(as.numeric(geno_data_up[,i+1])) - mean(as.numeric(geno_data_down[,i+1]))
            }

            summary(MecoScoreCh)

        # Inflammation MeCo refined

            geno_data_down <- geno_data[which(geno_data$Ensembl_id %in% Inflammatory_path_Down), ]
            geno_data_up <- geno_data[which(geno_data$Ensembl_id %in% Inflammatory_path_Up), ]

            MecoScoreInf <- data.frame(sequenced_patients) 
            MecoScoreInf$MeSc <- NA

            for (i in 1:nrow(MecoScoreInf)){
            MecoScoreInf$MeSc[i] <- mean(as.numeric(geno_data_up[,i+1])) - mean(as.numeric(geno_data_down[,i+1]))
            }

            summary(MecoScoreInf)

        # Antitumoral mechanism MeCo refined 
            geno_data_down <- geno_data[which(geno_data$Ensembl_id %in% Antitumoral_path_Down), ]
            geno_data_up <- geno_data[which(geno_data$Ensembl_id %in% Antitumoral_path_Up), ]

            MecoScoreAnt <- data.frame(sequenced_patients) 
            MecoScoreAnt$MeSc <- NA

            for (i in 1:nrow(MecoScoreAnt)){
            MecoScoreAnt$MeSc[i] <- mean(as.numeric(geno_data_up[,i+1])) - mean(as.numeric(geno_data_down[,i+1]))
            }

            summary(MecoScoreAnt)
-------------------------------------------------------------------------------------------------------------------------------------------
