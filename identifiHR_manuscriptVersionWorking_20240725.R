# background --------------------------------------------------------------

# feature engineering for IdentifiHR

# Pooling features of gene expression to improve the generalisability of the IdentifiHR classifier.
# using normalised gene expression values, and an elastic net logistic regression model did not support generalisability.

# loadPackages ------------------------------------------------------------

library(glmnet)
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(caret)
library(foreach)
library(limma)
library(edgeR)
library(readxl)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(ggridges)
library(glmnet)
library(ggrepel)

# loadData ----------------------------------------------------------------

# tcgaData ----------------------------------------------------------------

# TCGA counts
tcgaCounts <- read_csv("~/CountMatrix_TCGA_OC_20220826.csv")
dim(tcgaCounts)
# [1] 60660   382

# TCGA clinical
tcgaClin <- read.csv("~/identifiHR/benchmarking/predictingHrInOtherCancers/tcga-ov_clinical_20231109.csv")
addTcgaClin <- read.csv("~/TCGA/clinicalTcgaCompleteData_20230619.csv") %>% 
  dplyr::select( "Sample.ID.1", "Case.ID", "subtypeJoin" , "tp53Mutcomplete", "germlineBRCA2", "somaticBRCA1", "somaticBRCA2", "brcaStatus", "Germline.or.somatic.mutation.in.HR.gene",
                 "TSS", "participant", "sampleVial","portionAnalyte", "plate", "center") %>% 
  dplyr::rename(Sample_Id = Sample.ID.1)
tcgaClin <- left_join(tcgaClin, addTcgaClin, by = "Sample_Id") %>% 
  dplyr::filter(!duplicated(Sample_Id))
# write_csv(tcgaClin, "~/identifiHR/identifiHR_tcgaClinicalComplete_20230304.csv")

# aocsData ----------------------------------------------------------------

# Between 500 ng and 2 μg of total RNA (RNA integrity number >7, except for autopsy samples) was used for library preparation using the TruSeq RNA Sample Preparation v2 kit (Illumina) as per the manufacturer’s instructions using the Low-Throughput protocol. 
# All libraries were sequenced as paired-end 100-bp on a HiSeq2000 instrument (Illumina).
# **Raw data uploaded to European Genome-phenome Archive (EGA) repository under the accession code EGAD00001000877**
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

# AOCS counts
countsAocs <- read.delim("~/identifiHR/validation/identifiHR_PMID36456881_longTermSurv/GSE209964_analysis_14_08_2019.HTSeq.all.raw.txt")

# AOCS clinical
clinAocs <- read_csv("~/identifiHR/validation/identifiHR_PMID36456881_longTermSurv/GSE209964_series_matrix.csv")
chord <- read_excel("~/identifiHR/validation/identifiHR_PMID36456881_longTermSurv/NIHMS1894596-supplement-S__Tables_3.xlsx")
mutAocs <- read_excel("~/identifiHR/validation/identifiHR_PMID36456881_longTermSurv/NIHMS1894596-supplement-S__Tables_5.xlsx")

# wranglingTcga -----------------------------------------------------------

# run this regardless of the sample annotation section (ie. even if using the alternateData...)
colnames(tcgaCounts)[1] <- "X"
tcgaCounts$X <- gsub("\\..*","", tcgaCounts$X) # remove last part of ensembl ID
tcgaCounts <- tcgaCounts[!duplicated(tcgaCounts$X), ]
tcgaCounts <- tcgaCounts %>%
  column_to_rownames(var = "X")

colnames(tcgaCounts) <- str_replace_all(colnames(tcgaCounts), "\\.", "-") # correct sampleIDs
colnames(tcgaCounts) <- str_sub(colnames(tcgaCounts), start = 1, end = 15)

tcgaCountsT <- t(tcgaCounts)

# subset counts nad clinical to only cases that have complete HR status
dim(tcgaCounts)
# [1] 60616   381
tcgaCases <- colnames(tcgaCounts)
tcgaClinSub <- tcgaClin[tcgaClin$Sample_Id %in% tcgaCases, ]
rownames(tcgaClinSub) <- tcgaClinSub$Sample_Id
tcgaClinSub <- tcgaClinSub[tcgaCases, ]
tcgaCountsFil <- tcgaCounts[ ,colnames(tcgaCounts) %in% rownames(tcgaClinSub)]
tcgaCountsFilVec <- colnames(tcgaCountsFil)
tcgaClinSub <- tcgaClinSub[tcgaClinSub$Sample_Id %in% tcgaCountsFilVec, ]

dim(tcgaCountsFil)
# [1] 60616   361
dim(tcgaClinSub)
#[1] 361 182

table(colnames(tcgaCountsFil) == rownames(tcgaClinSub))
# n = 361 TCGA cohort

# get HR status labels for TCGA
tcgaClinClean <- tcgaClinSub %>% 
  dplyr::filter(!is.na(Sample_Id)) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Sample_Id") 

tcgaLabels <- tcgaClinSub %>% 
  dplyr::filter(!is.na(Sample_Id)) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Sample_Id") %>% 
  dplyr::select(hrStatus)

table(colnames(tcgaCountsFil) == rownames(tcgaLabels))
# TRUE 
# 361 

table(tcgaLabels$hrStatus)
# HRD HRP 
# 194 167

# wranglingAocs -----------------------------------------------------------

rownames(countsAocs) <- countsAocs$GeneID
clinAocsT <- t(clinAocs)

colnames(clinAocsT) <- clinAocsT[1,]
clinAocsT <- clinAocsT[-1,] %>% 
  as.data.frame()
colnames(clinAocsT)[7] <- "tissueSource"
colnames(clinAocsT)[9] <- "CaseId"
colnames(clinAocsT)[12] <- "histotype"
colnames(clinAocsT)[13] <- "runBatch"
colnames(clinAocsT)[14] <- "libraryPrepBatch"
clinAocsT$SampleIdCounts <- rownames(clinAocsT)
clinAocsT <- clinAocsT %>% 
  dplyr::select(1, 7, 9, 12, 13, 14, 47)
clinAocsT$Case_ID <- str_sub(clinAocsT$CaseId, start = 10, end = 17)

chord <- chord %>% 
  mutate(Case_ID = case_when(str_detect(Sample_ID, "AOCS") ~ str_sub(Sample_ID, start = 1, end = 8),
                             TRUE ~ str_sub(Sample_ID, start = 1, end = 9)))

colnames(mutAocs)[1] <- "Case_ID"

chordMutAocs <- left_join(chord, mutAocs, by = "Case_ID") %>% 
  mutate(hrStatus = case_when(Gene == "BRCA1" | Gene == "BRCA2" ~ "HRD",
                              hr_status == "HR_deficient" ~ "HRD",
                              hr_status == "HR_proficient" ~ "HRP")) %>% 
  dplyr::filter(!is.na(hrStatus))

clinHrAocs <- left_join(clinAocsT, chordMutAocs, by = "Case_ID", relationship = "many-to-many") %>% 
  mutate(hrStatus = case_when(tissueSource == "FT - normal control" ~ "HRP",
                              TRUE ~ hrStatus)) %>% 
  dplyr::filter(!is.na(hrStatus)) %>% 
  mutate(Sample = SampleIdCounts)

# filterAocsCounts --------------------------------------------------------

sampleClinAocs <- clinHrAocs$SampleIdCounts
countsAocsFil <- countsAocs[, colnames(countsAocs) %in% sampleClinAocs]
countsAocsFilIdx <- countsAocsFil[ ,sampleClinAocs]
countsAocsEns <- countsAocsFilIdx
rownames(countsAocsEns) <- str_sub(rownames(countsAocsEns), start = 1, end = 15)

# makeAocsLabels ----------------------------------------------------------

aocsLabelsClin <- clinHrAocs %>% 
  dplyr::select(Sample, hrStatus) %>% 
  column_to_rownames(var = "Sample")

table(colnames(countsAocsEns) == rownames(aocsLabelsClin))
# TRUE 
# 99 
aocsLabels <- aocsLabelsClin$hrStatus

# clinHrAocs <- clinHrAocs %>% 
#   dplyr::filter(!is.na(Tumor_cell_puritya)) %>% 
#   mean(Tumor_cell_puritya)
# 
# mean(trainCases$age_at_index)

# subsetCountsToOverlap ---------------------------------------------------

aocsGenes <- rownames(countsAocsEns)
tcgaGenes <- rownames(tcgaCountsFil)

countsAocsEns <- countsAocsEns[rownames(countsAocsEns) %in% tcgaGenes, ]
tcgaCountsFil <- tcgaCountsFil[rownames(tcgaCountsFil) %in% aocsGenes, ]

# splitIntoTrainingTesting ------------------------------------------------

set.seed(NULL)
set.seed(220310)
set.seed(220310)  # For reproducibility
# 220310
#5937
trainIndices <- sample(ncol(tcgaCountsFil), ncol(tcgaCountsFil) * 0.8)  # 80% for training
# write_csv(as.data.frame(trainIndices), "~/identifiHR/trainingAnalysis/identifiHR_trainIndices.csv")

trainData <- tcgaCountsFil[ ,trainIndices]
dim(trainData)
# [1] 52125    288
trainLabels <- tcgaLabels[trainIndices, ]
length(trainLabels)
# [1] 288
trainCases <- tcgaClinClean[colnames(trainData), ]

testData <- tcgaCountsFil[ ,-trainIndices]
dim(testData)
# [1]  52125    73
testLabels <- tcgaLabels[-trainIndices, ]
length(testLabels)
# [1] 73
testCases <- tcgaClinClean[colnames(testData), ]

# hrDivsionInTrainAndTest -------------------------------------------------

# training cases
table(trainLabels)
# trainLabels
# HRD HRP 
# 160 128

mean(trainCases$HRD)
# [1] 47.10069

mean(trainCases$age_at_index)
# [1] 59.91667

stageTrain <- table(trainCases$figo_stage, trainCases$hrStatus)
#               HRD HRP
#   '--          1   2
#   Stage IC     1   0
#   Stage IIA    0   3
#   Stage IIB    2   1
#   Stage IIC    9   3
#   Stage IIIA   3   2
#   Stage IIIB   7   5
#   Stage IIIC 110  92
#   Stage IV    27  20

raceTrain <- table(trainCases$race, trainCases$hrStatus)
#                                           HRD HRP
# american indian or alaska native            1   0
# asian                                       4   4
# black or african american                  12   7
# native hawaiian or other pacific islander   1   0
# not reported                                5   2
# white                                     137 115

min(trainCases$age_at_index)
# 34
max(trainCases$age_at_index)
# 87
# testing cases
table(testLabels)
# testLabels
# HRD HRP 
# 34  39 

mean(testCases$HRD)
# [1] 42.0411

mean(testCases$age_at_index)
# [1] 59.45205

stageTest <- table(testCases$figo_stage, testCases$hrStatus)
#             HRD HRP
# Stage IIC    1   1
# Stage IIIA   1   0
# Stage IIIB   1   1
# Stage IIIC  28  32
# Stage IV     3   5

raceTest <- table(testCases$race, testCases$hrStatus)
#                                   HRD HRP
# american indian or alaska native   1   0
# asian                              2   0
# black or african american          2   4
# not reported                       2   2
# white                             27  33

min(testCases$age_at_index)
# 30
max(testCases$age_at_index)
# 87
# statisticalTesting ------------------------------------------------------

t.test(trainCases$age_at_index, testCases$age_at_index, alternative = "two.sided")
# Welch Two Sample t-test
# data:  trainCases$age_at_index and testCases$age_at_index
# t = 0.30338, df = 107.3, p-value = 0.7622
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.571230  3.500454
# sample estimates:
#   mean of x mean of y 
# 59.91667  59.45205 
chisq.test(stageTrain)
chisq.test(stageTest)
chisq.test(raceTrain)
chisq.test(raceTest)

# training ----------------------------------------------------------------

# annotateGenes -----------------------------------------------------------

ensTcga <- rownames(trainData)
useEns <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")
ensTcgaLoc <- getBM(attributes = c( "ensembl_gene_id","hgnc_symbol","entrezgene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "band"),
                    filters = "ensembl_gene_id",
                    values = ensTcga, 
                    mart = useEns)
colnames(ensTcgaLoc)[1] <- "ENSEMBL"

ensTcgaAnno <- ensTcgaLoc

# filter so no Y or MT chr
# filter so only protein-coding genes are retained

ensKeep <- ensTcgaAnno %>% 
  dplyr::filter(chromosome_name != "Y") %>% 
  dplyr::filter(chromosome_name != "MT") %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  dplyr::filter(!duplicated(ENSEMBL))

# make vector of filtered ENSEMBL IDs for count filtering
ensKeepVec <- ensKeep$ENSEMBL
length(ensKeepVec)
# [1] 19152

# subsetGenesToRemoveNonCodingYMt -----------------------------------------

trainData <- trainData[rownames(trainData) %in% ensKeepVec, ]
trainData <- trainData[ensKeepVec, ]
dim(trainData)
# [1] 19152   288

# filterLowlyExpressedGenesInTraining -------------------------------------

# make the DGEList
dgeTrain <- DGEList(counts = trainData, genes = ensKeep, samples = trainLabels, group = trainCases$hrStatus) 
# Number of rows in 'samples' must equal number of columns in 'counts'

# calculate TMM normalization factors
dgeTrain <- calcNormFactors(dgeTrain)

# get the normalized counts
cpms <- cpm(dgeTrain, log=FALSE)

(min(dgeTrain[["samples"]][["lib.size"]]))/1000000
# [1] 13.73634

cpmsSum <- cbind(cpms, rowSums(cpms))
cpmsSum <- cpmsSum %>% 
  as.data.frame() %>% 
  dplyr::rename(sum = V289) %>% 
  dplyr::filter(sum >= 10/(min(dgeTrain[["samples"]][["lib.size"]])))

notLowlyExpGenes <- rownames(cpmsSum)
length(notLowlyExpGenes)
# [1] 18882

# subset training data counts
countsFilNotLow <- trainData[rownames(trainData) %in% notLowlyExpGenes, ]
ensTrainKeepVec <- rownames(countsFilNotLow)

# subset gene annotations
ensTcgaAnnoFil <- ensKeep[ensKeep$ENSEMBL %in% ensTrainKeepVec, ]
rownames(ensTcgaAnnoFil) <- ensTcgaAnnoFil$ENSEMBL

table(rownames(countsFilNotLow) == ensTcgaAnnoFil$ENSEMBL)
# TRUE 
# 18882  

# deAnalysisForFeatureSelection -------------------------------------------

# make the DGEList
dgeTrainDe <- DGEList(counts = countsFilNotLow, genes = ensTcgaAnnoFil, samples = trainCases, group = trainCases$hrStatus) 
keep <- filterByExpr(dgeTrainDe, group = trainCases$hrStatus) 
dgeTrainDeFil <- dgeTrainDe[keep, ,keep.lib.sizes = FALSE ]
dim(dgeTrainDeFil)
# [1] 15901   288
dgeTrainDeFil <- calcNormFactors(dgeTrainDeFil)
group <-  interaction(dgeTrainDeFil$samples$group)
modelMatrixTcgaOv <- model.matrix(~ 0 + group, na.action = 'na.pass')
colnames(modelMatrixTcgaOv) <- gsub(" ", ".", colnames(modelMatrixTcgaOv))
voomTcgaOv <- voom(dgeTrainDeFil, modelMatrixTcgaOv, plot = TRUE)
fitTcgaOv <- lmFit(voomTcgaOv, modelMatrixTcgaOv)
contr <- makeContrasts(groupHRD - groupHRP, levels = modelMatrixTcgaOv)
tmp <- contrasts.fit(fitTcgaOv, contrasts = contr)
tmp <- eBayes(tmp, robust = TRUE)
topTable <- topTable(tmp, sort.by = "P", n = Inf)
length(which(topTable$adj.P.Val < 0.05))
# [1] 2604
# write_csv(topTable, "~/identifiHR/identifiHR_featureSelection_deAnalysisTopTable_20240317.csv")
dt <- decideTests(tmp, p.value = 0.05)
summary(dt)
#         groupHRD - groupHRP
# Down                  1315
# NotSig               13297
# Up                    1289

# volcanoPlot -------------------------------------------------------------

topTable %>% 
  dplyr::mutate(diffExp = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Upregulated",
                                    adj.P.Val < 0.05 & logFC < 0 ~ "Downregulated",
                                    TRUE ~ "Not significant")) %>% 
  dplyr::mutate(label = case_when(adj.P.Val < 6.044680e-06 ~ hgnc_symbol,
                                  TRUE ~ NA)) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col = diffExp)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("Downregulated" = "steelblue2","Not significant" =  "grey","Upregulated" =  "indianred1",
                                name = "Differential expression")) +
  theme_minimal() +
  geom_label_repel(aes(label = label),
                   box.padding   = 0.4, 
                   point.padding = 0.2,
                   segment.color = 'grey50') +
  xlab(bquote(log[2]*FC)) +
  ylab(bquote(-log[10]*"(adjusted p value)"))

topTable <- read_csv("~/identifiHR/identifiHR_featureSelection_deAnalysisTopTable_20240317.csv")
deGenesDf <- topTable %>% 
  dplyr::filter(adj.P.Val < 0.05)
deGenes <- deGenesDf$ENSEMBL
deGenesEnt <- deGenesDf$entrezgene_id

 # deResultAnalysis --------------------------------------------------------

limma::plotMA(tmp, coef = 1, status = dt[ ,"groupHRD - groupHRP"], hl.col = c("steelblue2","indianred1"), hl.pch = 20, bg.col = "grey", main = "HRD - HRP")
o <- order(tmp$p.value[,"groupHRD - groupHRP"])
x <- tmp$Amean
y <- tmp$coefficients[,"groupHRD - groupHRP"]
G <- tmp$genes$hgnc_symbol
text(x[o[c(1:5, 8)]], y[o[c(1:5, 8)]], labels = G[o[c(1:5, 8)]]) # to show the top 8 DE gene symbols

# geneOntology ------------------------------------------------------------

goDeSimple <- goana(deGenesEnt) %>% 
  rownames_to_column(var = "GO") %>% 
  mutate(geneRatio = DE / N) %>% 
  dplyr::filter(P.DE < 0.05) %>% 
  arrange(P.DE) %>% 
  dplyr::filter(Ont == c("BP", "MF"))
goDeSimpleSub <- goDeSimple[1:10, ]
ggplot(goDeSimpleSub, aes(x = geneRatio, y = Term, 
                   fill = P.DE)) + 
  geom_col() +
  scale_fill_gradient(low = "indianred1", high = "steelblue1", name = "P value") +
  theme_minimal() + 
  ylab("") + 
  xlab("Gene ratio") + 
  theme() +
  ggtitle("GO enrichment analysis")

goDe <- goana(tmp, geneid = tmp[["genes"]]$entrezgene_id, FDR = 0.05)
goDe <- goDe %>% 
  rownames_to_column(var = "GO") %>% 
  mutate(geneRatioUp = Up / N) %>% 
  mutate(geneRatioDown = Down / N) 
goDeUp <- goDe %>% 
  dplyr::filter(P.Up < 0.05) %>% 
  arrange(P.Up) 
goDeUp <- goDeUp[1:20, ]
goDeDown <- goDe %>% 
  dplyr::filter(P.Down < 0.05) %>% 
  arrange(P.Down)
goDeDown <- goDeDown[1:20, ]
# write_csv(goDe, "~/identifiHR/deFeatureSelection/identifiHR_featureSelection_geneOntology_20240318.csv")

ggplot(goDeUp, aes(x = geneRatioUp, y = Term, 
                        fill = P.Up)) + 
  geom_col() +
  scale_fill_gradient(low = "indianred1", high = "steelblue1", name = "P value") +
  theme_minimal() + 
  ylab("") + 
  xlab("Gene ratio") + 
  theme() +
  ggtitle("GO enrichment analysis")

ggplot(goDeDown, aes(x = geneRatioUp, y = Term, 
                   fill = P.Down)) + 
  geom_col() +
  scale_fill_gradient(low = "indianred1", high = "steelblue1", name = "P value") +
  theme_minimal() + 
  ylab("") + 
  xlab("Gene ratio") + 
  theme() +
  ggtitle("GO enrichment analysis")

# frequencyDeGenes --------------------------------------------------------

sigLocTopTable <- topTable %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  mutate(de = case_when(logFC > 0 ~ "up",
                        TRUE ~ "down")) %>% 
  dplyr::filter(chromosome_name != "NA")
chrOrder <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
              "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X" )

ggplot(sigLocTopTable, aes(x = factor(chromosome_name, levels = chrOrder), fill = de)) +
  geom_bar(colour = "black") +
  theme_minimal() +
  xlab("Chromosome") +
  ylab("Frequency") +
  scale_fill_manual(values = c("down" = "steelblue1",
                               "up"="indianred1")) 

# frequencyByChrArm -------------------------------------------------------

ensTop <- sigLocTopTable$ENSEMBL
useEns <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")
ensAnno <- getBM(attributes = c( "ensembl_gene_id", "gene_biotype", "band"),
                 filters = "ensembl_gene_id",
                 values = ensTop, 
                 mart = useEns) %>% 
  dplyr::rename(ENSEMBL = ensembl_gene_id) %>% 
  left_join(., sigLocTopTable, by = "ENSEMBL") %>% 
  dplyr::mutate(arm = str_sub(band, start = 1, end = 1),
                chrArm = paste0(chromosome_name, sep = ".", arm)) 

chrArmOrder <- c("1.p",  "1.q", "2.p",  "2.q", "3.p",  "3.q",  "4.p",  "4.q",  "5.p",  "5.q",  "6.p",  "6.q",
                 "7.p",  "7.q", "8.p", "8.q",  "9.p",  "9.q", "10.p", "10.q", "11.p", "11.q", "12.p", "12.q",
                 "13.q", "14.q", "15.q", "16.p", "16.q", "17.p", "17.q", "18.p", "18.q", "19.p",
                 "19.q",  "20.p", "20.q", "21.q", "22.q",  "X.p",  "X.q")
ensAnno %>% 
  ggplot(aes(x = factor(chrArm, levels = chrArmOrder), fill = de)) +
           geom_bar(colour = "black") +
           theme_minimal() +
           xlab("Chromosome arm") +
           ylab("Frequency") +
           scale_fill_manual(values = c("down" = "steelblue1",
                                        "up"="indianred1"))

# locationDeGene ----------------------------------------------------------

chrOrder <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
              "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
locationOrder <- topTable %>% 
  dplyr::mutate(chr = case_when(chromosome_name == "X" ~ "23",
                                TRUE ~ chromosome_name)) %>% 
  dplyr::mutate(chr = as.numeric(chr)) %>% 
  mutate(location = paste0(chr,  sep = ".", start_position, sep = "-", ENSEMBL)) %>% 
  mutate(locationAlone = paste0(chr,  sep = ".", start_position)) %>% 
  arrange(chr, start_position) 
locationOrderFull <- locationOrder$location
locationOrderVec <- locationOrder %>% 
  group_by(chr) %>% 
  summarise(minStart = min(start_position)) %>% 
  mutate(locationAlone = paste0(chr, sep = ".", minStart))
locationOrderVec <- left_join(locationOrderVec, locationOrder, by = "locationAlone")
locationOrderVec <- locationOrderVec$location
locationOrder %>%
  mutate(de = case_when(adj.P.Val < 0.01 & logFC > 0 ~ "up",
                        adj.P.Val < 0.01 & logFC < 0 ~ "down")) %>% 
  dplyr::mutate(label = case_when(adj.P.Val < 6.144680e-06 ~ hgnc_symbol,
                                  TRUE ~ NA)) %>% 
  ggplot(aes(x = factor(location, levels = locationOrderFull), y = -log10(adj.P.Val), col = de, alpha = de)) +
  geom_point() +
  theme_minimal()  +
  xlab("Location by chromosome") +
  ylab(bquote(-log[10]*"(adjusted p value)")) +
  scale_color_manual(values = c("down" = "steelblue1",
                                "up"="indianred1")) +
  scale_alpha_manual(values = c("down" = 1, "up" = 1),
                     na.value = 0.4) +
  scale_x_discrete(breaks = locationOrderVec, labels = chrOrder, guide = guide_axis(n.dodge = 2)) +
  geom_label_repel(aes(label = label),
                   box.padding   = 0.4, 
                   point.padding = 0.2,
                   segment.color = 'grey50') 

locationOrder %>%
  dplyr::filter(chromosome_name == "8") %>% 
  mutate(de = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "up",
                        adj.P.Val < 0.05 & logFC < 0 ~ "down")) %>% 
  ggplot(aes(x = start_position, y = -log10(adj.P.Val), col = de, alpha = de)) +
  geom_point() +
  theme_minimal()  +
  xlab("Location on chromosome 8") +
  ylab(bquote(-log[10]*"(adjusted p value)")) +
  scale_color_manual(values = c("down" = "steelblue1",
                                "up"="indianred1")) +
  scale_alpha_manual(values = c("down" = 1, "up" = 1),
                     na.value = 0.4)
locationOrder %>%
  dplyr::filter(chromosome_name == "5") %>% 
  mutate(de = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "up",
                        adj.P.Val < 0.05 & logFC < 0 ~ "down")) %>% 
  ggplot(aes(x = start_position, y = -log10(adj.P.Val), col = de, alpha = de)) +
  geom_point() +
  theme_minimal()  +
  xlab("Location on chromosome 5") +
  ylab(bquote(-log[10]*"(adjusted p value)")) +
  scale_color_manual(values = c("down" = "steelblue1",
                                "up"="indianred1")) +
  scale_alpha_manual(values = c("down" = 1, "up" = 1),
                     na.value = 0.4)

locationOrder %>%
  dplyr::filter(chromosome_name == "19") %>% 
  mutate(de = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "up",
                        adj.P.Val < 0.05 & logFC < 0 ~ "down")) %>% 
  ggplot(aes(x = start_position, y = -log10(adj.P.Val), col = de, alpha = de)) +
  geom_point() +
  theme_minimal()  +
  xlab("Location on chromosome 19") +
  ylab(bquote(-log[10]*"(adjusted p value)")) +
  scale_color_manual(values = c("down" = "steelblue1",
                                "up"="indianred1")) +
  scale_alpha_manual(values = c("down" = 1, "up" = 1),
                     na.value = 0.4)

locationOrder %>% 
  mutate(de = case_when(adj.P.Val < 0.01 & logFC > 0 ~ "up",
                        adj.P.Val < 0.01 & logFC < 0 ~ "down")) %>% 
  ggplot(aes(x = factor(location, levels = locationOrderFull), y = logFC, col = de, alpha = de)) +
  geom_point() +
  theme_minimal()  +
  xlab("Location by chromosome") +
  ylab("logFC") +
  scale_alpha_manual(values = c("down" = 1, "up" = 1),
                     na.value = 0.4) +
  scale_color_manual(values = c("down" = "steelblue1",
                                "up"="indianred1")) +
  scale_x_discrete(breaks = locationOrderVec, labels = chrOrder) 


# subsetCountsToDeGenes ---------------------------------------------------

# also transform to log2counts-per-million
trainData <- tcgaCountsFil[ ,trainIndices]
dim(trainData)
# [1] 52125   288
trainCases <- tcgaClinClean[colnames(trainData), ]
table(colnames(trainData) == rownames(trainCases))
# TRUE 
# 288
trainLabels <- trainCases$hrStatus
trainDataDe <- trainData[deGenes, ]
trainDataCpm <- cpm(trainDataDe, log = TRUE)

# zScore ------------------------------------------------------------------

# calculate mean abundance per gene and store as a vector

meanGene <- rowMeans(trainDataCpm)
meanGeneDf <- as.data.frame(meanGene)
meanGeneIdx <- rownames(meanGeneDf)
  
# calculate standard deviation in abundance per gene and store as a vector

sdGene <- apply(trainDataCpm, 1, sd)
sdGeneDf <- as.data.frame(sdGene)
sdGeneIdx <- rownames(sdGeneDf)

table(rownames(meanGeneDf) == rownames(sdGeneDf))
# TRUE 
# 2604 

# apply zscore

trainDataCpmZ <- (trainDataCpm - meanGene) / sdGene
trainDataCpmZ <- t(trainDataCpmZ)
zDist <- trainDataCpmZ %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = "Chr",
               values_to = "z")

hist(zDist$z)
summary(zDist$z)
trainPca <- PCA(trainDataCpmZ)
fviz_eig(trainPca, 
         addlabels = TRUE, 
         ylim = c(0,100))

fviz_pca_ind(trainPca,
             axes = c(1, 2),
             geom = "point",
             repel = TRUE,
             col.ind = trainLabels,
             legend.title = "HR status"
) +
  labs("") +
  theme_minimal()

# scDeResultTest ----------------------------------------------------------

overlapModelTopEns <- overlapModelTop %>% 
  + dplyr::filter(sig == "sig")
overlapModelTopEns <- overlapModelTopEns$ENSEMBL

scOverlap <- trainDataCpmZ[ ,colnames(trainDataCpmZ) %in% overlapModelTopEns]

trainPca <- PCA(scOverlap)
fviz_eig(trainPca, 
         addlabels = TRUE, 
         ylim = c(0,100))

fviz_pca_ind(trainPca,
             axes = c(1, 2),
             geom = "point",
             repel = TRUE,
             col.ind = trainLabels,
             legend.title = "HR status"
) +
  labs("") +
  theme_minimal()

# hyperparameterTuning ----------------------------------------------------

# Tune parameters with glmnet only

set.seed(NULL)
set.seed(220310)
set.seed(220310)

alphasOfInterest <- seq(0, 1, 0.01)

# do all cross validations for each alpha
cvs <- lapply(alphasOfInterest, function(curAlpha) {
  cv.glmnet(x = trainDataCpmZ, 
            y = trainLabels, 
            alpha = curAlpha, 
            nfolds = 5,
            keep = TRUE,
            type.measure = "auc",
            family = "binomial")
})

# selectOptimumAlpha ------------------------------------------------------

set.seed(NULL)
set.seed(220310)
set.seed(220310)
# collect the measure (auc) and optimum lambda.min for each alpha
optimumPerAlpha <- sapply(seq_along(alphasOfInterest), function(curi) {
  curcvs <- cvs[[curi]]
  curAlpha <- alphasOfInterest[curi]
  max.auc <- cvs[[curi]]$cvm[cvs[[curi]]$lambda == cvs[[curi]]$lambda.min]
  sd.auc <- cvs[[curi]]$cvsd[cvs[[curi]]$lambda == cvs[[curi]]$lambda.min]
  c(lambda = curcvs$lambda.min, alpha = curAlpha, cvm = max.auc, cvsd = sd.auc)
}) %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()

ggplot(as.data.frame(optimumPerAlpha), aes(x = alpha, y = cvm)) +
  geom_point() +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd)) +
  theme_minimal() +
  xlab("Alpha") +
  ylab("Cross validation measure: AUC")

# filter to values that are within 1 standard deviation of the optimal AUC, then select the highest alpha and its corresponding lambda, to ensure sparsity.
maxAuc <- max(optimumPerAlpha$cvm)
threshold <- optimumPerAlpha %>% 
  dplyr::filter(cvm == maxAuc) %>% 
  dplyr::mutate(lowerAuc = cvm - cvsd) %>% 
  dplyr::select(lowerAuc) 
threshold <- threshold$lowerAuc
optimumPerAlphaSelect <- optimumPerAlpha %>% 
  dplyr::filter(cvm >= threshold)
bestAlpha <- max(optimumPerAlphaSelect$alpha)

modelParameters <- optimumPerAlphaSelect %>% 
  dplyr::filter(alpha == max(alpha))
bestAlpha <- modelParameters$alpha
print(bestAlpha)
# [1] 0.53
bestLambda <- modelParameters$lambda
print(bestLambda)
# [1] 0.004211567
print(modelParameters$cvm)
# [1] 0.9266607
print(modelParameters$cvsd)
# [1] 0.01390478

set.seed(NULL)
set.seed(220310)
set.seed(220310)
optimalValue <- cv.glmnet(x = trainDataCpmZ, 
                          y = trainLabels, 
                          alpha = 0.53, 
                          nfolds = 5,
                          keep = TRUE,
                          type.measure = "auc",
                          family = "binomial")

# to get mean AUC across all folds
mean(optimalValue[["cvm"]])
# [1] 0.8557531

set.seed(NULL)
set.seed(220310)
set.seed(220310)
optimalAccValue <- cv.glmnet(x = trainDataCpmZ, 
                          y = trainLabels, 
                          alpha = 0.53, 
                          nfolds = 5,
                          keep = TRUE,
                          type.measure = "class",
                          family = "binomial")

# to get mean misclassification error rate across all folds
misErr <- mean(optimalAccValue[["cvm"]])
# [1] 0.2182639

1 - misErr
# [1] 0.7817361

# extractCoefficients -----------------------------------------------------

# lmHrSig <- glmnet( x = trainDataCpmZ,
#                    y = trainLabels,
#                    family = "binomial",
#                    alpha = bestAlpha, # lasso (1), ridge (0), elastic net (1-E)
#                    lamda = lmHrSig$lambda[99],
#                    keep = TRUE)
lmHrSig <- readRDS("~/identifiHR/identifiHR_model209Genes.RDS")

# extract coefficients of lm 
betaCoefHr <- coef(lmHrSig, s = lmHrSig$lambda[99]) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL") %>% 
  dplyr::filter(s1 != 0) %>% 
  dplyr::filter(ENSEMBL != "(Intercept)")
dim(betaCoefHr)
# [1] 209   2

# annotate with gene symbols
betaCoefHrAnno <- left_join(betaCoefHr, ensTcgaAnno, by = "ENSEMBL") %>% 
  dplyr::filter(!duplicated(ENSEMBL))
coefHrVec <- betaCoefHr$ENSEMBL

# examineModelGenes -------------------------------------------------------

plot(lmHrSig)

chrOrder <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
              "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
betaCoefHrAnno <- betaCoefHrAnno %>% 
  dplyr::rename(betaCoef = s1) %>% 
  dplyr::mutate(colour = ifelse(betaCoef < 0, "blue","red"))

ggplot(betaCoefHrAnno, aes(x = start_position, y = (betaCoef*-1), color = as.factor(colour))) +
  geom_segment(aes(x = start_position, xend = start_position, y = 0, yend = (betaCoef*-1)), color = "grey") +
  geom_point(size = 1) +
  theme_light() +
  theme_minimal() +
  xlab("Location") +
  ylab("Beta coefficient") +
  facet_wrap(~ factor(chromosome_name, 
                      levels = chrOrder), 
             scales = "free_x",
             ncol = 3) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("blue" = "indianred1",
                                "red"="steelblue1"),
                     guide = "none")

betaCoefHrAnnoEnt <- betaCoefHrAnno$entrezgene_id
goModel <- goana(betaCoefHrAnnoEnt) %>% 
  dplyr::filter(Ont == c("BP", "MF")) %>% 
  mutate(geneRatio = DE / N) %>% 
  dplyr::filter(P.DE < 0.05) %>% 
  arrange(P.DE) 
goModel <- goModel[1:20, ]
ggplot(goModel, aes(x = geneRatio, y = Term, 
                   fill = P.DE)) + 
  geom_col() +
  scale_fill_gradient(low = "indianred1", high = "steelblue1", name = "P value") +
  theme_minimal() + 
  ylab("") + 
  xlab("Gene ratio") + 
  theme() +
  ggtitle("GO enrichment analysis")

# frequencyModelGenes -----------------------------------------------------

sigLocBetaGene <- betaCoefHrAnno %>% 
  mutate(betaCoefficient = case_when(betaCoefHrAnno$betaCoef > 0 ~ "up",
                        TRUE ~ "down")) %>% 
  dplyr::filter(chromosome_name != "NA")
chrOrder <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
              "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X" )

ggplot(sigLocBetaGene, aes(x = factor(chromosome_name, levels = chrOrder), fill = betaCoefficient)) +
  geom_bar(colour = "black") +
  theme_minimal() +
  xlab("Chromosome") +
  ylab("Frequency") +
  scale_fill_manual(values = c("down" = "steelblue1",
                               "up"="indianred1")) 


# # retrainWithRidgeRegression ----------------------------------------------
# 
# trainDataRidge <- trainDataCpmZ[ ,colnames(trainDataCpmZ) %in%  modelEns]
# 
# # Tune parameters with glmnet only
# 
# set.seed(NULL)
# set.seed(220310)
# set.seed(220310)
# 
# # do all cross validations for each alpha
# lmHrSig <-  cv.glmnet(x = trainDataRidge, 
#             y = trainLabels, 
#             alpha = 0, 
#             nfolds = 5,
#             keep = TRUE,
#             type.measure = "auc",
#             family = "binomial")
# plot(lmHrSig)
# max(lmHrSig[["cvm"]])
# 
# # best lambda --> assoicated with highest cvm aka AUC
# # lmHrSig[["lambda"]][99]
# bestLambda <- 0.26267981 

# testTcgaHoldOut ---------------------------------------------------------

testCases <- tcgaClinClean[colnames(testData), ]
table(rownames(testCases) == colnames(testData))
# TRUE 
# 73 
testLabels <- testCases$hrStatus

# log2cpm transform test data: TCGA
testDataDe <- testData[deGenes, ]
testDataDeIdx <- testDataDe[meanGeneIdx, ]
testDataCpm <- cpm(testDataDeIdx, log = TRUE)

# apply zscore
testDataCpmZ <- (testDataCpm - meanGene) / sdGene
testDataCpmZ <- t(testDataCpmZ)

# subset to only model genes
# testDataCpmZ <- testDataCpmZ[ ,colnames(testDataCpmZ) %in%  modelEns]


# testing trained model in the validation hold out subset from TCGA
fitAssess <- assess.glmnet(lmHrSig,
                           newx = testDataCpmZ,
                           newy = testLabels,
                           s =  bestLambda)
fitAssess$auc
# [1] 0.8574661

plot(lmHrSig[["lambda"]], lmHrSig[["nzero"]])

confMatTcga <- confusion.glmnet(lmHrSig,
                            newx = testDataCpmZ,
                            newy = testLabels,
                            s =  bestLambda)

confMatTcga
# True
# Predicted HRD HRP Total
    # HRD    29   6    35
    # HRP     5  33    38
    # Total  34  39    73
# 
# Percent Correct:  0.8493

# accuracy:
accTcgaOvTest <- (33 + 29) / (33 + 6 + 29 + 5)
accTcgaOvTest
# [1] 0.8493151

# precision:
precisionTcgaOvTest <- (33) / (33 + 5)
precisionTcgaOvTest
# [1] 0.8684211

# recall:
recallTcgaOvTest <- (33) / (33 + 6)
recallTcgaOvTest
# [1] 0.8461538

# misclassification error:
missClassTcgaOvTest <- (6 + 5) / (33 + 6 + 29 + 5)
missClassTcgaOvTest
# [1] 0.1506849

predictionHrTcga <- stats::predict(lmHrSig,
                                   newx = testDataCpmZ,
                                   s = bestLambda,
                                   type = "class") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Sample") 

predictedProbTcga <- stats::predict(lmHrSig,
                                    newx = testDataCpmZ,
                                    s = bestLambda,
                                    type = "response") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Sample") 

predictionHrTcgaDf <- left_join(predictionHrTcga, predictedProbTcga, by = "Sample")
tcgaClinCleanJoin <- tcgaClinClean %>% 
  rownames_to_column(var = "Sample")

predictionHrTcgaDfAcc <- left_join(predictionHrTcgaDf, tcgaClinCleanJoin, by = "Sample") %>% 
  mutate(correctLabel = case_when(s1.x == hrStatus ~ "correct", 
                                  TRUE ~ "incorrect")) %>% 
  dplyr::rename(IdentifiHR = s1.x) %>% 
  dplyr::rename(probability = s1.y)
# write_csv(predictionHrTcgaDfAcc, "/home/users/allstaff/weir.a/identifiHR/incorrectTcgaOvCases/tcgaOvTest_identifiHRPrediictions_20240703.csv")

colCorr <- as.factor(predictionHrTcgaDfAcc$correctLabel)
ggplot(predictionHrTcgaDfAcc, aes(x = HRD, y = probability)) +
  geom_point(aes(colour = colCorr)) +
  geom_vline(xintercept = 42, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_minimal() +
  ylab("Probability") +
  xlab("HRD score") +
  scale_colour_manual(values = c("correct" = "yellowgreen",
                                 "incorrect"= "darkorange"),
                      name = "Model prediction") +
  stat_cor(method = "pearson",
           label.x.npc = 0.6)

# testAocsCohort ----------------------------------------------------------

# log2cpm transform test data: AOCS
testDataCpmAocs <- countsAocsEns[rownames(countsAocsEns) %in% deGenes, ]
testDataCpmAocsIdx <- testDataCpmAocs[meanGeneIdx, ]
testDataCpmAocs <- cpm(testDataCpmAocsIdx, log = TRUE)
table(rownames(testDataCpmAocs) == rownames(meanGeneDf))
# TRUE 
# 2604 
# apply zscore
testDataCpmZAocs <- (testDataCpmAocs - meanGene) / sdGene
testDataCpmZAocs <- t(testDataCpmZAocs)

table(rownames(testDataCpmZAocs) == clinHrAocs$Sample)
aocsLabels <- clinHrAocs$hrStatus
pcaAllNormCounts <- PCA(testDataCpmZAocs[ ,colnames(testDataCpmZAocs) %in% coefHrVec])
fviz_eig(pcaAllNormCounts, 
         addlabels = TRUE, 
         ylim = c(0,100))

fviz_pca_ind(pcaAllNormCounts,
             axes = c(1, 2),
             geom = "point",
             repel = TRUE,
             col.ind = aocsLabels,
             legend.title = "HR status"
) +
  labs("") +
  theme_minimal()

#testDataCpmZAocs <- testDataCpmZAocs[ ,colnames(testDataCpmZAocs) %in%  modelEns]

# testing trained model in the validation hold out subset from TCGA
fitAssess <- assess.glmnet(lmHrSig,
                           newx = testDataCpmZAocs,
                           newy = aocsLabels,
                           s =  bestLambda,
                           family = "binomial")
fitAssess$auc
# [1] 0.9472177

confMatAocs <- confusion.glmnet(lmHrSig,
                            newx = testDataCpmZAocs,
                            newy = aocsLabels,
                            s =  bestLambda,
                            family = "binomial")
confMatAocs
# True
# Predicted HRD HRP Total
    # HRD    49   6    55
    # HRP     3  41    44
    # Total  52  47    99
# 
# Percent Correct:  0.9091 

# accuracy:
accAocsTest <- (41 + 49) / (41 + 6 + 49 + 3)
accAocsTest
# [1] 0.9090909

# precision:
precisionAocsTest <- (41) / (41 + 3)
precisionAocsTest
# [1] 0.9318182

# recall:
recallAocsTest <- (41) / (41 + 6)
recallAocsTest
# [1] 0.8723404

# misclassification error:
missClassAocsTest <- (6 + 3) / (41 + 6 + 49 + 3)
missClassAocsTest
# [1] 0.09090909

predictionHr <- stats::predict(lmHrSig,
                               newx = testDataCpmZAocs,
                               s = bestLambda,
                               type = "class") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Sample") 

predictedProb <- stats::predict(lmHrSig,
                                newx = testDataCpmZAocs,
                                s = bestLambda,
                                type = "response") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Sample") 

predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")

predictionHrDfAcc <- left_join(predictionHrDf, clinHrAocs, by = "Sample") %>% 
  mutate(correctLabel = case_when(s1.x == hrStatus ~ "correct", 
                                  TRUE ~ "incorrect")) %>% 
  dplyr::rename(IdentifiHR = s1.x)

table(predictionHrDfAcc$correctLabel)
# correct incorrect 
# 90         9 

table(predictionHrDfAcc$tissueSource, predictionHrDfAcc$correctLabel)
# correct incorrect
# Ascites                  12         0
# AutopsyTumour             6         0
# FT - normal control       7         0
# PrimaryTumour            65         9

table(predictionHrDfAcc$IdentifiHR, predictionHrDfAcc$tissueSource)
#       Ascites AutopsyTumour FT - normal control PrimaryTumour
# HRD      10             6                   0            39
# HRP       2             0                   7            35

table(predictionHrDfAcc$hrStatus, predictionHrDfAcc$tissueSource)
#       Ascites AutopsyTumour FT - normal control PrimaryTumour
# HRD      10             6                   0            36
# HRP       2             0                   7            38

incorrectAocs <- predictionHrDfAcc %>% 
  dplyr::filter(correctLabel == "incorrect")

correctAocs <- predictionHrDfAcc %>% 
  dplyr::filter(correctLabel == "correct") %>% 
  dplyr::filter(!is.na(Tumor_cell_puritya))

kruskal.test(Tumor_cell_puritya ~ correctLabel, data = predictionHrDfAcc)

# Kruskal-Wallis rank sum test
# 
# data:  Tumor_cell_puritya by correctLabel
# Kruskal-Wallis chi-squared = 1.0923, df = 1, p-value = 0.296

# makeMetricsDf -----------------------------------------------------------

tcgaTestMetrics <- c(accTcgaOvTest, precisionTcgaOvTest, recallTcgaOvTest, missClassTcgaOvTest)
aocsTestMetrics <- c(accAocsTest, precisionAocsTest, recallAocsTest, missClassAocsTest)
metrics <- c("accuracy", "precision", "recall", "misclassificationError")
identifiHRMetrics <- data.frame(metrics, tcgaTestMetrics, aocsTestMetrics) 

# pcaAllGenesTraining -----------------------------------------------------

zScore <- function(x, na.rm = TRUE) {
  return((x - mean(x)) / (sd(x) + 0.00001))
}

trainDataAllGene <- tcgaCountsFil[ ,trainIndices]
trainDataAllGeneCpmZ <- trainDataAllGene %>% 
  cpm(., log = TRUE) %>% 
  apply(., 1, zScore)

table(rownames(trainDataAllGeneCpmZ) == rownames(trainCases))
# TRUE 
# 288

allGenePca <- PCA(trainDataAllGeneCpmZ)
fviz_pca_ind(allGenePca,
             axes = c(1, 2),
             geom = "point",
             col.ind = trainCases$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(allGenePca,
             axes = c(1, 2),
             geom = "point",
             col.ind = trainCases$hrStatus,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

# pcaDeGene ---------------------------------------------------------------

trainDataAllGene <- tcgaCountsFil[ ,trainIndices]
trainDataDeGeneCpmZ <- trainDataAllGene[rownames(trainDataAllGene) %in% deGenes$ENSEMBL, ] %>% 
  cpm(., log = TRUE) %>% 
  apply(., 1, zScore)
table(rownames(trainDataDeGeneCpmZ) == rownames(trainCases))
# TRUE 
# 288

deGenePca <- PCA(trainDataDeGeneCpmZ)
fviz_pca_ind(deGenePca,
             axes = c(1, 2),
             geom = "point",
             col.ind = trainCases$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(deGenePca,
             axes = c(1, 2),
             geom = "point",
             col.ind = trainCases$hrStatus,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

# pcaModelGenes -----------------------------------------------------------

trainDataAllGene <- tcgaCountsFil[ ,trainIndices]
trainDataModelGeneCpmZ <- trainDataAllGene[rownames(trainDataAllGene) %in% betaCoefHrAnno$ENSEMBL, ] %>% 
  cpm(., log = TRUE) %>% 
  apply(., 1, zScore)
table(rownames(trainDataModelGeneCpmZ) == rownames(trainCases))
# TRUE 
# 288

modelGenePca <- PCA(trainDataModelGeneCpmZ)

fviz_pca_ind(modelGenePca,
             axes = c(1, 2),
             geom = "point",
             col.ind = trainCases$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(modelGenePca,
             axes = c(1, 2),
             geom = "point",
             col.ind = trainCases$hrStatus,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

# testCohortPcas ----------------------------------------------------------


# tcgaTestPca -------------------------------------------------------------

table(rownames(testDataCpmZ) == predictionHrTcgaDfAcc$Sample)
# TRUE 
# 73
testDataCpmZModelGene <- testDataCpmZ[ ,colnames(testDataCpmZ) %in% betaCoefHrAnno$ENSEMBL]
dim(testDataCpmZModelGene)
# [1]  73 209  
tcgaTestPca <- PCA(testDataCpmZModelGene)
fviz_pca_ind(tcgaTestPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = predictionHrTcgaDfAcc$hrStatus,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

fviz_pca_ind(tcgaTestPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = predictionHrTcgaDfAcc$IdentifiHR,
             pointshape = 19,
             palette = c("HRP" = "yellowgreen", "HRD" = "olivedrab"),
             legend.title = "Model status"
) +
  theme_minimal() 

# correlation between Dim 1 and Dim 2 with HRD score
tcgaTestPcaInd <- get_pca_ind(tcgaTestPca)
tcgaTestCoord <- tcgaTestPcaInd[["coord"]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% 
  left_join(., predictionHrTcgaDfAcc, by = "Sample")

ggplot(tcgaTestCoord, aes(x = HRD, y = Dim.1)) +
  geom_point(aes(colour = colCorr)) +
  geom_vline(xintercept = 42, linetype = "dashed") +
  theme_minimal() +
  ylab("PCA Dim1 coordinate") +
  xlab("HRD score") +
  scale_colour_manual(values = c("correct" = "yellowgreen",
                                 "incorrect"= "darkorange"),
                      name = "Model prediction") +
  stat_cor(method = "pearson",
           label.x.npc = 0.6)

# aocsTestPca -------------------------------------------------------------

table(rownames(testDataCpmZAocs) == predictionHrDfAcc$Sample)
# TRUE 
# 99
testDataCpmZAocsModelGene <- testDataCpmZAocs[ ,colnames(testDataCpmZAocs) %in% betaCoefHrAnno$ENSEMBL]
dim(testDataCpmZAocsModelGene)
# [1]  99 209  
aocsTestPca <- PCA(testDataCpmZAocsModelGene)

fviz_pca_ind(aocsTestPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = predictionHrDfAcc$hrStatus,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

fviz_pca_ind(aocsTestPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = predictionHrDfAcc$IdentifiHR,
             pointshape = 19,
             palette = c("HRP" = "yellowgreen", "HRD" = "olivedrab"),
             legend.title = "Model status"
) +
  theme_minimal() 

# probabilityDistributionIdentifiHR ---------------------------------------

predictionHrTcgaDfAccSub <- predictionHrTcgaDfAcc %>% 
  dplyr::select(Sample, IdentifiHR, s1.y, hrStatus, correctLabel)

predictionHrDfAccSub <- predictionHrDfAcc %>% 
  dplyr::select(Sample, IdentifiHR, s1.y, hrStatus)

predictTests <- rbind(predictionHrTcgaDf, predictionHrDfAccSub) %>% 
  dplyr::mutate(cohort = str_sub(Sample, start = 1, end = 4))

ggplot(predictTests, aes(x = s1.y, y = cohort, fill = cohort)) +
  geom_violin() +
  geom_point(position = "jitter") +
  ylab("Testing cohort") +
  xlab("IdentifiHR probability") +
  theme_minimal() +
  scale_fill_manual(values = c("burlywood1", "darkseagreen2"))+
  theme(legend.position = "none")

# # # checkMintie -------------------------------------------------------------
# # Nadia Davidson checking mintie results
# # # loadMintie
# mintie <- read_csv("~/mintie_TCGA_OV_boundOutputs_20230109.csv")
# fileName <- read.delim("~/TCGA/gdc_sample_sheet.2023-01-11.tsv")
# clinFreqMintie <- read_csv("~/allClinCasesFracHrd_TCGAOv_20230421.csv")
# 
# # joinMintieWithClin
# colnames(fileName)[1] <- "sample"
# namedMintie <- left_join(mintie, fileName, by = "sample") %>%
#   mutate(Sample.ID.1 = str_sub(Sample.ID, start = 1, end = 15))
# namedMintieClin <- left_join(namedMintie, clinFreqMintie, by = "Sample.ID.1", relationship = "many-to-many")
# 
# namedMintieClin <- namedMintieClin %>%
#   dplyr::rename(Sample = Sample.ID.1)
# namedMintieClinIncorrect <- left_join(namedMintieClin, predictionHrTcgaDf, by = "Sample")
# namedMintieClinIncorrect <- namedMintieClinIncorrect %>% 
#   mutate(correctLabel = case_when(s1.x == "HRD" & hrStatus == "HRD" ~ "correct",
#                                   s1.x == "HRP" & hrStatus == "HRP" ~ "correct",
#                                   s1.x == "HRP" & hrStatus == "HRD" ~ "incorrect",
#                                   s1.x == "HRD" & hrStatus == "HRP" ~ "incorrect",
#                                   TRUE ~ NA)) %>%
#   dplyr::filter(correctLabel == "incorrect") %>% 
#   dplyr::rename(IdentifiHR_prediction = s1.x) %>% 
#   dplyr::rename(IdentifiHR_probability = s1.y)
# write_csv(namedMintieClinIncorrect, "/vast/scratch/users/weir.a/identifiHR_incorrectTcgaTestCases_mintieResults_20240318.csv")
# rm(mintie)
# rm(namedMintie)
# rm(namedMintieClin)

# saveModel ---------------------------------------------------------------

saveRDS(lmHrSig, "~/identifiHR/identifiHR_model209Genes.RDS")

# saveCohortsClin ---------------------------------------------------------

tcgaClinClean %>%
  rownames_to_column(var = "Sample") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvClin.csv")

trainCases %>% 
  rownames_to_column(var = "Sample") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvClinTrainingCases.csv")

testCases %>%
  rownames_to_column(var = "Sample") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvClinTestingCases.csv")

# saveCohortsRawData ------------------------------------------------------

tcgaCountsFil %>% 
  rownames_to_column(var = "ENSEMBL") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvRawCountsAllCases.csv")

trainData %>% 
  rownames_to_column(var = "ENSEMBL") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvRawCountsTrainingCases.csv")

testData %>%
  rownames_to_column(var = "ENSEMBL") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvRawCountsTestingCases.csv")

# saveCohortsProcessedData ------------------------------------------------

trainDataCpmZ %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvProcessedCountsTrainingCases.csv")

testDataCpmZ %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>%
  write_csv(., "~/identifiHR/cleanCohorts/tcgaOvProcessedCountsTestingCases.csv")

testDataCpmZAocs %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>%
  write_csv(., "~/identifiHR/cleanCohorts/aocsOvProcessedCountsTestingCases.csv")

# saveZScoreVectors -------------------------------------------------------

meanGene %>% 
  as.data.frame() %>% 
  dplyr::rename(mean = 1) %>% 
  rownames_to_column(var = "ENSEMBL") %>% 
  write_csv("~/identifiHR/cleanCohorts/identifiHR_meanModelGenes.csv")

sdGene %>% 
  as.data.frame() %>% 
  dplyr::rename(sd = 1) %>% 
  rownames_to_column(var = "ENSEMBL") %>% 
  write_csv("~/identifiHR/cleanCohorts/identifiHR_sdModelGenes.csv")

# saveAnnotatedModelGenes -------------------------------------------------

betaCoefHrAnno %>% 
  write_csv("~/identifiHR/cleanCohorts/identifiHR_modelGenes209Annotated.csv")

# saveTestingMetrics ------------------------------------------------------

identifiHRMetrics %>% 
  write_csv("~/identifiHR/benchmarking/identifiHR_testMetricsTcgaAocs.csv")

# savePredictions ---------------------------------------------------------

predictionHrTcgaDfAcc %>% 
  write_csv("~/identifiHR/benchmarking/predicitons/identifiHR_predictionsTcgaOv_20240317.csv")

predictionHrDfAcc %>% 
  write_csv("~/identifiHR/benchmarking/predicitons/identifiHR_predictionsAocs_20240317.csv")
