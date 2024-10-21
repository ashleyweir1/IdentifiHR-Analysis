

# installPackages ---------------------------------------------------------

# remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("Glimma")
# BiocManager::install("biomaRt")
# install.packages("FactoMineR")
# install_github("kassambara/factoextra")

# loadPackages ------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(edgeR)
library(limma)
library(Glimma)
library(readxl)
library(ggplot2)
library(biomaRt)
library(FactoMineR)
library(devtools)
library(factoextra)
library(readr)

# loadData ----------------------------------------------------------------

# cells <- readRDS("/vast/projects/lab_davidson/weir.a/Cohort_processed_filtered_annotated_release-2.rds")
# the slot "images" was empty - filled with an empty list
cells <- readRDS("/vast/projects/lab_davidson/weir.a/Ovarian.cancer.super_processed_filtered_annotated_release.rds")
cells@images <- list()
meta <- cells@meta.data

# loadMetadata ------------------------------------------------------------

hrStatusMeta <- read_excel("identifiHR/singleCellTesting/41586_2022_5496_MOESM4_ESM (1).xlsx")

meta <- left_join(meta, hrStatusMeta, by = "patient_id")
cells@meta.data <- left_join(cells@meta.data,hrStatusMeta, by = "patient_id")

# plotMeta ----------------------------------------------------------------

ggplot(meta, aes(x = tumor_type)) +
  geom_histogram(stat = "count") +
  theme_minimal()
# try various colnames to look at frequency of cells per site/type etc

ggplot(meta, aes(x = consensus_hr_status)) +
  geom_histogram(stat = "count") +
  theme_minimal()

# pseudobulking -----------------------------------------------------------

Idents(cells) <- cells@meta.data$patient_id
pbAll <- AggregateExpression(object = cells, slot = "counts")
pbAll <- as.data.frame(pbAll[["RNA"]])
samp <- colnames(pbAll)

rownames(hrStatusMeta) <- hrStatusMeta$patient_id
hrStatusMetaSamp <- hrStatusMeta[rownames(hrStatusMeta) %in% samp, ]
rownames(hrStatusMetaSamp) <- hrStatusMetaSamp$patient_id
hrStatusMetaSamp <- hrStatusMetaSamp[samp, ]
table(colnames(pbAll) == hrStatusMetaSamp$patient_id)

sampSite <- meta %>% 
  dplyr::select(patient_id, tumor_megasite, tumor_supersite, tumor_site, tumor_subsite, tumor_type, cell_type_super) %>% 
  dplyr::filter(!duplicated(patient_id))

primSamp <- sampSite %>% 
  dplyr::filter(tumor_type == "Primary")
primSamp <- primSamp$patient_id  

# createFreqTableOfSampSites ----------------------------------------------

megasite <- as.data.frame(table(meta$patient_id, meta$tumor_megasite)) %>% 
  pivot_wider(names_from = Var2,
              values_from = Freq) %>% 
  dplyr::rename(patient_id = Var1) 
colnames(megasite)[-1] <- paste0(colnames(megasite)[-1],sep = "_", "megasite")

supersite <- as.data.frame(table(meta$patient_id, meta$tumor_supersite)) %>% 
  pivot_wider(names_from = Var2,
              values_from = Freq) %>% 
  dplyr::rename(patient_id = Var1) 
colnames(supersite)[-1] <- paste0(colnames(supersite)[-1],sep = "_", "supersite")

tumortype <- as.data.frame(table(meta$patient_id, meta$tumor_type)) %>% 
  pivot_wider(names_from = Var2,
              values_from = Freq) %>% 
  dplyr::rename(patient_id = Var1) %>% 
  mutate(tumorTypeRatioPtoM = (Primary + 1) / Metastasis)
colnames(tumortype)[-1] <- paste0(colnames(tumortype)[-1],sep = "_", "tumortype")

supercell <- as.data.frame(table(meta$patient_id, meta$cell_type_super)) %>% 
  pivot_wider(names_from = Var2,
              values_from = Freq) %>% 
  dplyr::rename(patient_id = Var1) 
colnames(supercell)[-1] <- paste0(colnames(supercell)[-1],sep = "_", "supercell")

freqSite <- left_join(megasite, supersite, by = "patient_id")
freqSite <- left_join(freqSite, tumortype, by = "patient_id")
freqSite <- left_join(freqSite, supercell, by = "patient_id")
hrStatusMetaSampFreq <- left_join(hrStatusMetaSamp, freqSite, by = "patient_id")

# annotateGenes -----------------------------------------------------------

geneSym <- rownames(pbAll)
useSym <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")
annoSym <- getBM(attributes = c( "ensembl_gene_id","hgnc_symbol","entrezgene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "band"),
                    filters = "hgnc_symbol",
                    values = geneSym, 
                    mart = useSym)

annoSym <- annoSym %>% 
  dplyr::filter(hgnc_symbol %in% geneSym)
rownames(annoSym) <- annoSym$hgnc_symbol

# createDGE ---------------------------------------------------------------

dge <- DGEList(counts = pbAll, samples = hrStatusMetaSampFreq, group = hrStatusMetaSampFreq$consensus_hr_status)

pca <- PCA(dge$counts)
fviz_screeplot(pca, ncp=10)
plot.PCA(pca, axes = c(1,2))

fviz_pca_ind(pca,
             axes = c(1, 2),
             # axes = c(2, 3),
             #axes = c(3, 4),
             geom = "point",
             repel = TRUE,
             col.ind = dge$samples$Ovarian.cancer.super_supercell,
             legend.title = "HR status"
) +
  labs("") +
  theme_minimal()

# DE using primary samples only
# dgePrim <- dge[ ,primSamp]
# dge <- dgePrim


# filterLowlyExpressedGenes -----------------------------------------------

keep <- filterByExpr(dge, group = dge$samples$consensus_hr_status) 
dgeFil <- dge[keep, ,keep.lib.sizes = FALSE ]
dim(dgeFil)
# [1] 5000   41 --> all cells
#  21534    41 --> only HGSC cancer cells
# normalise ---------------------------------------------------------------

dgeFil <- calcNormFactors(dgeFil)
group <-  interaction(dge$samples$group)

# modelMatrix -------------------------------------------------------------

modelMatrix <- model.matrix(~ 0 + group, na.action = 'na.pass')

# voomAndFit --------------------------------------------------------------

voom <- voom(dgeFil, modelMatrix, plot = TRUE)
fit <- lmFit(voom, modelMatrix)

# contrast ----------------------------------------------------------------

contr <- makeContrasts(groupHRD - groupHRP, levels = modelMatrix)
tmp <- contrasts.fit(fit, contrasts = contr)

# eBayes ------------------------------------------------------------------

tmp <- eBayes(tmp)

# topTable ----------------------------------------------------------------

topTable <- topTable(tmp, sort.by = "P", n = Inf) %>% 
  rownames_to_column(var = "hgnc_symbol")
length(which(topTable$adj.P.Val < 0.05))
# [1] 562
sigTop <- topTable %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::select(hgnc_symbol)
sigTop$sig <- "sig"

# compareWithModel --------------------------------------------------------

modelGenes <- read_csv("identifiHR/cleanCohorts/identifiHR_modelGenes209Annotated.csv")
modelTopTable <- read_csv("identifiHR/identifiHR_featureSelection_deAnalysisTopTable_20240317.csv") %>% 
  dplyr::filter(adj.P.Val < 0.05)

overlapModel <- left_join(modelGenes, sigTop, by = "hgnc_symbol")
table(overlapModel$sig)
# sig 
# 8 

overlapModelTop <- left_join(modelTopTable, sigTop, by = "hgnc_symbol")
table(overlapModelTop$sig)
# sig 
# 194 

# analysisToDo ------------------------------------------------------------

# where are the DE genes located? See which regions are TME-dependent?
# include met v primary in model matrix
# DE of other cell types
# number of sig genes in DE and model
# for the 194 genes that are sig in both bulk and scRNAseq DE analysis: 
# --> can they separate scRNAseq cohort on PCA?
# --> can they separate the bulk TCGA cohort (maybe try only training?)
# if yes --> make these the features of the model and do a ridge regression, then test and validate


# extra -------------------------------------------------------------------


##Pseudobulk

#Change idents to *sample* for grouping.
Idents(allcells) <- allcells$Sample
pseudobulk_matrix <- AggregateExpression( allcells,  slot = 'counts', assays='SCT' )[['SCT']]
pseudobulk_matrix[1:35,]
pseudobulk_matrix <- pseudobulk_matrix[,c("A1","B1","C1", "A3","B3","C3")]
replicates_lookup <- c("A1"="HC","B1"="HC","C1"="HC", "A3"="INF","B3"="INF","C3"="INF")
dge <- DGEList(pseudobulk_matrix)
logcounts <- cpm(dge, log = TRUE)
plotMDS(dge)
#interactive MDS plot with Glimma
glMDSPlot(dge, groups=colnames(dge), main="Allcells_HCvINF_DE", gene.selection = c("pairwise"),html = "JiyiPang_Allcells_HCvINF_PairwiseDE_MDS_20221018")

dge <- calcNormFactors(dge)
condition <- as.factor(replicates_lookup)
levels(condition) #checking everything worked
design <- model.matrix(~condition)

vm  <- voom(dge, design = design,plot = TRUE)
plotMDS(vm,top = 1000, pch = 19, cex = 1.5)
glMDSPlot(vm, groups=colnames(vm), main="Allcells_HCvINF_voom_DE", gene.selection = c("pairwise"),html = "JiyiPang_Allcells_HCvINF_voomDE_MDS_20221018")

fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
dim(fit)
#> dim(fit)
#[1] 20755     2
glimmaMA(fit, dge=dge,  main="AllCells_DE_HCvINF")
glimmaVolcano(fit, dge=dge, main="AllCells_DE_HCvINF")

summa.fit <- decideTests(fit)
summary(summa.fit)
de_result_pseudobulk <- topTable(fit, n = Inf, adjust.method = "BH")
#> Removing intercept from test coefficients
de_result_pseudobulk <- arrange(de_result_pseudobulk , adj.P.Val)
de_result_pseudobulk$genes <- rownames(de_result_pseudobulk)
write.csv(de_result_pseudobulk, "./JiyiPang_AllCells_HCvINF_DE_20221012.csv")

