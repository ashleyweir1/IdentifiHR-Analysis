# background --------------------------------------------------------------

# testing the BRCAness transcriptional signature as a classifier of HR status in HGSC.

# loadPackages ------------------------------------------------------------

library(readxl)
library(tidyverse)
library(stats)
library(factoextra)
library(FactoMineR)
library(ComplexHeatmap)
library(biomaRt)
library(edgeR)
library(fossil)

# tcgaData ----------------------------------------------------------------

countsTcgaBrcaness <- read.csv("~/identifiHR/identifiHR_finalCountsTcgaOv.csv")

# TCGA clinical
clinTcgaBrcaness <- read.csv("~/identifiHR/identifiHR_tcgaClinicalComplete_20230304.csv") %>% 
  dplyr::rename(Sample = Sample_Id)

# TCGA test case identifiers
tcgaTest <- read_csv("~/identifiHR/IdentifiHR/data/cleanCohorts/tcgaOvClinTestingCases.csv")

# subsetToOnlyTestCases ---------------------------------------------------

colnames(countsTcgaBrcaness) <- str_replace_all(colnames(countsTcgaBrcaness),  "\\.", "-")
rownames(countsTcgaBrcaness) <- countsTcgaBrcaness$ENSEMBL
countsTcgaBrcaness <- countsTcgaBrcaness[ ,colnames(countsTcgaBrcaness) %in% tcgaTest$Sample]
countsTcgaBrcaness <- countsTcgaBrcaness[ ,tcgaTest$Sample]
clinTcgaBrcaness <- clinTcgaBrcaness[clinTcgaBrcaness$Sample %in% tcgaTest$Sample, ] 
rownames(clinTcgaBrcaness) <- clinTcgaBrcaness$Sample
clinTcgaBrcaness <- clinTcgaBrcaness[tcgaTest$Sample, ]

# test ordering of clin samples is the same
table(colnames(countsTcgaBrcaness) == clinTcgaBrcaness$Sample)
# TRUE 
# 73 

# loadAocs ----------------------------------------------------------------

countsAocsBrcaness <- read.csv("~/identifiHR/identifiHR_finalCountsAocs.csv")

# AOCS clinical
clinAocsBrcaness <- read.csv("~/identifiHR/clincialAocsComplete_20240312.csv")

# loadBrcanessGenes -------------------------------------------------------

brcaness <- read_excel("~/identifiHR/benchmarking/brcanessGeneSig/brcanessSig_36497135.xlsx", col_names = FALSE)

# annotateWithEns ---------------------------------------------------------

# annotate BRCAness signatture with ensembl IDs

brcanessSymDf <- brcaness %>% 
  dplyr::rename(hgnc_symbol = 1)
brcanessSym <- brcanessSymDf$hgnc_symbol
useEns <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")
ensBrcaness <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id","entrezgene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"),
                          filters = "hgnc_symbol",
                          values = brcanessSym, 
                          mart = useEns) %>% 
  dplyr::filter(!duplicated(hgnc_symbol))
colnames(ensBrcaness)[2] <- "ENSEMBL"

missingGene <- left_join(brcanessSymDf, ensBrcaness, by = "hgnc_symbol")
# hgnc_symbol
# MRE11A  
missingGene[15, 2] <- "ENSG00000020922"
missingGene[missingGene$hgnc_symbol == "MRE11A",]  # MRE11A was not annotated, so ENSEMBL ID was manually added
# hgnc_symbol ensembl_gene_id entrezgene_id chromosome_name start_position end_position strand gene_biotype
# <chr>       <chr>                   <int> <chr>                    <int>        <int>  <int> <chr>       
#   1 MRE11A      ENSG00000020922            NA NA                          NA           NA     NA NA  

# update PTEN ensembl id --> took id assoicated with patch?
missingGene[1, 2] <- "ENSG00000171862"
missingGene[missingGene$hgnc_symbol == "PTEN",] 
# hgnc_symbol ENSEMBL         entrezgene_id chromosome_name start_position end_position strand gene_biotype  
# <chr>       <chr>                   <int> <chr>                    <int>        <int>  <int> <chr>         
#   1 PTEN        ENSG00000171862          5728 HG2334_PATCH             78462       183684      1 protein_coding
ensBrcanessVec <- missingGene$ENSEMBL

# joinCounts --------------------------------------------------------------

countsTcgaBrcaness <- countsTcgaBrcaness %>% 
  rownames_to_column(var = "ENSEMBL")
countsBrcaness <- left_join(countsTcgaBrcaness, countsAocsBrcaness, by = "ENSEMBL") %>% 
  column_to_rownames(var = "ENSEMBL")
dim(countsBrcaness)
# [1] 52125   1z72

# sampleIndexing ----------------------------------------------------------

clinTcgaBrcaness <- clinTcgaBrcaness %>% 
  dplyr::select(Sample, hrStatus) %>% 
  mutate(cohort = "TCGA")

clinAocsBrcaness <- clinAocsBrcaness %>% 
  dplyr::select(Sample, hrStatus) %>% 
  mutate(cohort = "AOCS")
  
clinBrcaness <- rbind(clinTcgaBrcaness, clinAocsBrcaness)
dim(clinBrcaness)
# [1] 172   3

clinBrcaness <- clinBrcaness[clinBrcaness$Sample %in% colnames(countsBrcaness), ] %>% 
  dplyr::filter(!duplicated(Sample)) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Sample")
clinBrcaness <- clinBrcaness[ colnames(countsBrcaness),]

table(rownames(clinBrcaness) == colnames(countsBrcaness))
# TRUE 
# 172

# logCpmNormalise ---------------------------------------------------------

# get the normalized counts
cpmCountsBrcaness <- cpm(countsBrcaness, log=TRUE)

# subsetTcgaCounts --------------------------------------------------------

countsBrcanessSub <- cpmCountsBrcaness[ensBrcanessVec, ]  
dim(countsBrcanessSub)
# [1]  40 172

# optionalSubsetToTcga ----------------------------------------------------

countsBrcanessSub <- countsBrcanessSub[,grepl("TCGA", colnames(countsBrcanessSub))]
clinBrcaness <- clinBrcaness[colnames(countsBrcanessSub),]

table(rownames(clinBrcaness) == colnames(countsBrcanessSub))
# TRUE 
# 73

# zScore ------------------------------------------------------------------

zScore <- function(x, na.rm = TRUE) {
  return((x - mean(x)) / (sd(x)))
}

cpmZCountsBrcaness <- apply(countsBrcanessSub, 1, zScore) 

# pca ---------------------------------------------------------------------

pcaBrcaness <- PCA(cpmZCountsBrcaness)

fviz_eig(pcaBrcaness, 
         addlabels = TRUE, 
         ylim = c(0,100))

fviz_pca_ind(pcaBrcaness,
             axes = c(1, 2),
             # axes = c(2, 3),
             # axes = c(3, 4),
             geom = "point",
             col.ind = clinBrcaness$hrStatus,
             # col.ind = clinBrcaness$cohort,
             # palette = c("mediumpurple3", "olivedrab"),
             pointshape = 20,
             legend.title = "HR status",
             #legend.title = "Cohort"
) +
  theme_minimal()

# kMeansClustering --------------------------------------------------------

set.seed(NULL)
set.seed(220310)
kMeans <- kmeans(cpmZCountsBrcaness, centers = 2, iter.max = 10, nstart = 1)

kMeansClusVec <- as.data.frame(kMeans$cluster) %>% 
  rownames_to_column(var = "Sample")
clinBrcaness$Sample = rownames(clinBrcaness)

brcanessClinClus <- left_join(clinBrcaness, kMeansClusVec, by = "Sample") %>% 
  mutate(kMeansCluster = as.factor(`kMeans$cluster`))

fviz_pca_ind(pcaBrcaness,
             axes = c(1, 2),
             geom = "point",
             col.ind = brcanessClinClus$kMeansCluster,
             palette = c("lightsalmon1", "lightskyblue"), 
             legend.title = "K means cluster",
             pointshape = 20
) +
  theme_minimal()

# geneContributionsToEachCluster ------------------------------------------

brcanessSigContrib <- kMeans[["centers"]] %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL") %>% 
  left_join(., missingGene, by = "ENSEMBL") %>% 
  dplyr::rename(kClus1 = 2) %>% 
  dplyr::rename(kClus2 = 3)

# makePcaWithAllGenes -----------------------------------------------------

dim(countsBrcaness)
#[1] 52125   172

# get the normalized counts
cpmCountsBrcanessAll <- cpm(countsBrcaness, log=TRUE)

# z score
cpmZCountsBrcanessAll <- apply(cpmCountsBrcanessAll, 1, zScore) 

# pca
pcaBrcanessAll <- PCA(cpmZCountsBrcanessAll)

fviz_eig(pcaBrcanessAll, 
         addlabels = TRUE, 
         ylim = c(0,100))

fviz_pca_ind(pcaBrcanessAll,
             axes = c(1, 2),
             # axes = c(2, 3),
             # axes = c(3, 4),
             geom = "point",
             col.ind = clinBrcaness$hrStatus,
             # col.ind = clinBrcaness$cohort,
             # palette = c("mediumpurple3", "olivedrab"),
             pointshape = 20,
             legend.title = "HR status",
             # legend.title = "Cohort"
) +
  theme_minimal()

fviz_pca_ind(pcaBrcanessAll,
             axes = c(1, 2),
             geom = "point",
             col.ind = brcanessClinClus$kMeansCluster,
             palette = c("lightsalmon1", "lightskyblue"), 
             legend.title = "K means cluster",
             pointshape = 20
) +
  theme_minimal()

# pcaOnOnlyTcga -----------------------------------------------------------

cpmZCountsBrcanessTcga <- cpmZCountsBrcaness %>% 
  as.data.frame() %>% 
  dplyr::filter(str_detect(rownames(cpmZCountsBrcaness), "TCGA"))
  
table(rownames(cpmZCountsBrcanessTcga) == brcanessClinClusTcga$Sample)
TRUE 
73

# pca
pcaBrcanessTcga <- PCA(cpmZCountsBrcanessTcga)

fviz_eig(pcaBrcanessTcga, 
         addlabels = TRUE, 
         ylim = c(0,100))

fviz_pca_ind(pcaBrcanessTcga,
             axes = c(1, 2),
             # axes = c(2, 3),
             # axes = c(3, 4),
             geom = "point",
             col.ind = brcanessClinClusTcga$kMeansCluster,
             pointshape = 19,
             palette = c("1" = "pink", "2" = "palevioletred"),
             legend.title = "Kmeans cluster",
             addEllipses = TRUE
             # legend.title = "Cohort"
) +
  theme_minimal()

fviz_pca_ind(pcaBrcanessTcga,
             axes = c(1, 2),
             # axes = c(2, 3),
             # axes = c(3, 4),
             geom = "point",
             col.ind = brcanessClinClusTcga$BRCAness,
             pointshape = 19,
             palette = c("HRP" = "pink", "HRD" = "palevioletred"),
             # legend.title = "Cohort"
) +
  theme_minimal()

fviz_cluster(kMeans, data = cpmZCountsBrcaness)

# heatmap -----------------------------------------------------------------

brcanessSigContribMat <- brcanessSigContrib %>% 
  dplyr::select(hgnc_symbol, kClus1, kClus2) %>% 
  column_to_rownames(var = "hgnc_symbol")

Heatmap(t(brcanessSigContribMat), 
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        name = " Cluster center",
        row_labels = c("1", "2"))

# violinPlots -------------------------------------------------------------

longTcgaCounts <- cpmZCountsBrcanessTcga %>% 
  rownames_to_column(var = "Sample") %>% 
  pivot_longer(cols = -Sample,
               names_to = "ENSEMBL",
               values_to = "transScaleCount") %>% 
  left_join(., brcanessClinClus, by = "Sample")

# BRCA1
longTcgaCounts %>% 
  filter(ENSEMBL == "ENSG00000012048") %>% 
  ggplot(., aes(x = kMeansCluster, y = transScaleCount, colour = kMeansCluster)) +
  geom_violin() +
  geom_point() +
  xlab("K means cluster") +
  ylab("BRCA1: logCPM + z-score 
       scaled counts") +
  labs(colour="K means") +
  scale_color_manual(values =  c("1" = "pink","2" = "palevioletred"),
                     guide = "none")+
  theme_minimal() +
  stat_kruskal_test(label.x.npc = 0.6,
                    label.y.npc = 0.9)

# BRCA2
longTcgaCounts %>% 
  filter(ENSEMBL == "ENSG00000139618") %>% 
  ggplot(., aes(x = kMeansCluster, y = transScaleCount, colour = kMeansCluster)) +
  geom_violin() +
  geom_point() +
  xlab("K means cluster") +
  ylab("BRCA2: logCPM + z-score 
       scaled counts") +
  labs(colour="K means") +
  scale_color_manual(values =  c("1" = "pink","2" = "palevioletred"),
                     guide = "none")+
  theme_minimal() +
  stat_kruskal_test(label.x.npc = 0.5,
                    label.y.npc = 0.96)

# quantifyIncorrectLabels -------------------------------------------------

table(brcanessClinClus$kMeansCluster, brcanessClinClus$hrStatus)

# HRD HRP
# 1  15  22
# 2  19  17

brcanessClinClusTcga <- brcanessClinClus %>% 
  dplyr::filter(cohort == "TCGA") %>% 
  mutate(BRCAness = case_when(kMeansCluster == 1 ~ "HRP",
                              kMeansCluster == 2 ~ "HRD")) %>% 
  mutate(hrStatusNum = case_when(hrStatus == "HRP" ~ 1,
                                 hrStatus == "HRD" ~ 2)) %>% 
  dplyr::select(!`kMeans$cluster`) %>% 
  dplyr::select(!cohort)
table(brcanessClinClusTcga$BRCAness, brcanessClinClusTcga$hrStatus)
# HRD HRP
# HRD  19  17
# HRP  15  22

# metrics -----------------------------------------------------------------

table(brcanessClinClus$hrStatus, brcanessClinClus$kMeansCluster)
#     1  2
# HRD 15 19
# HRP 22 17
# Therefore...
#     HRD  HRP
# HRD 19 15
# HRP 17 22

# accuracy: 
accBRCAnessTcgaOvTest <- (19 + 22) / (15 + 17 + 19 + 22)
accBRCAnessTcgaOvTest
# [1] 0.5616438

# precision:
precisionBRCAnessTcgaOvTest <- 19 / (19 +15)
precisionBRCAnessTcgaOvTest
# [1] 0.5588235

# recall:
recallBRCAnessTcgaOvTest <- 19 / (19 + 17)
recallBRCAnessTcgaOvTest
# [1] 0.5277778

# misclassification error:
missClassBRCAnessTcgaOvTest <- (17 + 15) / (15 + 17 + 19 + 22)
missClassBRCAnessTcgaOvTest
# [1] 0.4383562

# adjRandIndex ------------------------------------------------------------

# calculate the adjusted rand index
adj.rand.index(brcanessClinClusTcga$hrStatusNum, brcanessClinClusTcga$kMeansCluster)
# [1] 0.001397143

# saveTcgaTestResults -----------------------------------------------------

brcanessClinClusTcga %>% 
  write_csv("~/identifiHR/benchmarking/brcanessGeneSig/identifiHR_BRCAnessGeneSigPredictions_tcgaOvTest_20240731.csv")

cpmZCountsBrcaness %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% 
  write_csv("~/identifiHR/benchmarking/brcanessGeneSig/brcanessGeneCountsProcessed_20240731.csv")
