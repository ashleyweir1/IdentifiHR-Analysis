# background --------------------------------------------------------------

# contrasting the results of IdentifiHR, with those of the BRCAness and MtuliscaleHRD gene expresison signatures in the TCGA-OV testing cohort.
# UPDATE to this script: 
# added expHRD results
# updated testing metrics so that HRD is the positive class

# loadPackages ------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(FactoMineR)
library(ggpubr)
library(caret)
library(ComplexHeatmap)
library(VennDiagram)

# loadData ----------------------------------------------------------------

# IdentifiHRResult --------------------------------------------------------

gBRCA1 <- read_csv("TCGA/brcaStatus/germlineBrca1_tcgaOV_integPaperSupp8_20230301.csv") %>% 
  dplyr::rename(Case.ID = "Case.ID.x")
predictTcgaOvTest <- read_csv("identifiHR/benchmarking/predicitons/identifiHR_predictionsTcgaOv_20240317.csv") %>% 
  dplyr::rename(identifiHR = s1.x) %>% 
  left_join(., gBRCA1, by = "Case.ID")

tcgaOvProcessedCountsTestingCases <- read.csv("~/identifiHR/IdentifiHR/data/cleanCohorts/tcgaOvProcessedCountsTestingCases.csv") %>% 
  column_to_rownames(var = "Sample")
identifiHRModelGenes <- read_csv("~/identifiHR/modelBackup/identifiHR_modelGenes209Annotated.csv")
identifiHRModelGenesEns <- identifiHRModelGenes$ENSEMBL
identifiHRCounts <- tcgaOvProcessedCountsTestingCases[ ,colnames(tcgaOvProcessedCountsTestingCases) %in% identifiHRModelGenesEns]

# identifiHRMetrics -------------------------------------------------------

table(predictTcgaOvTest$identifiHR, predictTcgaOvTest$hrStatus)
#     HRD HRP
# HRD  29   6
# HRP   5  33

# accuracy: 
accIdentifiHRTcgaOvTest <- (29 + 33) / (29 + 6 + 5 + 33)
accIdentifiHRTcgaOvTest
# [1] 0.8493151

# precision:
precisionIdentifiHRTcgaOvTest <- 29 / (29 + 6)
precisionIdentifiHRTcgaOvTest
# [1]  0.8285714

# recall:
recallIdentifiHRTcgaOvTest <- 29 / (29 + 5)
recallIdentifiHRTcgaOvTest
# [1] 0.8529412

# misclassification error:
missClassIdentifiHRTcgaOvTest <- (5 + 6) / (29 + 6 + 5 + 33)
missClassIdentifiHRTcgaOvTest
# [1] 0.1506849

# correlationWithHrdScoreIdentifiHR ---------------------------------------

colCorr <- as.factor(predictTcgaOvTest$correctLabel)
ggplot(predictTcgaOvTest, aes(x = HRD, y = (1 - s1.y))) +
  geom_point(aes(colour = colCorr)) +
  geom_vline(xintercept = 42, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_minimal() +
  ylab("Probability") +
  xlab("HRD score") +
  scale_colour_manual(values = c("correct" = "yellowgreen",
                                 "incorrect"= "darkorange3"),
                      name = "Model prediction") +
  stat_cor(method = "pearson",
           label.x.npc = 0.6)

# BRCAnessResult ----------------------------------------------------------

brcanessPredict <- read.csv("~/identifiHR/benchmarking/brcanessGeneSig/identifiHR_BRCAnessGeneSigPredictions_tcgaOvTest_20240731.csv")
brcanessCounts <- read_csv("identifiHR/benchmarking/brcanessGeneSig/brcanessGeneCountsProcessed_20240731.csv") %>% 
  column_to_rownames(var = "Sample")

table(brcanessPredict$BRCAness, brcanessPredict$hrStatus)
#     HRD HRP
# HRD  19  17
# HRP  15  22

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

# MultiscaleHRDResult -----------------------------------------------------

multiscaleHRD <- read_csv("identifiHR/benchmarking/multiscaleHRD/identifiHR_benchmarkingMultiscaleHRD_20240317.csv")
multiscaleHRDCounts <- read_csv("identifiHR/benchmarking/multiscaleHRD/multiscaleHRDGeneCountsProcessed_20240315.csv")%>% 
  column_to_rownames(var = "Sample")

table(multiscaleHRD$multiscaleHrd, multiscaleHRD$hrStatus) 
#           HRD HRP
# HRDmulti  26  14
# HRPmulti   8  25

# accuracy: 
accMutliscaleHRDTcgaOvTest <- (25 + 26) / (25 + 8 + 14 +26)
accMutliscaleHRDTcgaOvTest
# [1] 0.6986301

# precision:
precisionMutliscaleHRDTcgaOvTest <- (26) / (26 + 14)
precisionMutliscaleHRDTcgaOvTest
# [1] 0.65

# recall:
recallMutliscaleHRDTcgaOvTest <- (26) / (26 + 8)
recallMutliscaleHRDTcgaOvTest
# [1] 0.7647059

# misclassification error:
missClassMutliscaleHRDTcgaOvTest <- (8 + 14) / (25 + 8 + 14 +26)
missClassMutliscaleHRDTcgaOvTest
# [1] 0.3013699

multiscaleHRDSub <- multiscaleHRD %>% 
  dplyr::select(c(Sample,corHRD, multiscaleHrd)) %>% 
  mutate(MultiscaleHRD = case_when(multiscaleHrd == "HRDmulti" ~ "HRD",
                                   multiscaleHrd == "HRPmulti" ~ "HRP"))

# expHRDResult ------------------------------------------------------------

expHRD <- read_delim("identifiHR/benchmarking/expHRD/tcgaOvTest/tcgaOvTest.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE) %>% 
  dplyr::rename(Sample = 'Patient ID') %>% 
  left_join(., tcgaTestClin, by = "Sample") %>% 
  dplyr::mutate(hrStatusExpHrd = case_when(predicted_HRD < 42 ~ "HRP",
                                           TRUE ~ "HRD")) %>% 
  dplyr::mutate(correctLabel = case_when(hrStatusExpHrd == hrStatus ~ "correct",
                                         hrStatusExpHrd != hrStatus ~ "incorrect"))

expHRDCounts <- read_csv("/home/users/allstaff/weir.a/identifiHR/benchmarking/expHRD/tcgaOvTest/identifiHR_benchmarkingAgainstExpHRD_deSeq2NormCountsTcgaOvTest_20240724.csv")

# subset counts to only expHRD model genes - note that only n = 348 were present in the data
expHrdGene <- read_csv("identifiHR/benchmarking/expHRD/FINAL_98_GENE_ENSG_ID_repo.csv", 
                       col_names = FALSE) %>% 
  dplyr::filter(X1 == "HRD_GENE_TOTAL") %>% 
  column_to_rownames(var = "X1") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(expHRD = "expHRD") %>% 
  dplyr::rename(ENSEMBL = HRD_GENE_TOTAL) %>% 
  left_join(., expHRDCounts, by = "ENSEMBL") %>% 
  dplyr::filter(expHRD == "expHRD") %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  dplyr::select(-expHRD) %>% 
  column_to_rownames(var = "ENSEMBL") %>% 
  t()

table(expHRD$hrStatusExpHrd, expHRD$hrStatus) 
#     HRD HRP
# HRD  34  35
# HRP   0   4

# accuracy: 
accExpHRDTcgaOvTest <- (34 + 4) / (34 + 35 + 0 + 4)
accExpHRDTcgaOvTest
# [1] 0.5205479

# precision:
precisionExpHRDTcgaOvTest <- (34) / (34 + 35)
precisionExpHRDTcgaOvTest
# [1] 0.4927536

# recall:
recallExpHRDTcgaOvTest <- (34) / (34 + 0)
recallExpHRDTcgaOvTest
# [1] 1

# misclassification error:
missClassExpHRDTcgaOvTest <- (0 + 35) / (34 + 35 + 0 + 4)
missClassExpHRDTcgaOvTest
# [1] 0.4794521

# joinTables --------------------------------------------------------------

metrics <- as.data.frame(c("accuracy", "precision", "recall", "misclassificationError")) 
colnames(metrics)[1] <- "metrics"
metrics$IdentifiHRMetrics <- c(accIdentifiHRTcgaOvTest, 
                               precisionIdentifiHRTcgaOvTest,
                               recallIdentifiHRTcgaOvTest,
                               missClassIdentifiHRTcgaOvTest)

metrics$BRCAnessMetrics <- c(accBRCAnessTcgaOvTest, 
                             precisionBRCAnessTcgaOvTest,
                             recallBRCAnessTcgaOvTest,
                             missClassBRCAnessTcgaOvTest)

metrics$multiscaleHRDMetrics <- c(accMutliscaleHRDTcgaOvTest, 
                                  precisionMutliscaleHRDTcgaOvTest,
                                  recallMutliscaleHRDTcgaOvTest,
                                  missClassMutliscaleHRDTcgaOvTest)

metrics$expHRDMetrics <- c(accExpHRDTcgaOvTest, 
                           precisionExpHRDTcgaOvTest,
                           recallExpHRDTcgaOvTest,
                           missClassExpHRDTcgaOvTest)

# reshapeMetrics ----------------------------------------------------------

modelOrder <- c("IdentifiHR", "BRCAness", "MultiscaleHRD", "ExpHRD")
metricOrder <- c("accuracy", "precision", "recall" , "misclassificationError")
metricsLong <- pivot_longer(metrics,
                            names_to = "model",
                            values_to = "metricValue",
                            cols = c("IdentifiHRMetrics", "BRCAnessMetrics",  "multiscaleHRDMetrics", "expHRDMetrics")) %>% 
  mutate(Model = case_when(model == "IdentifiHRMetrics" ~ "IdentifiHR",
                           model == "BRCAnessMetrics" ~ "BRCAness", 
                           model == "multiscaleHRDMetrics" ~ "MultiscaleHRD",
                           model == "expHRDMetrics" ~ "ExpHRD")) %>% 
  mutate(Model = factor(Model, levels = modelOrder)) %>% 
  mutate(Metrics = factor(metrics, levels = metricOrder))

# plotMetrics -------------------------------------------------------------

ggplot(metricsLong, aes(fill = Model, y = metricValue, x = Metrics)) + 
  geom_bar(position="dodge", stat="identity", colour = "black") +
  ylim(c(0,1.0)) +
  ylab("Measure") +
  xlab("Metric") +
  theme_minimal() +
  scale_fill_manual(values = c("yellowgreen", "palevioletred", "orchid4", "dodgerblue4"))+
  theme(legend.position = "bottom")

# joinPredictions ---------------------------------------------------------

predictions <- left_join(predictTcgaOvTest, brcanessPredict, by = "Sample") %>% 
  left_join(., multiscaleHRDSub, by = "Sample") %>% 
  left_join(., expHRD, by = "Sample") %>% 
  mutate(germlineBRCA1 = case_when(is.na(germlineBRCA2.x) & is.na(somaticBRCA1.x) & is.na(somaticBRCA2.x) & brcaStatus.x == "mutatedBRCA" ~ "gBRCA1"))

sampPredictions <- predictions %>% 
  dplyr::select(c(Sample, HRD.x, hrStatus.x, germlineBRCA1, germlineBRCA2.x, somaticBRCA1.x, somaticBRCA2.x, brcaStatus.x,  identifiHR, BRCAness, MultiscaleHRD, hrStatusExpHrd)) %>% 
  arrange(HRD.x) %>% 
  dplyr::rename(IdentifiHR = identifiHR) %>%  
  dplyr::rename(ExpHRD = hrStatusExpHrd) %>% 
  dplyr::mutate(germlineBRCA = case_when(germlineBRCA1 == "gBRCA1"~ "germlineBRCA",
                                         germlineBRCA2.x == "gBRCA2" ~ "germlineBRCA",
                                         TRUE ~ "wildtypeBRCA"),
                somaticBRCA = case_when(somaticBRCA1.x == "sBRCA1"~ "somaticBRCA",
                                        somaticBRCA2.x == "sBRCA2" ~ "somaticBRCA",
                                         TRUE ~ "wildtypeBRCA"))

# heatmap -----------------------------------------------------------------

hrdColFun <- colorRamp2(c(0, 42, 100), c("turquoise", "lightgrey", "salmon"))
haIdentifiHR <- HeatmapAnnotation('HRD Score' = sampPredictions$HRD.x,
                                  'HR Status' = sampPredictions$hrStatus.x,
                                  # 'BRCA Status' = sampPredictions$brcaStatus.x,
                                  # 'gBRCA1' = sampPredictions$germlineBRCA1,
                                  # 'gBRCA2' = sampPredictions$germlineBRCA2.x,
                                  # 'sBRCA1' = sampPredictions$somaticBRCA1.x,
                                  # 'sBRCA2' = sampPredictions$somaticBRCA2.x,
                                  "gBRCA1/2" = sampPredictions$germlineBRCA,
                                  "sBRCA1/2" = sampPredictions$somaticBRCA,
                                  IdentifiHR = sampPredictions$IdentifiHR,
                                  BRCAness = sampPredictions$BRCAness,
                                  MultiscaleHRD = sampPredictions$MultiscaleHRD,
                                  ExpHRD = sampPredictions$ExpHRD,
                                  col = list('HRD Score' = hrdColFun,
                                             'HR Status' = c("HRP" = "turquoise", "HRD" = "salmon"), 
                                             'gBRCA1/2' = c("germlineBRCA" = "goldenrod1", "wildtypeBRCA" = "navajowhite1"),
                                             'sBRCA1/2' = c("somaticBRCA" = "lightsteelblue4", "wildtypeBRCA" = "lightsteelblue1"),
                                             IdentifiHR = c("HRP" = "yellowgreen", "HRD" = "olivedrab"),
                                             BRCAness = c("HRP" = "pink", "HRD" = "palevioletred"),
                                             MultiscaleHRD = c("HRP" = "thistle", "HRD" = "orchid4"),
                                             ExpHRD = c("HRP" = "slategray3", "HRD" = "dodgerblue4")),
                                  gp = gpar(col = "black"))


draw(haIdentifiHR)

# pcaIdentifiHR -----------------------------------------------------------

sampPredictionsSor <- sampPredictions %>% 
  column_to_rownames(var = "Sample")
identifiHRCountsIdx <- rownames(identifiHRCounts)
identifiHRSampPred <- sampPredictionsSor[identifiHRCountsIdx, ]
table(rownames(identifiHRCounts) == rownames(identifiHRSampPred))
# TRUE 
# 73 
identifiHRPca <- PCA(identifiHRCounts)
fviz_pca_ind(identifiHRPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = identifiHRSampPred$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(identifiHRPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = identifiHRSampPred$hrStatus.x,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

fviz_pca_ind(identifiHRPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = identifiHRSampPred$IdentifiHR,
             pointshape = 19,
             palette = c("HRP" = "yellowgreen", "HRD" = "olivedrab"),
             legend.title = "Model status"
) +
  theme_minimal() 

# pcaBRCAness -------------------------------------------------------------

brcanessCountsIdx <- rownames(brcanessCounts)
brcanessSampPred <- sampPredictionsSor[brcanessCountsIdx, ]
table(rownames(brcanessCounts) == rownames(brcanessSampPred))
# TRUE 
# 73 
options(ggrepel.max.overlaps = Inf)
brcanessPca <- PCA(brcanessCounts)
fviz_pca_ind(brcanessPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = brcanessSampPred$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(brcanessPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = brcanessSampPred$hrStatus.x,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

fviz_pca_ind(brcanessPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = brcanessSampPred$BRCAness,
             pointshape = 19,
             palette = c("HRP" = "pink", "HRD" = "palevioletred"),
             legend.title = "Model status"
) +
  theme_minimal() 

# pcaMultiscaleHRD --------------------------------------------------------

multiscaleHrdCountsIdx <- rownames(multiscaleHRDCounts)
multiscaleHrdSampPred <- sampPredictionsSor[multiscaleHrdCountsIdx, ]
table(rownames(multiscaleHRDCounts) == rownames(multiscaleHrdSampPred))
# TRUE 
# 73 
options(ggrepel.max.overlaps = Inf)
multiscaleHRDPca <- PCA(multiscaleHRDCounts)
fviz_pca_ind(multiscaleHRDPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = multiscaleHrdSampPred$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(multiscaleHRDPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = multiscaleHrdSampPred$hrStatus.x,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

fviz_pca_ind(multiscaleHRDPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = multiscaleHrdSampPred$MultiscaleHRD,
             pointshape = 19,
             palette = c("HRP" = "thistle", "HRD" = "orchid4"),
             legend.title = "Model status"
) +
  theme_minimal() 

# pcaExpHRD ---------------------------------------------------------------

expHRDCountsIdx <- rownames(expHrdGene)
expHRDSampPred <- sampPredictionsSor[expHRDCountsIdx, ]
table(rownames(expHrdGene) == rownames(expHRDSampPred))
# TRUE 
# 73 
options(ggrepel.max.overlaps = Inf)
expHRDPca <- PCA(expHrdGene)
fviz_pca_ind(expHRDPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = expHRDSampPred$HRD,
             pointshape = 19,
             legend.title = "HRD score"
) +
  theme_minimal() +
  scale_color_gradient2(low = "darkturquoise", 
                        high = "tomato",
                        mid = "salmon",
                        midpoint = 70) 

fviz_pca_ind(expHRDPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = expHRDSampPred$hrStatus.x,
             pointshape = 19,
             palette = c("tomato", "darkturquoise"),
             legend.title = "HR status"
) +
  theme_minimal() 

fviz_pca_ind(expHRDPca,
             axes = c(1, 2),
             geom = "point",
             col.ind = expHRDSampPred$ExpHRD,
             pointshape = 19,
             palette = c("HRP" = "slategray3", "HRD" = "dodgerblue4"),
             legend.title = "Model status"
) +
  theme_minimal() 

# upsetPlot ---------------------------------------------------------------

identifiHRGenes <- colnames(identifiHRCounts)
brcanessGenes <- colnames(brcanessCounts)
multiscaleHRDGenes <- colnames(multiscaleHRDCounts)
expHRDGenes <- colnames(expHrdGene)

Reduce(intersect,list(identifiHRGenes, brcanessGenes, multiscaleHRDGenes, expHRDGenes))
# [1] "ENSG00000140525"

binMat <- list(identifiHRGenes = identifiHRGenes,
               brcanessGenes = brcanessGenes,
               multiscaleHRDGenes = multiscaleHRDGenes,
               expHRDGenes = expHRDGenes)
library(ComplexHeatmap)
set.seed(NULL)
set.seed(220310)
combMat <- make_comb_mat(identifiHR = binMat$identifiHRGenes, 
                         BRCAness = binMat$brcanessGenes,
                         multiscaleHRD = binMat$multiscaleHRDGenes,
                         expHRD = binMat$expHRDGenes)

combMat
# Top 8 combination sets are:
#   identifiHR BRCAness multiscaleHRD expHRD code size
#                                         x 0001  321
#                                  x        0010  204
#           x                               1000  180
#                    x                      0100   34
#           x                             x 1001   18
#                                  x      x 0011   11
#           x                      x        1010    9
#                    x                    x 0101    3
# 
# Sets are:
#           set size
# identifiHR  209
# BRCAness   40
# multiscaleHRD  228
# expHRD  356

UpSet(combMat, 
      set_order = c("identifiHR", "BRCAness", "multiscaleHRD", "expHRD"), 
      comb_order = order(comb_size(combMat)),
      right_annotation = upset_right_annotation(combMat, 
                                                gp = gpar(fill = c("yellowgreen", "palevioletred", "orchid4", "dodgerblue4"),
                                                          annotation_name_side = "top",
                                                          axis_param = list(side = "top"))),
      top_annotation = upset_top_annotation(combMat,
                                            gp = gpar(fill = c("black","black", "black", "black","black","black", "black", "black", "yellowgreen", "palevioletred", "orchid4", "dodgerblue4")),
                                            annotation_name_rot = 90,
                                            annotation_name_side = "right",
                                            axis_param = list(side = "right")))

# vennDiagram -------------------------------------------------------------

venn <- venn.diagram(binMat, 
                     filename = NULL,
                     fill = c("yellowgreen", "palevioletred", "orchid4", "dodgerblue4"))
grid.draw(venn)
