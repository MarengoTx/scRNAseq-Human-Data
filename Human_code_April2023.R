
#### ------------------ R code for the Human manuscript ------------------ ####

options(stringsAsFactors = F, max.print = 100)
library(SingleCellExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(ggExtra)
library(SingleR)
library(Matrix)
library(patchwork)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
library(scales)
library(scater)
library(dplyr)
library(cowplot)
library(circlize)
library(enrichplot)
library(ggridges)
library(Rmagic)
library(ggpmisc)
library(glmnet)
library(foreach)
library(ggrepel)
library(EnvStats)
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)
library(ggplot.multistats)
library(ggpubr)
library(survminer)
library(harmony)
library(knitr)
library(dplyr)
library(data.table)
library(kableExtra)
library(scRepertoire)
library(SCINA)
library(ProjecTILs)
library(SAVER)
library(future)

c25 <- c("dodgerblue2",
         "#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2",
         "#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# load rds data ####
load("sc_data_human_as_of_03072023.RData")


sc_ctr <- Read10X(data.dir = "filtered_feature_bc_matrix/")

sc_data <- CreateSeuratObject(counts = sc_ctr$`Gene Expression`, )
sc_data[["Protein"]] <- CreateAssayObject(counts = sc_ctr$`Antibody Capture`)


sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-", assay = "RNA")
# map correct IDs ####
# 1 = 4
# 2 = 5
# 3 = 9
# 4 = 64
# 5 = 89
# 6 = 93
old_IDs <- sapply(rownames(sc_data@meta.data), function(xx) strsplit(xx, "-")[[1]][2])
sc_data@meta.data$PatientID <- ifelse(old_IDs == "1", "Pt4",
                                      ifelse(old_IDs == "2", "Pt5",
                                             ifelse(old_IDs == "3", "Pt9",
                                                    ifelse(old_IDs == "4", "Pt64",
                                                           ifelse(old_IDs == "5", "Pt89",
                                                                  ifelse(old_IDs == "6", "Pt93","other"))))))
sc_data@meta.data$barcode_ptID <- paste0(sapply(rownames(sc_data@meta.data), function(xx) strsplit(xx, "-")[[1]][1]), "-1_", sc_data$PatientID)
table(sc_data@meta.data$PatientID)

# import UMI-treatment data ####
UMI_files <- list.files(".", 
                        pattern = "HashtagAssignment.csv")
UMI_data <- rbindlist(lapply(UMI_files, function(xx) read.table(xx, sep = ",", header = T) %>% 
                               mutate(PatientID = paste0("Pt",strsplit(xx, "-")[[1]][1]))))
UMI_data$barcode_ptID <- paste0(UMI_data$barcode, "_", UMI_data$PatientID)
table(UMI_data$PatientID,
      UMI_data$AHH)
head(UMI_data)
all(UMI_data$barcode == row.names(sc_data@meta.data))
which(UMI_data$barcode != row.names(sc_data@meta.data))


df_test <- data.frame(cbind(UMI_data, 
                            barcode_sc_data=gsub("-[0-9]$","-1",row.names(sc_data@meta.data)),
                            sc_data@meta.data[,c("PatientID", "barcode_ptID")]))

df_test <- data.frame(merge(data.frame(cbind(rownames(sc_data@meta.data),
                                             sc_data@meta.data)), 
                            UMI_data, 
                            by="barcode_ptID"))
df_test$barcode_sc_data <- gsub("-[0-9]$","-1",colnames(sc_data))

df_test$matching <- ifelse(df_test$barcode %in% df_test$barcode_sc_data, T, F)
table(df_test$matching, df_test$AHH)
rownames(df_test) <- df_test$rownames.sc_data.meta.data.
df_test <- df_test[match(colnames(sc_data@assays$RNA@counts), row.names(df_test)),]
all(colnames(sc_data@assays$RNA@counts) == colnames(sc_data@assays$Protein@counts))


sc_data@meta.data <- df_test

table(sc_data@meta.data$PatientID.x,
      sc_data@meta.data$AHH)

rownames(sc_data)[grep("il2", rownames(sc_data), ignore.case = T)]


# read ID from cloupe ####
lib_ID_cloupe <- read.table("~/Downloads/LibraryID.csv", sep=",", header = T)
head(lib_ID_cloupe)
agg_cloupe <- read.table("../cloupe_files_01172022/t-SNE-Projections_agg.csv", sep=",", header = T, row.names = 1)
agg_cloupe <- merge(agg_cloupe, sc_data@meta.data, by="row.names")
head(agg_cloupe)
table(agg_cloupe$PatientID,
      agg_cloupe$AHH)
FeatureScatter(sc_data,
               feature1 = "ADT-CD25-TotalSeqC",
               feature2 = "ADT-CD8-TotalSeqC",
               group.by = "AHH")
View(cbind(colnames(sc_data@assays$Protein@data),
           sc_data@meta.data))
ggplot(data.frame(cbind(sc_data@meta.data,
                        CD8 = log10(sc_data@assays$Protein@data["ADT-CD8-TotalSeqC",]),
                        IL2RA = log10(sc_data@assays$Protein@data["ADT-CD25-TotalSeqC",]))) %>%
         filter(percent.mt < 5) %>%
         filter(AHH == "AHH04") %>%
         arrange(percent.mt),
       aes(x = CD8,
           y = IL2RA)) +
  geom_point(size=0.125) +
  theme_bw() +
  geom_vline(xintercept = 2, col="red") +
  geom_hline(yintercept = 2, col="red") +
  facet_wrap(~ PatientID.x) +
  stat_quadrant_counts(xintercept = 2, yintercept = 2)


ggplot(data.frame(cbind(sc_data@meta.data,
                        CD8 = log10(sc_data@assays$Protein@data["ADT-CD8-TotalSeqC",]),
                        IL2RA = log10(sc_data@assays$Protein@data["ADT-CD25-TotalSeqC",]))) %>%
         filter(percent.mt < 5) %>%
         filter(AHH == "AHH04") %>%
         arrange(percent.mt),
       aes(x = CD8,
           y = IL2RA)) +
  geom_point(size=0.125) +
  theme_bw() +
  geom_vline(xintercept = 2, col="red") +
  geom_hline(yintercept = 2, col="red") +
  facet_wrap(~ PatientID.x) +
  stat_quadrant_counts(xintercept = 2, yintercept = 2)

# filter the data based on UMI count and mitochondrial DNA content
sc_data_2 <- subset(sc_data, cells = colnames(sc_data)[which(sc_data@meta.data$nFeature_RNA > 300 & sc_data@meta.data$nFeature_RNA < 4000 & sc_data@meta.data$percent.mt < 5)])
table(sc_data_2$PatientID.x,
      sc_data_2$AHH)

ggplot(data.frame(cbind(sc_data_2@meta.data,
                        CD8 = log10(sc_data_2@assays$Protein@data["ADT-CD8-TotalSeqC",]),
                        IL2RA = log10(sc_data_2@assays$Protein@data["ADT-CD25-TotalSeqC",]))) %>%
         filter(percent.mt < 5) %>%
         filter(AHH == "AHH04") %>%
         arrange(percent.mt),
       aes(x = CD8,
           y = IL2RA)) +
  geom_point(size=0.125) +
  theme_bw() +
  geom_vline(xintercept = 2, col="red") +
  geom_hline(yintercept = 2, col="red") +
  facet_wrap(~ PatientID.x) +
  stat_quadrant_counts(xintercept = 2, yintercept = 2)

table(sc_data_no_int$PatientID.1,
      sc_data_no_int$AHH)
ggplot(data.frame(cbind(sc_data_no_int@meta.data,
                        CD8 = log10(sc_data_no_int@assays$Protein@counts["ADT-CD8-TotalSeqC",]),
                        IL2RA = log10(sc_data_no_int@assays$Protein@counts["ADT-CD25-TotalSeqC",]))) %>%
         filter(percent.mt < 5) %>%
         filter(AHH == "AHH04") %>%
         arrange(percent.mt),
       aes(x = CD8,
           y = IL2RA)) +
  geom_point(size=0.125) +
  theme_bw() +
  geom_vline(xintercept = 2, col="red") +
  geom_hline(yintercept = 2, col="red") +
  facet_wrap(~ PatientID.1) +
  stat_quadrant_counts(xintercept = 2, yintercept = 2)
View(cbind(colnames(sc_data_no_int@assays$Protein@data),
           sc_data_no_int@meta.data))
View(sc_data_2@meta.data)
sc_data_2@meta.data$barcode_sc_data <- NULL
sc_data_2@meta.data$matching <- NULL

#   -    ####
# >> from here re-run using corrected annotation << ####
#   -    ####
sc_data_2 <- subset(sc_data_2, cells=colnames(sc_data_2)[which(sc_data_2$AHH %in% c("AHH01", "AHH03", "AHH04"))])
nrow(sc_data_2@meta.data) == nrow(sc_data_no_int@meta.data)

sc_data_no_int <- sc_data_2 %>%
  SCTransform(vars.to.regress = "percent.mt", verbose = T, variable.features.n = 15000) %>%
  RunPCA(verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.25) %>%
  FindClusters(resolution = 0.125) %>%
  FindClusters(resolution = 0.05)

qc_plot <- DropletUtils::barcodeRanks(sc_data@assays$RNA@counts[,
                                                                which(colnames(sc_data@assays$RNA@counts) %in% rownames(sc_data@meta.data)[which(sc_data$PatientID.x == "Pt64")])
])
plot(qc_plot$rank, qc_plot$total, log="xy")

empty_plot <- DropletUtils::emptyDrops(sc_data@assays$RNA@counts)

DimPlot(sc_data_no_int,
        group.by = "AHH")

DimPlot(sc_data_no_int,
        group.by = "SCT_snn_res.0.25",
        cols = c25) +
  guides(color=guide_legend(ncol=2, 
                            override.aes = list(size=5)))

DimPlot(sc_data_no_int,
        group.by = "PatientID.x")

FeaturePlot(sc_data_no_int,
            features = c("CD8A", "ADT-CD8-TotalSeqC", 
                         "CD4","ADT-CD4-TotalSeqC",
                         "IL7R", "KLF2", "LEF1", "GZMB", "CCL5", "LAG3", "NKG7", "IL2RA"),
            min.cutoff = "q25",
            order = T)

sc_data_no_int$AHH_cluster <- paste0(sc_data_no_int$seurat_clusters, "_", sc_data_no_int$AHH)
Idents(sc_data_no_int) <- sc_data_no_int$SCT_snn_res.0.25
VlnPlot(sc_data_no_int, 
        features = c(genes2test, "FOXP3", "TRBV13", "ZAP70", "SELL", "CCR7"),
        slot = "scale.data",
        stack = T,
        pt.size = 0) 

set.seed(1)
input <- sc_data_no_int@assays$RNA@counts
ind_heat <- sample(1:ncol(input), 5000, replace = F)

ind_rows <- which(rownames(input) %in% VariableFeatures(sc_data_no_int)[1:20])
heat <- input[ind_rows,ind_heat]

rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
heat[1:10,1:10]
colnames(heat) <- cls
rownames(heat) <- rws
df_anno_cols <- data.frame(Clusters = sc_data_no_int$SCT_snn_res.0.25[ind_heat],
                           AHH = sc_data_no_int$AHH[ind_heat])
df_anno_cols <- df_anno_cols[order(as.numeric(as.character(df_anno_cols$Clusters)), decreasing = F),]
df_anno_cols[order()]
Heatmap(heat[,match(rownames(df_anno_cols),colnames(heat))],
        show_row_names = T,
        show_column_names = F,
        clustering_distance_rows = "pearson",
        cluster_columns = F,
        col = colorRamp2(breaks = c(-2,0,1),
                         colors = c(muted("blue"), scales::alpha("lightgray",0.3), muted("red"))),
        top_annotation = HeatmapAnnotation(df = df_anno_cols,
                                           col = list(AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        use_raster = T)


DoHeatmap(sc_data_no_int,
          features = genes2test,
          group.by = "seurat_clusters",
          # group.by = "AHH_cluster",
          disp.min = -1, disp.max = 1,
          draw.lines = F) + 
  scale_fill_gradientn(colors = c("cyan", "white", "tomato"))
table(sc_data$AHH)
table(sc_data$PatientID)


# split CD8 and CD4 positive objects ####
VlnPlot(sc_data_no_int, 
        assay = "SCT",
        slot = "scale.data",
        log = T,
        features = c("CD8A", "CD4"),
        pt.size = 0)
VlnPlot(sc_data_no_int, 
        assay = "Protein",
        slot = "data",
        log = T,
        features = c("ADT-CD8-TotalSeqC", "ADT-CD4-TotalSeqC"),
        pt.size = 0)

FeatureScatter(sc_data_no_int, 
               slot = "data",
               # shuffle = T,
               feature1 = "ADT-CD3-TotalSeqC",
               feature2 = "ADT-CD8-TotalSeqC",
               pt.size = 1) +
  xlim(0,100)


# use protein data ####
rownames(sc_data_no_int@assays$Protein@counts)
ggplot(NULL, 
       aes(x=log1p(
         as.numeric(
           sc_data_no_int@assays$Protein@counts["ADT-CD3-TotalSeqC",])),
         y=log1p(
           as.numeric(
             sc_data_no_int@assays$Protein@counts["ADT-CD8-TotalSeqC",])),
         col=log1p(
           as.numeric(
             # sc_data_no_int@assays$Protein@counts["ADT-CD8-TotalSeqC",])),
             sc_data_no_int@assays$SCT@data["CD8B",] + sc_data_no_int@assays$SCT@data["CD8A",])),
       )) +
  geom_point(size=0.25) +
  theme_classic() +
  labs(color="CD8A+CD8B") + 
  scale_color_gradient(low = "cyan", high = "red") +
  xlab("CD3") + ylab("CD8") +
  geom_hline(yintercept = 5.1, linetype="dashed")


ggplot(NULL, 
       aes(x=log1p(
         as.numeric(
           sc_data_no_int@assays$Protein@counts["ADT-CD3-TotalSeqC",])),
         y=log1p(
           as.numeric(
             sc_data_no_int@assays$Protein@counts["ADT-CD4-TotalSeqC",])),
         col=log1p(
           as.numeric(
             sc_data_no_int@assays$SCT@data["CD8B",] + sc_data_no_int@assays$SCT@data["CD8A",])),
       )) +
  geom_point(size=0.25) +
  theme_classic() +
  labs(color="CD8A+CD8B") + 
  scale_color_gradient(low = "cyan", high = "red") +
  xlab("CD3") + ylab("CD4") +
  geom_hline(yintercept = 5.1, linetype="dashed")


#prep CD8 ####
sc_cd8 <- subset(sc_data_no_int, cells = colnames(sc_data_no_int)[which(sc_data_no_int@assays$Protein@counts["ADT-CD8-TotalSeqC",] > 5.1)]) %>%
  SCTransform(vars.to.regress = "percent.mt", verbose = T, variable.features.n = 10000) %>%
  RunPCA(verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.145) %>%
  FindClusters(resolution = 0.045) 


# start from here after maping UMIs ####
genes2test <- c("IL2RA", "KLF2", "TBX21", "EOMES", "TOX", "TOX2", "BATF3", "BATF", "IRF4", "IRF8", "PRF1", "MKI67")
genes2test <- unique(c(genes2test, c("CD8A", "CD4", "KLRB1", "IL7R", "CCR7", "SELL", "CD27", "CD69", "ITGAL", "CD44", "GZMB", "GZMA", "GZMK", "IFNG", "TNF", "CD28")))

rownames(sc_data)[grep("ki67", rownames(sc_data), ignore.case = T)]
Heatmap(table(sc_data_no_int$SCT_snn_res.0.25,sc_data_no_int$AHH),
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(breaks = c(0,50,400,2000),
                         colors = c("white", "#F8F8FF", "#AFEEEE", "#4169E1")), 
        border = T, 
        row_title = "Clusters", column_title = "Treatment Groups", name="# cells")

VlnPlot(sc_data, 
        group.by = "seurat_clusters",
        # features = c("GZMB", "IFNG", "PRF1"),
        features = genes2test,
        split.by = "AHH", 
        log=T,
        stack = T, flip = F,
        pt.size = 0) + 
  geom_hline(yintercept = c(1.5,2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), col="#DC143C", size=1, linetype="dashed")

sc_data@assays$Protein@data[,1:4]
FeatureScatter(sc_data_no_int,
               feature1 = "ADT-CD8-TotalSeqC",
               feature2 = "ADT-CD4-TotalSeqC",
               # group.by = "seurat_clusters",
               plot.cor = F) +
  facet_wrap(~ sc_data_no_int$SCT_snn_res.0.25, ncol = 4)

DoHeatmap(sc_data_no_int,
          features = genes2test, 
          group.bar = T, 
          group.by = "AHH", 
          slot = "data",
          disp.min = -1, disp.max = 1,
          draw.lines = F) + 
  scale_fill_gradientn(colors = c("cyan", "white", "tomato"))

sc_data_no_int$ClusterTreatment <- paste0(sc_data_no_int$SCT_snn_res.0.25, "_", sc_data_no_int$AHH)
avgexp <- AverageExpression(sc_data_no_int, 
                            features = genes2test,
                            group.by = "ClusterTreatment",
                            return.seurat = F)
avgexp <- avgexp$SCT
cls <- colnames(avgexp)
rws <- rownames(avgexp)
heat <- t(apply(avgexp, 1, scale))
colnames(heat) <- cls
heat <- t(heat)
order_rows <- as.numeric(sapply(rownames(heat), function(xx) strsplit(xx, "_")[[1]][1]))
heat <- heat[order(order_rows, decreasing = F),]
Heatmap(heat, 
        cluster_rows = F,
        cluster_columns = F,
        border = T, 
        col = colorRamp2(breaks = c(-2,0,2), colors = c("#40E0D0", "white", "#DC143C")),
        split = order_rows[order(order_rows, decreasing = F)],
        name = "Expression")

# score cell type (grossly) with SCINA ####
rownames(sc_data_no_int)[grep("IL12", rownames(sc_data_no_int), ignore.case = T)]
cd8_tcm <- c("CD8A", "CD8B","CCR7", "SELL", "CD27")
cd8_tem <- c("CD8A", "CD8B","CX3CR1", "IL7R")
cd4_tcm <- c("CD4", "SELL", "CCR7", "IL2RA")
treg <- c("CD4", "FOXP3")
cd8_eff <- cd8t <- c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "CD69")
cd8_exhausted <- c("CD8A", "CD8B", "PDCD1", "CTLA4", "LAG3")
th1 <- c("CD4", "IFNG", "TNF", "STAT4", "IL12RB2", "IL12A")
cd8_trm <- c("CD8A", "CD8B", "CXCR3", "CXCR6", "CCR5", "CCR6", "VCAM1", "ITGAE")
th17 <- c("CD4", "IL17A", "IL22", "IL25")
th2 <- c("CD4", "IL4", "IL5", "IL13", "STAT6")
max_length <- 9
df_signatures <- data.frame(CD8_TCM = c(cd8_tcm,
                                        rep(NA, max_length - length(cd8_tcm))),
                            CD8_TEM = c(cd8_tem,
                                        rep(NA, max_length - length(cd8_tem))),
                            CD4_TCM = c(cd4_tcm,
                                        rep(NA, max_length - length(cd4_tcm))),
                            TREG = c(treg,
                                     rep(NA, max_length - length(treg))),
                            CD8_EFF = c(cd8_eff,
                                        rep(NA, max_length - length(cd8_eff))),
                            CD8_EXHAUSTED = c(cd8_exhausted,
                                              rep(NA, max_length - length(cd8_exhausted))),
                            TH1 = c(th1,
                                    rep(NA, max_length - length(th1))),
                            CD8_TRM = c(cd8_trm,
                                        rep(NA, max_length - length(cd8_trm))),
                            TH17 = c(th17,
                                     rep(NA, max_length - length(th17))),
                            TH2 = c(th2,
                                    rep(NA, max_length - length(th2))))

write.table(df_signatures, "SCINA_cell_signatures.csv", na = "", sep=",", row.names = F, quote = F)
ref_immune_markers <- preprocess.signatures("SCINA_cell_signatures.csv")

rna4scina <- sc_data_no_int@assays$SCT@data  
rna4scina <- data.matrix(rna4scina[which(rownames(rna4scina) %in% unique(unlist(ref_immune_markers))),])
# rna4scina <- log2(rna4scina+1)
# rna4scina[1:10,1:5]
exp_scina <- preprocessCore::normalize.quantiles(rna4scina)
rownames(exp_scina) <- rownames(rna4scina)
colnames(exp_scina) <- colnames(rna4scina)

res_scina <- SCINA(exp_scina, 
                   ref_immune_markers, 
                   max_iter = 100, 
                   convergence_n = 10, 
                   convergence_rate = 0.999, 
                   sensitivity_cutoff = 0.9, 
                   rm_overlap=TRUE, 
                   allow_unknown=TRUE, 
                   log_file='SCINA.log')
table(res_scina$cell_labels)
head(res_scina$cell_labels)
colnames(sc_data_no_int@meta.data)
sc_data_no_int$SCINA <- res_scina$cell_labels

sc_data_no_int@meta.data[,32:41] <- t(res_scina$probabilities)


Heatmap(table(sc_data_no_int$SCINA, sc_data_no_int$TILPRED_compiled),
        col = colorRamp2(breaks = c(0,50,100,1000),
                         colors = c("white", "white", "lightgray", "black")),
        border = T)

head(sc_data_no_int)
colnames(sc_data@meta.data)

DimPlot(sc_data_no_int,
        reduction = "aUMAP",
        group.by = "SCINA", 
        cols = c25, 
        pt.size = 0.5,
        split.by = "AHH")


Heatmap(table(sc_data$seurat_clusters,sc_data$SCINA),
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(breaks = c(0,10,500,2000),
                         colors = c("white", "#F8F8FF", "#AFEEEE", "#4169E1")), 
        border = T, 
        row_title = "Clusters", column_title = "Predicted Cell Type")


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

ggplot(gg_t,
       aes(x=TBX21*BATF, 
           col=Group)) +
  stat_ecdf(geom = "step", size=1.25) +
  theme_bw() +
  labs(col = "Group") +
  ylab("cumulative frequency") +
  ggtitle("") +
  facet_wrap(~ CellType)


# run projecTILs ####
options(future.globals.maxSize = 8000 * 1024^2)
plan("multiprocess", workers = 6)
print(plan())

ref <- load.reference.map(ref = "~/Google Drive/Protocols_bfx/ref_TILAtlas_mouse_v1.rds")
table(ref$functional.cluster)
table(ref$Subclass_Cell_Identity)
# ref$functional.cluster <- ref$Subclass_Cell_Identity
DefaultAssay(sc_data_no_int) <- "SCT"
query.projected <- make.projection(sc_data_no_int, ref = ref, ncores = 6, future.maxSize = 8000 * 1024^2)
plot.projection(ref, query.projected, pointsize = 0.5,linesize = 0.7)
query.projected <- cellstate.predict(ref = ref, query = query.projected)
table(query.projected$functional.cluster)
pred_tils <- query.projected@meta.data
colnames(pred_tils)
pred_tils <- pred_tils[,c(19:23)]
head(pred_tils)

head(cbind(sc_data_no_int@meta.data, pred_tils[match(rownames(sc_data_no_int@meta.data), rownames(pred_tils)),]))

sc_data_no_int@meta.data <- data.frame(cbind(sc_data_no_int@meta.data, 
                                             row_pred_tils = rownames(pred_tils)[match(rownames(sc_data_no_int@meta.data), rownames(pred_tils))],
                                             pred_tils[match(rownames(sc_data_no_int@meta.data), rownames(pred_tils)),4:5]))
colnames(sc_data_no_int@meta.data)[20:21] <- c("TIL_ref.cluster", "TIL_ref.cluster_conf")
DimPlot(sc_data_no_int,
        group.by = "TIL_ref.cluster",
        cells = colnames(sc_data_no_int)[which(sc_data_no_int$TIL_ref.cluster_conf > 0.9)],
        cols=c25) +
  facet_wrap(~ sc_data_no_int$AHH[which(sc_data_no_int$TIL_ref.cluster_conf > 0.9)])


# >> normalize protein data ####
DefaultAssay(sc_data_no_int) <- 'Protein'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
sc_data_no_int@assays$Protein@counts[1:6,1:4]
ggplot(NULL, aes(x=sc_data_no_int@assays$Protein@data["HTO-AHH03-TotalSeqC",],
                 y=sc_data_no_int@assays$Protein@data["HTO-AHH04-TotalSeqC",],
                 z=sc_data_no_int@assays$Protein@data["HTO-AHH01-TotalSeqC",],
                 col=sc_data_no_int$AHH)) +
  geom_point(size=0.2) +
  theme_half_open()

VariableFeatures(sc_data_no_int) <- rownames(sc_data_no_int[["Protein"]])[-c(1:4)]
sc_data_no_int <- NormalizeData(sc_data_no_int, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'prot_PCA') %>%
  RunUMAP(features = VariableFeatures(sc_data_no_int), reduction.name = 'prot_UMAP', reduction = "prot_PCA") %>%
  FindNeighbors(reduction = "prot_PCA", dims = 1:20) %>% 
  FindClusters(resolution = 0.25)
Loadings(sc_data_no_int[["prot_PCA"]])
DimPlot(sc_data_no_int,
        reduction  = "prot_UMAP",
        cols = c25)  
DimPlot(sc_data_no_int,
        reduction  = "prot_UMAP",
        group.by = "AHH") 

Heatmap(table(sc_data_no_int$SCT_snn_res.0.25,sc_data_no_int$Protein_snn_res.0.25),
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(breaks = c(0,25,100,1000),
                         colors = c("white", "#F8F8FF", "#AFEEEE", "#4169E1")), 
        border = T, 
        row_title = "Clusters RMA", column_title = "Clusters Proteins", name="# cells")
DimPlot(sc_data_no_int %>%
          subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$TIL.cluster_conf > 0.9)]),
        group.by = "TIL.cluster",
        cols = c25)

DefaultAssay(sc_data_no_int) <- "SCT"
rownames(sc_data_no_int)[grep("tnf",rownames(sc_data_no_int), ignore.case = T)]

FeaturePlot(sc_data_no_int,
            reduction = "prot_UMAP",
            features = c("CD8A", "CD4"), 
            order = T,
            min.cutoff = "q10")

markers_prot <- FindAllMarkers(sc_data_no_int,
                               min.pct = 0.25, 
                               features = rownames(sc_data_no_int)[-c(1:4)],
                               logfc.threshold = 0.25,
                               assay = "Protein",
                               slot = "data", 
                               only.pos = T)


rownames(sc_data@meta.data)
VlnPlot(sc_data_no_int, 
        group.by = "TIL.cluster",
        assay = "Protein", slot = "scale.data",
        features = rownames(sc_data_no_int@assays$Protein@data)[-c(1:4)],
        split.by = "AHH", 
        log=F,
        same.y.lims = F,
        stack = T, flip = F,
        pt.size = 0) + 
  geom_hline(yintercept = c(1.5,2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), col="#DC143C", size=1, linetype="dashed")

# make heatmap with protein data ####
input <- sc_data_no_int@assays$Protein@data
ind_heat <- sample(1:ncol(input), 2000, replace = F)
heat <- input[,ind_heat]
heat <- heat[-c(1:4),]
rownames(heat) <- gsub("-TotalSeqC","",gsub("ADT-","",rownames(heat)))
heat <- as.data.frame(heat)

rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
heat[1:10,1:10]
colnames(heat) <- cls
rownames(heat) <- rws
df_anno_cols <- data.frame(Clusters = sc_data_no_int$Protein_snn_res.0.25[ind_heat],
                           AHH = sc_data_no_int$AHH[ind_heat])
df_anno_cols <- df_anno_cols[order(as.numeric(as.character(df_anno_cols$Clusters)), decreasing = F),]
table(df_anno_cols$Clusters, df_anno_cols$AHH)
table(df_anno_cols$Clusters, df_anno_cols$AHH)
for(ii in 0:9){
  temp_clust <- ii
  temp_df <- df_anno_cols[which(df_anno_cols$Clusters == ii),]
  temp_df <- temp_df[order(temp_df$AHH),]
  temp_df$UMI <- rownames(temp_df)
  if(ii == 0){
    df_anno_cols_2 <- temp_df
  } else {
    df_anno_cols_2 <- data.frame(rbind(df_anno_cols_2,
                                       temp_df))
  }
}
rownames(df_anno_cols_2) <- df_anno_cols_2$UMI
df_anno_cols_2$UMI <- NULL
Heatmap(heat[,match(rownames(df_anno_cols_2),colnames(heat))],
        show_row_names = T,
        show_column_names = F,
        clustering_distance_rows = "pearson",
        cluster_columns = F,
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c(scales::alpha("blue", 0.2), scales::alpha("lightgray",0.3), scales::alpha("red", 0.3))),
        top_annotation = HeatmapAnnotation(df = df_anno_cols_2,
                                           col = list(AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        use_raster = T,
        name = "Abundance")



#clean yp CD8 Naive-like
# keep only CD8+
FeatureScatter(sc_data_no_int, 
               slot = "data", 
               group.by = "seurat_clusters", 
               cells = colnames(sc_data_no_int)[which(sc_data_no_int$TILPRED_compiled == "CD8_NaiveLike")],
               feature1 = "ADT-CD8-TotalSeqC", 
               feature2 = "ADT-CD4-TotalSeqC",
               # feature2 = "CD8A",
               plot.cor = F,
               pt.size = 0.5, 
               col=c25) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5))) +
  geom_vline(xintercept = 1.75, linetype="dashed", size=0.75) +
  xlim(0,5)

# check for CCR7 adn CD62L
FeaturePlot(sc_data_no_int,
            reduction = "aUMAP2",
            slot = "data",
            features = c("ADT-CD8-TotalSeqC","ADT-CD197-TotalSeqC", 
                         "ADT-CD45RA-TotalSeqC", "SELL"),
            order = T, 
            min.cutoff = "q50", 
            ncol = 2)

ggplot(NULL, aes(x= sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",],
                 y = sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",])) +
  geom_vline(xintercept = 1.8, linetype="dashed", col="red") +
  geom_hline(yintercept = 0.8, linetype="dashed", col="red") +
  geom_point(size=0.25) +
  theme_classic() +
  xlab("") +
  ylab("") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  stat_quadrant_counts()

sc_data_no_int$CD8_CD45RA <- ifelse(sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] >= 1.8 &
                                      sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",] >= 0.8, 1, 0)

table(sc_data_no_int$CD8_CD45RA)
rownames(sc_data_no_int@assays$Protein@data)
ggplot(NULL, 
       aes(x= sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",which(sc_data_no_int$Protein_snn_res.0.25== 4)],
           y = sc_data_no_int@assays$SCT@scale.data["SELL",which(sc_data_no_int$Protein_snn_res.0.25== 4)])) +
  geom_point(size=0.25) +
  theme_classic() +
  xlab("") +
  ylab("") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  stat_quadrant_counts()
DimPlot(sc_data_no_int,
        reduction = "aUMAP2",
        cols=c25)
table(sc_data_no_int$TILPRED_compiled)
FeaturePlot(sc_data_no_int,
            
            reduction = "aUMAP2",
            slot = "data",
        
            features = c("ADT-CD279-TotalSeqC", "ADT-CD366-TotalSeqC","ADT-CD223-TotalSeqC",
                         "IL2RA", "IFNG", "ZAP70", "IL2", "ADT-CD127-TotalSeqC", "CD27","EOMES"),
            order = T, 
            
            min.cutoff = "q50", 
            ncol = 2)
VlnPlot(sc_data_no_int %>%
          subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$TILPRED_compiled == "CD8_NaiveLike")]),
        
        idents = c(3,8),
        group.by = "KNN",
        slot = "data",
        features = c("ADT-CD8-TotalSeqC", "ADT-CD197-TotalSeqC", "ADT-CD45RA-TotalSeqC", "SELL","GZMB",
                     "ADT-CD279-TotalSeqC", "ADT-CD366-TotalSeqC","ADT-CD223-TotalSeqC",
                     "IFNG", "ZAP70", "ADT-CD127-TotalSeqC", "CD27", "EOMES"),
        stack = T, flip = T,
        pt.size = 0)

# cells sit in 3 adn 8, run DGE for Naive Memory ####
Idents(sc_data_no_int) <- sc_data_no_int$KNN
diff_abundance_c3_c8 <- FindMarkers(sc_data_no_int, 
                                    subset.ident = c(3,8),
                                    # slot = "data",
                                    min.cells.group = 5,
                                    ident.1 = "8", ident.2 = "3",
                                    logfc.threshold = 0.2, 
                                    min.pct = 0.3)
diff_abundance_c3_c8$Gene <- rownames(diff_abundance_c3_c8)
diff_abundance_c3_c8$Direction <- ifelse(diff_abundance_c3_c8$avg_log2FC > 0 & diff_abundance_c3_c8$p_val_adj < 0.05, "UP",
                                         ifelse(diff_abundance_c3_c8$avg_log2FC < 0 & diff_abundance_c3_c8$p_val_adj < 0.05, "DOWN", "NS"))
diff_abundance_c3_c8$Labels <- ifelse(diff_abundance_c3_c8$Direction != "NS" & 
                                        diff_abundance_c3_c8$pct.1 > 0.4 , 
                                      diff_abundance_c3_c8$Gene,"")
diff_abundance_c3_c8$Labels <- ifelse(diff_abundance_c3_c8$Labels == "ACTB" | grepl("^MT-",diff_abundance_c3_c8$Labels), "", diff_abundance_c3_c8$Labels)
head(diff_abundance_c3_c8)
ggplot(diff_abundance_c3_c8, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_abundance_c3_c8$Labels,
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

# add CD8 naive to metadata ####
sc_data_no_int$Manual_compiled <- ifelse(sc_data_no_int$KNN == 3, "CD8_Naive_Memory", 
                                         ifelse(sc_data_no_int$KNN == 8, "CD8_EffectorLike_Memory","other"))

# check CD4 Naive ####
FeaturePlot(sc_data_no_int,
            reduction = "aUMAP2",
            slot = "data",
            features = c("ADT-CD4-TotalSeqC", "ADT-CD197-TotalSeqC", "ADT-CD45RA-TotalSeqC", "SELL"),
            order = T, 
            min.cutoff = "q50", 
            ncol = 2)

# use Kmeans clustering to get the optimal K ####
library(factoextra)
tot_withinss <- map_dbl(1:20,  function(k){
  model <- kmeans(x = t(sc_data_no_int@assays$Protein@scale.data), centers = k)
  model$tot.withinss
})

# Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:20,
  tot_withinss = tot_withinss
)

# Plot the elbow plot
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:20)

sc_data_no_int$KMeans <- kmeans(x = t(sc_data_no_int@assays$Protein@scale.data), centers = 15)$cluster
DimPlot(sc_data_no_int,
        reduction = "prot_UMAP",
        group.by = "KMeans", 
        label = T,label.box = T,repel = T,label.color = "white",
        cols=c25) +
  guides(color=guide_legend(ncol=2, 
                            override.aes = list(size=2)))


VlnPlot(sc_data_no_int,
        group.by = "KNN", 
        features = c("ADT-CD8-TotalSeqC","ADT-CD4-TotalSeqC"),
        stack = T,
        pt.size = 0)

# cluster KNN 5 and 7 seem to be the two populations of CD4 naive
table(sc_data_no_int$KNN, 
      sc_data_no_int$AHH)

Idents(sc_data_no_int) <- sc_data_no_int$KNN
diff_abundance_c14_c5 <- FindMarkers(sc_data_no_int, 
                                     subset.ident = c(14,5),
                                     min.cells.group = 5,
                                     ident.1 = "5", ident.2 = "14",
                                     logfc.threshold = 0.2, 
                                     min.pct = 0.3)
diff_abundance_c14_c5$Gene <- rownames(diff_abundance_c14_c5)
diff_abundance_c14_c5$Direction <- ifelse(diff_abundance_c14_c5$avg_log2FC > 0 & diff_abundance_c14_c5$p_val_adj < 0.05, "UP",
                                          ifelse(diff_abundance_c14_c5$avg_log2FC < 0 & diff_abundance_c14_c5$p_val_adj < 0.05, "DOWN", "NS"))
diff_abundance_c14_c5$Labels <- ifelse(diff_abundance_c14_c5$Direction != "NS" & 
                                         diff_abundance_c14_c5$pct.1 > 0.4 , 
                                       diff_abundance_c14_c5$Gene,"")
diff_abundance_c14_c5$Labels <- ifelse(diff_abundance_c14_c5$Labels == "ACTB" | grepl("^MT-",diff_abundance_c14_c5$Labels), "", diff_abundance_c14_c5$Labels)
head(diff_abundance_c14_c5)
ggplot(diff_abundance_c14_c5, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_abundance_c14_c5$Labels,
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

# make feature plot with top DEGs
FeaturePlot(sc_data_no_int,
            reduction = "aUMAP",
            slot = "data",
            features = c("GZMB", "IFNG", "IL7R", "KLF2", "CD74", "CD44"),
            order = T, 
            min.cutoff = "q75", 
            ncol = 2)
VlnPlot(sc_data_no_int,
        slot = "data",
        idents = c(5,14),
        features = c("GZMB", "IFNG","CD74", "IL7R", "KLF2",  "CD44"),
        group.by = "KNN",
        pt.size = 0,
        ncol=3)

# add labels for clusters 5 and 7 to metadata ####
sc_data_no_int$Manual_compiled <- ifelse(sc_data_no_int$KNN == 14, "CD4_Naive_Memory", 
                                         ifelse(sc_data_no_int$KNN == 5, "CD4_EffectorLike_Memory", sc_data_no_int$Manual_compiled))

# check Treg population ####
ggplot(NULL, 
       aes(x= sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",which(sc_data_no_int@assays$SCT@data["FOXP3",] >= 0.6931472)],
           y = sc_data_no_int@assays$Protein@data["ADT-CD25-TotalSeqC",which(sc_data_no_int@assays$SCT@data["FOXP3",] >= 0.6931472)],
           col=sc_data_no_int@assays$SCT@data["FOXP3",which(sc_data_no_int@assays$SCT@data["FOXP3",] >= 0.6931472)])) +
  geom_vline(xintercept = 1.65, linetype="dashed", col="red") +
  geom_hline(yintercept = 2, linetype="dashed", col="red") +
  geom_point(size=0.75) +
  theme_classic() +
  xlab("CD4prot") +
  ylab("CD25/IL2RAprot") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  stat_quadrant_counts(col="red") +
  scale_color_gradient2(name="FOXP3",low = "blue", mid = "gray",high = "red", midpoint = 1)

# see how many cells express FOXP3
sum(sc_data_no_int@assays$SCT@data["FOXP3",] > 0)
sc_data_no_int$Manual_compiled <- ifelse(sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.65 &
                                           sc_data_no_int@assays$Protein@data["ADT-CD25-TotalSeqC",] > 2 &
                                           sc_data_no_int@assays$SCT@data["FOXP3",] >= 0.6931472, "Treg", 
                                         ifelse(sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.65 &
                                                  sc_data_no_int@assays$Protein@data["ADT-CD25-TotalSeqC",] < 2 &
                                                  sc_data_no_int@assays$SCT@data["FOXP3",] >= 0.6931472, "Treg_NC",sc_data_no_int$Manual_compiled))

FeaturePlot(sc_data_no_int,
            reduction = "aUMAP2",
            slot = "data",
            features = c("ADT-CD127-TotalSeqC", "ADT-CD197-TotalSeqC","ADT-CD278-TotalSeqC","ADT-CD196-TotalSeqC", "MKI67"),
            order = T, 
            min.cutoff = "q75", 
            ncol = 2)
VlnPlot(sc_data_no_int,
        slot = "data",
        idents = c(4,14),
        split.by = "AHH",
        features = c("IL17A", "IL10", "IL12A"),
        group.by = "KNN",
        pt.size = 0,
        ncol=3)

# add labels for clusters 4 and 14 Treg to metadata ####
sc_data_no_int$Manual_compiled <- ifelse(sc_data_no_int$KNN == 4, "CD4_Exhausted", 
                                         ifelse(sc_data_no_int$KNN == 14, "CD8_Exhausted", sc_data_no_int$Manual_compiled))
table(sc_data_no_int$Manual_compiled,
      sc_data_no_int$AHH)

# check Effector cycling ####
ggplot(as.data.frame.matrix(table(sc_data_no_int$TILPRED_compiled,
                                  sc_data_no_int$KNN)) %>%
         mutate(CellType = rownames(.)) %>%
         melt() %>%
         filter(CellType == "Eff Cycling"),
       aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("KNN clusters") +
  ylab("num Eff Cycling cells")


# check CD8 Effector memory
# neg/low: ccr7, sell, cd27, cd45ra
# pos: cd8, cd127, cx3cr1
View(as.data.frame.matrix(table(sc_data_no_int$KNN,
                                sc_data_no_int$AHH)))
rownames(sc_data_no_int)[grep("il7",rownames(sc_data_no_int), ignore.case = T)]
VlnPlot(sc_data_no_int,
        slot = "data",
        idents = c(13, 10),
        split.by = "AHH",
        features = c("ADT-CD8-TotalSeqC", "ADT-CD279-TotalSeqC","ADT-CD223-TotalSeqC","ADT-CD366-TotalSeqC",
                     "ADT-CD45RA-TotalSeqC", "ADT-CD197-TotalSeqC", "CCR7", "SELL", "CD27", #neg
                     "IL7R", "ADT-CD127-TotalSeqC", #pos
                     "NFKB1", "IL7R", "IFNG", "CD69","KLRK1",
                     "PRF1", "GZMA", "GZMK"),
        group.by = "KNN",
        stack = T,
        pt.size = 0,
        ncol=4)


# start selecting CCR7/SELL negative CD8 cells
sc_data_no_int$CD8pCD45RAn <- ifelse(sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] > 2 &
                                       sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",] < 0.8, 1, 0)
ggplot(NULL, 
       aes(
         # x= sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",which(sc_data_no_int$CD8pCD45RAn == 1)],
         x= sc_data_no_int@assays$SCT@scale.data["SELL",which(sc_data_no_int$CD8pCD45RAn == 1)],
         y = sc_data_no_int@assays$Protein@data["ADT-CD197-TotalSeqC",which(sc_data_no_int$CD8pCD45RAn == 1)])
) +
  
  geom_point(size=0.125) +
  theme_classic() +
  xlab("CD62L/SELL") +
  ylab("CCR7") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  
  annotation_logticks() 
# 

ind_random <- sample(1:ncol(sc_data_no_int), 2000, replace = F)
heat <- sc_data_no_int@assays$SCT@scale.data[which(rownames(sc_data_no_int@assays$SCT@scale.data) %in% c("LAG3","HAVCR2","CTLA4","GZMB","CD28", "FAS",
                                                                                                         "CCR7", "SELL", "CD27", "IL2RA","ITGAL","CD58",
                                                                                                         "NFKB1", "IFNG", "CD69","KLRK1","B3GAT1",
                                                                                                         "PRF1", "ZAP70", "LCK")),
                                             which(sc_data_no_int$KNN == 14)]

heat[1:4,1:4]
rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat, 1, scale))
colnames(heat)<- cls
rownames(heat) <- gsub("ADT-","",gsub("-TotalSeqC","",rownames(heat)))
df_anno_cols <- data.frame(Clusters = sc_data_no_int$KNN[which(sc_data_no_int$KNN == 14)],
                           AHH = sc_data_no_int$AHH[which(sc_data_no_int$KNN == 14)])
df_anno_cols <- df_anno_cols[order(as.numeric(as.character(df_anno_cols$Clusters)), decreasing = F),]

Heatmap(heat, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        top_annotation = HeatmapAnnotation(df = df_anno_cols,
                                           col = list(AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        show_column_names = F,
        col=colorRamp2(breaks = c(-2,-0,2),
                       colors = c(alpha("cyan", 0.2),
                                  alpha("lightgray", 0.1),
                                  "#FA8072")),
        name = "Abundance")

# let's look for effectors separately ####
ggplot(NULL, 
       aes(
         x= sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",],
         y= sc_data_no_int@assays$SCT@scale.data["ZAP70",],
         # y = sc_data_no_int@assays$Protein@data["ADT-CD137L-TotalSeqC",],
         col = sc_data_no_int$AHH)
) +
 
  geom_point(size=0.125) +
  theme_classic() +
  xlab("") +
  ylab("") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))

FeaturePlot(sc_data_no_int,
            reduction = "aUMAP2",
            slot = "data",
            features = c("IFNG","GZMB", "PRF1", "ZAP70", "LAG3", "MKI67","CD8A", "PDCD1", "CTLA4"),
            order = T, 
            # pt.size = 1.2,
            min.cutoff = "q75", 
            ncol = 3)

# add labels for Effector to metadata ####
sc_data_no_int$Manual_compiled <- ifelse(sc_data_no_int$KNN == 13, "CD8_Effector",
                                         ifelse(sc_data_no_int$KNN == 10, "CD8_Effector_Exhausted",sc_data_no_int$Manual_compiled))
View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled,
                                sc_data_no_int$AHH)))
View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled,
                                sc_data_no_int$KNN)))
View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled,
                                sc_data_no_int$TIL.cluster)))

# check central memory cells ####
# neg: cd69, cd103, cd45ra
# pos: ccr7, cd62l, s1pr1, cd25,CD127/il7r, tcf1, eomes
# some: cxcr5 (tfh)
Idents(sc_data_no_int)
FeaturePlot(sc_data_no_int, 
            reduction = "aUMAP2",
            features = c("ADT-CD69-TotalSeqC", "ADT-CD103-TotalSeqC", "ADT-CD45RA-TotalSeqC",
                         "ADT-CD197-TotalSeqC", "SELL", "ADT-CD25-TotalSeqC", "ADT-CD127-TotalSeqC", "TCF7", "EOMES"),
            min.cutoff = "q25",
            ncol = 4)
RidgePlot(sc_data_no_int, 
          sort = F, 
          group.by = "KNN",
          # group.by = "Manual_compiled",
          stack = F,same.y.lims = T,
          features = c("ADT-CD45RA-TotalSeqC",
                       "ADT-CD197-TotalSeqC", 
                       "SELL",
                       "FAS"),
          slot = "data") +
  geom_vline(xintercept = 0.35)

# pick CD45RA negative cells
ggplot(NULL, 
       aes(
         x= log10(sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",]),
         # y= log10(sc_data_no_int@assays$Protein@data["ADT-CD69-TotalSeqC",]),
         y = log10(sc_data_no_int@assays$SCT@scale.data["SELL",]+5),
         # col = factor(sc_data_no_int$KNN)
       )
) +

  geom_point(size=0.025, alpha=0.5) +
  theme_classic() +
  xlab("") +
  ylab("")


ggplot(NULL, 
       aes(
         y= 0,
         x = sc_data_no_int@assays$SCT@data["GZMB",]
       )
) +
  geom_density_ridges(scale=1.5) +

  theme_classic() +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = 0.875)
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.25]", "(0.25, 0.75]", "(0.75, 1]")
  )

temp_cm <- ifelse(sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",] < 0.7 &
                    sc_data_no_int@assays$Protein@data["ADT-CD197-TotalSeqC",] > 0.35 &
                    sc_data_no_int@assays$SCT@scale.data["SELL",]+5 > 4.75 &
                    sc_data_no_int@assays$Protein@data["ADT-CD69-TotalSeqC",] < 1.3 &
                    sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] > 2 &
                    sc_data_no_int@assays$SCT@data["GZMB",] < 0.875, "CD8_CM",
                  ifelse(sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",] < 0.7 &
                           sc_data_no_int@assays$Protein@data["ADT-CD197-TotalSeqC",] > 0.35 &
                           sc_data_no_int@assays$SCT@scale.data["SELL",]+5 > 4.75 &
                           sc_data_no_int@assays$Protein@data["ADT-CD69-TotalSeqC",] < 1.3 &
                           sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.8 &
                           sc_data_no_int@assays$SCT@data["GZMB",] < 0.875, "CD4_CM", "others"))
table(temp_cm, sc_data_no_int$Manual_compiled)
table(temp_cm, sc_data_no_int$AHH)
sc_data_no_int$Manual_compiled_CM <- ifelse(temp_cm != "others", temp_cm, sc_data_no_int$Manual_compiled)

Idents(sc_data_no_int) <- sc_data_no_int$Manual_compiled_CM
VlnPlot(sc_data_no_int, 
        sort = F, 
        stack = F,
        features = c("ADT-CD45RA-TotalSeqC",
                     "ADT-CD197-TotalSeqC", 
                     "ADT-CD69-TotalSeqC",
                     "ADT-CD127-TotalSeqC",
                     "CD27","CD28",
                     "SELL",
                     "FAS",
                     "IFNG",
                     "GZMB",
                     "ZAP70",
                     "ADT-CD137L-TotalSeqC"),
        slot = "data",
        ncol = 4,
        pt.size = 0,
        cols=c25)
table(sc_data_no_int$AHH,
      sc_data_no_int$Manual_compiled_CM)

# check other cell types ####
table(sc_data_no_int$TILPRED_compiled)
View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled,
                                sc_data_no_int$Manual_compiled_CM)))

heat <- sc_data_no_int@assays$Protein@data[-c(1:4),which(sc_data_no_int$Manual_compiled_CM == "other")]
rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
rownames(heat) <- gsub("ADT-","",gsub("-TotalSeqC","",rownames(heat)))
heat[1:4,1:4]
df_anno_cols <- data.frame(Clusters = as.character(sc_data_no_int$KNN[which(sc_data_no_int$Manual_compiled_CM == "other")]),
                           AHH = sc_data_no_int$AHH[which(sc_data_no_int$Manual_compiled_CM == "other")])
df_anno_cols <- df_anno_cols[order(as.numeric(df_anno_cols$Clusters), decreasing = F),]

Heatmap(heat,
        show_column_names = F,
        clustering_distance_rows = "pearson",
        clustering_distance_columns =  "spearman",
        use_raster = T,
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c(scales::alpha("blue", 0.2), scales::alpha("lightgray",0.3), scales::alpha("red", 0.3))),
        top_annotation = HeatmapAnnotation(df = df_anno_cols,
                                           col = list(AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        name = "Abundance")


FeatureScatter(sc_data_no_int,
               slot = "data",
               feature1 = "ADT-CD183-TotalSeqC",
               feature2 = "ADT-CD8-TotalSeqC",
               group.by = "KNN")


VlnPlot(sc_data_no_int %>%
          subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$Manual_compiled_CM == "other")]),
        group.by = "KNN",
        slot = "scale.data",
        features = c("ADT-CD8-TotalSeqC", "ADT-CD4-TotalSeqC","ADT-CD183-TotalSeqC", "IL17A",
                     "ADT-CD352-TotalSeqC", "ADT-CD163-TotalSeqC","ADT-CD25-TotalSeqC", "ADT-CD56-TotalSeqC",
                     "ADT-CD319-TotalSeqC", "ADT-CD103-TotalSeqC",
                     "ADT-CD134-TotalSeqC","ADT-CD137L-TotalSeqC", "ADT-CD223-TotalSeqC",
                     "ADT-HLA-DR-TotalSeqC", "ADT-CD278-TotalSeqC", "ADT-CD279-TotalSeqC",
                     "ADT-CD38-TotalSeqC", "ADT-CD45RA-TotalSeqC","ADT-CD244-TotalSeqC", "ADT-CD197-TotalSeqC"),
        pt.size = 0)

# annotate the ohter ####
sc_data_no_int$Manual_compiled_CM <- ifelse(sc_data_no_int$Manual_compiled_CM == "other" &sc_data_no_int$KNN == 11, "NK_Exhausted",
                                            ifelse(sc_data_no_int$Manual_compiled_CM == "other" &sc_data_no_int$KNN == 6, "CD4_Early_Efector",
                                                   ifelse(sc_data_no_int$Manual_compiled_CM == "other" &sc_data_no_int$KNN == 15, "CD4_Exhausted",sc_data_no_int$Manual_compiled_CM)))
View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled_CM,
                                sc_data_no_int$KNN)))
View(as.data.frame.matrix(table(sc_data_no_int$TILPRED_compiled,
                                sc_data_no_int$KNN)))
# other stuff ####
rownames(sc_data_no_int@assays$Protein@scale.data)
DimPlot(sc_data_no_int,
        reduction = "aUMAP",
        group.by = "TILPRED_compiled",
        cols=c25) +
  guides(color=guide_legend(ncol=2, 
                            override.aes = list(size=2)))

Idents(sc_data_no_int) <- sc_data_no_int$KNN
VlnPlot(sc_data_no_int,
        slot = "scale.data",
        features = c("CD4", "ADT-CD4-TotalSeqC", "MKI67",
                     "GZMB", "ADT-CD279-TotalSeqC", "PRF1", 
                     "IL2RA", "ADT-CD134-TotalSeqC", "IFNG",
                     "LAMP1", "ADT-CD137L-TotalSeqC", "FOXP3"),
        ncol = 3,
        pt.size = 0)

# check TEMRA ####
sc_data_no_int$Manual_compiled_CM_TEMRA <- ifelse(sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",] > 0.7 &
                                                    sc_data_no_int@assays$Protein@data["ADT-CD197-TotalSeqC",] < 0.35 &
                                                    sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] > 2, "CD8_TEMRA",
                                                  ifelse(sc_data_no_int@assays$Protein@data["ADT-CD45RA-TotalSeqC",] > 0.7 &
                                                           sc_data_no_int@assays$Protein@data["ADT-CD197-TotalSeqC",] < 0.35 &
                                                           sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.8, "CD4_TEMRA",sc_data_no_int$Manual_compiled_CM))

View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled_CM_TEMRA,
                                sc_data_no_int$KNN)))
View(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled_CM_TEMRA,
                                sc_data_no_int$AHH)))



# CD183 = CXCR3
# CD223 = LAG3
# CD197 = CCR7
# CD366 = TIM3
# CD127 = IL7R
# CD163 > only monocytes?
# CD137L = 41BB
# CD196 = CCR6 - effector/memory T cells
# CD134 = OX40 - activates NFKB
# CD278 = ICOS
# CD319 = NK cells activation
# CD352 = promotes T differentiation into Th17

groups <- as.character(unique(sc_data_no_int$AHH))
cell_pop <- unique(sc_data_no_int$Manual_compiled_CM_TEMRA)
for(ii in 1:length(cell_pop)){
  if(ii == 1)
    res_deg_all <- list()
  
  temp.cells <- subset(sc_data_no_int, cells=which(sc_data_no_int$Manual_compiled_CM_TEMRA == cell_pop[ii]))
  DefaultAssay(temp.cells) <- "SCT"
  Idents(temp.cells) <- temp.cells$AHH
  
  tt_temp <- table(Idents(temp.cells))
  print(tt_temp)
  
  # find DEGs
  message("Finding DEGs")
  temp.res <- tryCatch(FindMarkers(temp.cells, 
                                   # slot = "data",
                                   min.cells.group = 5,
                                   ident.1 = "AHH04", ident.2 = "AHH03",
                                   logfc.threshold = 0.2, 
                                   min.pct = 0.3),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       })
  
  temp.res$Gene <- rownames(temp.res)
  temp.res$Comparison <- "AHH04_AHH03"
  res_deg_all[[ii]] <- data.frame(rbind(temp.res))
  names(res_deg_all)[ii] <- cell_pop[ii]
  
  message(" VVVVVVVVVVVVVVVVVVVVVV ")
  message(paste(" >> Done for",cell_pop[ii], "<<"))
  message(" ^^^^^^^^^^^^^^^^^^^^^^^ ")
}

# plot DEGs ####
res_all_AHH03_AHH01 <- data.table::rbindlist(res_deg_all, idcol = T)
head(res_all_AHH03_AHH01)
res_all_AHH04_AHH01 <- data.table::rbindlist(res_deg_all, idcol = T)
head(res_all_AHH04_AHH01)
res_all_AHH04_AHH03 <- data.table::rbindlist(res_deg_all, idcol = T)
head(res_all_AHH04_AHH03)
res_all <- data.frame(rbind(res_all_AHH03_AHH01,
                            res_all_AHH04_AHH01,
                            res_all_AHH04_AHH03))
colnames(res_all)[1] <- "CellType"

res_all$Significant <- ifelse(res_all$p_val_adj < 0.05, "YES", "NO")
res_all$Direction <- ifelse(res_all$Significant == "YES" & res_all$avg_log2FC > 0, "UP",
                            ifelse(res_all$Significant == "YES" & res_all$avg_log2FC < 0, "DOWN", "NS"))
res_all$size <- ifelse(res_all$Direction == "UP", 0.5,
                       ifelse(res_all$Direction == "DOWN", 0.5, 0.2))
res_all$Labels <- ifelse(!grepl("^MT-",res_all$Gene) & res_all$Direction != "NS", res_all$Gene, "") 
head(res_all)

res_all[-grep("^MT-", res_all$Gene),] %>%

  filter(Significant == "YES" & Comparison == "AHH04_AHH03") %>%
  group_by(CellType) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(Gene %in% c("CD74", "CCL5", "TIGIT", "IFNG", "TNFRSF4", "GZMB")) %>%
  mutate(size = NULL) %>%
  kbl() %>%
  kable_styling()

ggplot(res_all[-grep("^MT-", res_all$Gene),] %>%
         filter(Comparison == "AHH04_AHH03"),
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
  geom_point() + 
  theme_pubclean() +
  facet_wrap(~ CellType, scales = "free_y") +
  scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                "DOWN" = alpha("lightblue", 0.5),
                                "NS" = "lightgray")) +
  geom_vline(xintercept = 0, linetype="dashed")

table(res_all$CellType,
      res_all$Direction)

ggplot(res_all %>%
         filter(Comparison == "AHH04_AHH03" & CellType %in% c("CD8_TEMRA", "CD4_Early_Efector", "CD8_EffectorLike_Memory")),
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
  geom_point() + 
  theme_pubclean() +
  scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                "DOWN" = alpha("lightblue", 0.5),
                                "NS" = "lightgray")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(label=res_all$Labels[which(res_all$Comparison == "AHH04_AHH01" & res_all$CellType %in% c("CD8_TEMRA", "CD4_Early_Efector", "CD8_EffectorLike_Memory"))],
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black') +
  facet_wrap(~ CellType) +
  ggtitle("AHH04 vs. AHH01")

ggplot(res_all %>%
         filter(Comparison == "AHH04_AHH01" & CellType %in% c("CD4_TEMRA", "CD8_EffectorLike_Memory", "CD8_Exhausted")),
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
  geom_point() + 
  theme_pubclean() +
  scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                "DOWN" = alpha("lightblue", 0.5),
                                "NS" = "lightgray")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(label=res_all$Labels[which(res_all$Comparison == "AHH04_AHH01" & res_all$CellType %in% c("CD4_TEMRA", "CD8_EffectorLike_Memory", "CD8_Exhausted"))],
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black') +
  facet_wrap(~ CellType, scales = "free_y") +
  ggtitle("AHH04 vs. AHH01")

ggplot(res_all %>%
         filter(Comparison == "AHH03_AHH01" & CellType %in% c("CD8_TEMRA","CD4_TEMRA", "NK_Exhausted", "CD8_EffectorLike_Memory", "CD8_Naive_Memory", "CD8_Exhausted")),
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
  geom_point() + 
  theme_pubclean() +
  scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                "DOWN" = alpha("lightblue", 0.5),
                                "NS" = "lightgray")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(label=res_all$Labels[which(res_all$Comparison == "AHH03_AHH01" & res_all$CellType %in% c("CD8_TEMRA","CD4_TEMRA", "NK_Exhausted", "CD8_EffectorLike_Memory", "CD8_Naive_Memory", "CD8_Exhausted"))],
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black') +
  facet_wrap(~ CellType, scales = "free_y") +
  ggtitle("AHH03 vs. AHH01")


View(rbind(as.data.frame.matrix(table(res_all$CellType[which(res_all$Comparison == "AHH04_AHH03")],
                                      res_all$Direction[which(res_all$Comparison == "AHH04_AHH03")]))))

View(rbind(as.data.frame.matrix(table(sc_data_no_int$Manual_compiled_CM_TEMRA,
                                      sc_data_no_int$AHH))))




table(res_all$Direction[which(res_all$Comparison == "AHH03_AHH01")], 
      res_all$CellType[which(res_all$Comparison == "AHH03_AHH01")])


tt <- as.data.frame.matrix(table(res_all$Direction[which(res_all$Comparison == "AHH04_AHH03")], 
                                 res_all$CellType[which(res_all$Comparison == "AHH04_AHH03")]))[c(1,3),]
tt
Heatmap(t(tt), 
        col = colorRamp2(breaks = c(0,25,50,100),
                         colors = c(muted("cyan", l=30, c=100),
                                    muted(scales::alpha("lightgray", 0.1), l=100, c=20),
                                    muted(scales::alpha("yellow", 0.1), l=100, c=20),
                                    muted("tomato", l=5, c=50))),
        name = paste("NumDEGs","AHH04_AHH01", sep="\n"))

res_all[-grep("^MT-", res_all$Gene),] %>%
  filter(Significant == "YES" & Comparison == "AHH04_AHH03") %>%
  group_by(CellType) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(Gene %in% c("CD74", "CCL5", "TIGIT", "IFNG", "TNFRSF4", "GZMB")) %>%
  mutate(size = NULL) %>%
  kbl() %>%
  kable_styling()

VlnPlot(sc_data_no_int %>%
          subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$AHH %in% c("AHH04", "AHH03"))]),
        features = c("CD74", "CCL5", "TIGIT", "IFNG", "TNFRSF4", "GZMB"),
        group.by = "TILPRED_compiled",
        split.by = "AHH", split.plot = T,
        cols = c("lightgray", scales::alpha("red", 0.4)),
        pt.size=0)

tt <- as.data.frame.matrix(table(sc_data$TILPRED_compiled,
                                 paste0(sc_data$Protein_snn_res.0.25, "_", sc_data$AHH)))
tt$CellType <- rownames(tt)
tt <- melt(tt)
tt$Cluster <- sapply(as.character(tt$variable), function(xx) strsplit(xx, "_")[[1]][1])
tt$Treatment <- sapply(as.character(tt$variable), function(xx) strsplit(xx, "_")[[1]][2])
tt$PercentageCluster <- as.numeric(
  unlist(
    lapply(
      split(tt$value, tt$Cluster), function(xx){
        temp_perc <- round((xx/sum(xx))*100, 3)
        temp_perc
      }
    )))
head(tt)

ggplot(tt %>%
         filter(Cluster == 7), 
       aes(x=Treatment, y=PercentageCluster, fill=CellType)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values = c25) +
  facet_wrap(~ CellType, scales = "free_x", ncol=9) +
  theme(axis.text.x = element_text(angle=90),
        legend.position = "none") +
  ylim(0,35)


Idents(sc_data_no_int) <- sc_data_no_int$KNN
diff_abundance_c13_c10 <- FindMarkers(sc_data_no_int, 
                                      subset.ident = c(13,10),
                                      min.cells.group = 5,
                                      ident.1 = "10", ident.2 = "13",
                                      logfc.threshold = 0.2, 
                                      min.pct = 0.3)
diff_abundance_c13_c10$Gene <- rownames(diff_abundance_c13_c10)
diff_abundance_c13_c10$Direction <- ifelse(diff_abundance_c13_c10$avg_log2FC > 0 & diff_abundance_c13_c10$p_val_adj < 0.05, "UP",
                                           ifelse(diff_abundance_c13_c10$avg_log2FC < 0 & diff_abundance_c13_c10$p_val_adj < 0.05, "DOWN", "NS"))
diff_abundance_c13_c10$Labels <- ifelse(diff_abundance_c13_c10$Direction != "NS" & 
                                          diff_abundance_c13_c10$pct.1 > 0.4 , 
                                        diff_abundance_c13_c10$Gene,"")
diff_abundance_c13_c10$Labels <- ifelse(diff_abundance_c13_c10$Labels == "ACTB" | grepl("^MT-",diff_abundance_c13_c10$Labels), "", diff_abundance_c13_c10$Labels)
head(diff_abundance_c13_c10)
ggplot(diff_abundance_c13_c10, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_abundance_c13_c10$Labels,
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))



# make heatmap requested by Allart ####
# make heatmap with protein data ####
genes2test <- c("IL2RA", "KLF2", "TBX21", "TOX", "TOX2", "BATF3", "BATF", "IRF4", "IRF8", "PRF1", "MKI67")
genes2test <- unique(c(genes2test, c("CD8A", "CD4", "KLRB1", "IL7R", "CCR7", "SELL", "CD27", "CD69", "ITGAL", "CD44", "GZMB", "IFNG", "TNF", "CD28")))

ind_heat <- sample(1:ncol(sc_data_no_int), 5000, replace = F)
heat <- sc_data_no_int@assays$SCT@scale.data[which(rownames(sc_data_no_int@assays$SCT@scale.data) %in% genes2test),ind_heat]

rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
rownames(heat) <- rws
heat[1:10,1:10]
df_anno_cols <- data.frame(Clusters = sc_data_no_int$KNN[ind_heat],
                           AHH = sc_data_no_int$AHH[ind_heat])
table(df_anno_cols$Clusters, df_anno_cols$AHH)

for(ii in 1:max(as.numeric(df_anno_cols$Clusters))){
  temp_clust <- ii
  temp_df <- df_anno_cols[which(df_anno_cols$Clusters == ii),]
  temp_df <- temp_df[order(temp_df$AHH),
  if(ii == 1){
    df_anno_cols_2 <- temp_df
  } else {
    df_anno_cols_2 <- data.frame(rbind(df_anno_cols_2,
                                       temp_df))
  }
}
rws_df <- rownames(df_anno_cols_2)
df_anno_cols_2 <- data.frame(AHH = df_anno_cols_2[order(as.numeric(gsub("AHH","",df_anno_cols_2$AHH)), decreasing = F),2])
rownames(df_anno_cols_2) <- rws_df
head(df_anno_cols_2)

df_anno_cols_2[1:10,]

Heatmap(heat[,match(rownames(df_anno_cols_2), colnames(heat))],
        
        show_row_names = T,
        show_column_names = F,
        clustering_distance_rows = "spearman",
        cluster_columns = F,
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c(scales::alpha("blue", 0.8), scales::alpha("lightgray",0.3), scales::alpha("red", 0.8))),
        top_annotation = HeatmapAnnotation(df = df_anno_cols_2,
                                           col = list(AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        
        use_raster = F,
        name = "Abundance")

# prep table for KNN ####
tt_knn <- as.data.frame.matrix(table(paste0(sc_data_no_int$KNN),
                                     sc_data_no_int$AHH)) %>%
  mutate(PtID = sapply(rownames(.), function(xx) strsplit(xx, "_")[[1]][2])) %>%
  mutate(ClusterKNN = sapply(rownames(.), function(xx) strsplit(xx, "_")[[1]][1])) %>%
  
  kable(., format = "rst", align = "c")


head(tt_knn) 
ggplot(tt_knn, aes(x=factor(ClusterKNN, levels=1:15),
                   y=log10(value+1),
                   fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  scale_fill_manual(values = c("AHH01" = "gray",
                               "AHH03" = "green",
                               "AHH04" = "blue")) +
  facet_wrap(~ PtID) +
  theme(axis.text = element_text(size=15)) +
  xlab("KNN Cluster") + ylab("Log10(Number of Cells)") +
  labs(fill="Treatment")

# compare C10 vs. C13 ####
sc_data_no_int$KNN
sc_data_c10_c13 <- subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$KNN %in% c(10, 13))])
tt_c10_c13 <- as.data.frame.matrix(table(sc_data_c10_c13$KNN, sc_data_c10_c13$AHH)) %>%
  mutate(KNN = rownames(.)) %>% 
  melt()

VlnPlot(sc_data_c10_c13,
        features = c("ADT-CD278-TotalSeqC"),
        group.by = "KNN",
        cols = c25[14:15])

diff_abundance_c13_c10 <- FindMarkers(sc_data_c10_c13,
                                      min.cells.group = 3,
                                      ident.1 = "13", ident.2 = "10",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0.1)
diff_abundance_c13_c10$Gene <- rownames(diff_abundance_c13_c10)
diff_abundance_c13_c10$Direction <- ifelse(diff_abundance_c13_c10$avg_log2FC > 0 & diff_abundance_c13_c10$p_val_adj < 0.05, "UP",
                                           ifelse(diff_abundance_c13_c10$avg_log2FC < 0 & diff_abundance_c13_c10$p_val_adj < 0.05, "DOWN", "NS"))
diff_abundance_c13_c10$Labels <- ifelse(diff_abundance_c13_c10$Direction != "NS" & 
                                          diff_abundance_c13_c10$pct.1 > 0.4 , 
                                        diff_abundance_c13_c10$Gene,"")
diff_abundance_c13_c10$Labels <- ifelse(diff_abundance_c13_c10$Labels == "ACTB" | grepl("^MT-",diff_abundance_c13_c10$Labels), "", diff_abundance_c13_c10$Labels)
head(diff_abundance_c13_c10)
ggplot(diff_abundance_c13_c10, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_abundance_c13_c10$Labels,
                  box.padding = 0.3, 
                  force = 5,
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

rownames(sc_data_c10_c13)[grep("cx3cr", rownames(sc_data_c10_c13), ignore.case = T)]
kable(diff_abundance_c13_c10[which(diff_abundance_c13_c10$Gene %in% c("BATF3", "IRF4", "BCL2L11", "CX3CR1")),], "pipe", align = "c")


# run GSEA in groups ####
run_gsea_sc <- function(DGEresults, Signature){
  # prep signatures
  if(Signature == "TF"){
    sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c3.tft.v7.4.symbols.gmt")
  } else if(Signature == "CellType"){
    sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c8.all.v7.4.symbols.gmt")
  } else if (Signature == "allTF"){
    sig_all <- list(all_TF = read.table("~/Google Drive/Protocols_bfx/GSEA_signatures/all_TF.txt", sep="\t", header = F)[,1])
  } else {
    avail_sig <- list.files(path = "~/Google Drive/Protocols_bfx/GSEA_signatures/", 
                            pattern = ".gmt")
    avail_sig_full <- list.files(path = "~/Google Drive/Protocols_bfx/GSEA_signatures/", 
                                 pattern = ".gmt", full.names = T)
    avail_sig <- data.frame(ID = 1:length(avail_sig),
                            Pathway = avail_sig)
    print(kable(avail_sig, format = "rst", align = "c"))
    message("")
    path_sel <- as.numeric(readline("Select ID from above:"))
    sig_all <- qusage::read.gmt(avail_sig_full[path_sel])
  }
  
  # collect ranked signatures per each cell line and comparison
  temp_comparison <- DGEresults[-grep("^MT-", DGEresults$Gene, ignore.case = T),]
  
  temp_comparison <- temp_comparison[order(temp_comparison$avg_log2FC, decreasing = T),]
  
  message("Done setting up signatures")
  
  ranked_ii <- temp_comparison$avg_log2FC
  names(ranked_ii) <- temp_comparison$Gene
  
  message(paste("Running fgsea"))
  temp_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all,
                                      stats = ranked_ii),
                         error = function(xx){
                           message(xx)
                           message("\nadding NAs\n")
                           dummy_list <- NA
                           return(dummy_list)
                         })
  temp_fgsea <- temp_fgsea[order(temp_fgsea$padj, decreasing = F),]
  
  
  
  
  message(paste("  > Done running fgsea"))
  return(temp_fgsea)
}

fsgea_c10c13_tf <- run_gsea_sc(diff_abundance_c13_c10, "TF")
head(fsgea_c10c13_tf)
fsgea_c10c13_alltf <- run_gsea_sc(diff_abundance_c13_c10, "allTF")
head(fsgea_c10c13_alltf)
fsgea_c10c13_alltf$leadingEdge

fsgea_c10c13_immCells <- run_gsea_sc(diff_abundance_c13_c10, "Immune")
fsgea_c10c13_immCells %>% 
  subset(grepl("T_CELL", pathway)) %>%
  filter(padj < 0.05) %>%
  View()
head(fsgea_c10c13_immCells)

fsgea_c10c13_reactome <- run_gsea_sc(diff_abundance_c13_c10, "Reactome")
fsgea_c10c13_reactome %>% 
  filter(padj < 0.05) %>%
  View()

fgsea_cd4_c5c4_alltf <- run_gsea_sc(diff_abundance_c5_c4, "allTF")
head(fgsea_cd4_c5c4_alltf)
fgsea_cd4_c5c4_alltf$leadingEdge

fgsea_cd4_c5c4_reactome <- run_gsea_sc(diff_abundance_c5_c4, "Reactome")
fgsea_cd4_c5c4_reactome %>% 

  filter(padj < 0.05) %>%
  View()

# check TFs
all_TFs <- read.table("~/Google Drive/Protocols_bfx/GSEA_signatures/all_TF.txt", header = F, sep = "\t")[,1]
diff_TFs <- diff_abundance_c5_c4  %>%
  subset(Gene %in% all_TFs) %>%
  filter(p_val_adj < 0.05) %>%
  .$Gene
tf4heat <- unique(c(diff_TFs, unlist(fgsea_cd4_c5c4_alltf$leadingEdge)))
heat <- data.frame(sc_data_c5_c4@assays$SCT@data) %>%
  subset(rownames(.) %in% tf4heat)
cls <- colnames(heat)
rws <- rownames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- gsub("\\.","-",cls)
heat[1:10,1:4]

df_knn <- data.frame(KNN = sc_data_c5_c4$KNN)

Heatmap(heat,
        show_column_names = F,

        clustering_distance_columns = "spearman",
        clustering_distance_rows = "spearman",
        top_annotation = HeatmapAnnotation(df = df_knn,
                                           col = list(KNN=c("4" = "#9370DB",
                                                            "5" = "#00FFFF"))),
        col = colorRamp2(breaks = c(-1,0,1), 
                         colors = c("#AFEEEE", "white", "#B22222")),
        row_names_gp = gpar(fontsize=10),
        name="Expression")


# run UMAP from TF data
tf_umap <- umap(d = t(heat))
head(sc_data_c10_c13@meta.data)
head(tf_umap$layout)
all(rownames(tf_umap$layout) == rownames(sc_data_c10_c13@meta.data))
sc_data_c10_c13@meta.data <- data.frame(cbind(sc_data_c10_c13@meta.data, tf_umap$layout))
colnames(sc_data_c10_c13@meta.data)[33:34] <- c("UMAP1_TF", "UMAP2_TF")
colnames(sc_data_c10_c13@meta.data)[32] <- "Manual_compiled"
sc_data_c10_c13@meta.data$UMAP1_TF

pp <- prcomp(t(heat))
plot(pp$x[,1:2], pch=21, bg=factor(df_knn$KNN))
top_loadings_df <- data.frame(head(pp$rotation[order(abs(pp$rotation[,1]), decreasing = T),1:3], n = 30))
top_loadings_df$Gene <- rownames(top_loadings_df)
ggplot(top_loadings_df, 
       aes(x=factor(Gene, levels=rev(rownames(top_loadings_df))),
           y=abs(PC1))) +
  geom_point() +
  geom_bar(stat="identity", width = 0.1) +
  coord_flip() +
  theme_bw() +
  xlab("Top PCA TF") + ylab("PCA loadings") +
  theme(axis.text.y = element_text(size=12))

col_gene <- sc_data_c10_c13@assays$SCT@data["GZMB",]
rownames(sc_data_c10_c13@assays$Protein@data)
col_gene <- sc_data_c10_c13@assays$Protein@data["ADT-CD25-TotalSeqC",]
ggplot(NULL, aes(
  x=sc_data_c10_c13@meta.data$UMAP1_TF,
  y=sc_data_c10_c13@meta.data$UMAP2_TF,
  col=col_gene
)) +
  geom_point() +
  theme_classic() +
  labs(col="Expr") + xlab("UMAP1_TF") + ylab("UMAP2_TF") +
  scale_color_gradient2(low = "lightgray", 
                        mid = "white", 
                        high = "#DC143C", 
                        midpoint = median(col_gene))

# extrat GSEA results ####
gsea_res <- rbindlist(res_all_fgsea, idcol = T)
colnames(gsea_res)[1] <- "CellType"
gsea_res$Direction <- ifelse(gsea_res$padj < 0.05 & gsea_res$NES > 0, "UP", 
                             ifelse(gsea_res$padj < 0.05 & gsea_res$NES < 0, "DOWN", "NS"))
gsea_res$leadingEdge <- as.character(gsea_res$leadingEdge)

gsea_res <- read.table("all_GSEA_results_12152021.txt", sep = "\t", header = T)
old_gsea_res <- read.table("all_GSEA_results_use.txt", sep = "\t", header = T)

table(gsea_res$Direction[grep("CD8_Effector_gmt", gsea_res$CellType)])
table(old_gsea_res$Direction[grep("CD8_Effector_gmt", old_gsea_res$CellType)])

View(gsea_res[grep("CD8_Effector_gmt", gsea_res$CellType),])
View(old_gsea_res[grep("CD8_Effector_gmt", old_gsea_res$CellType),])

vv <- gplots::venn(list(NEW = gsea_res$pathway[grep("CD8_Effector_gmt", gsea_res$CellType)],
                        OLD= old_gsea_res$pathway[grep("CD8_Effector_gmt", old_gsea_res$CellType)]))
View(attributes(vv)[3])
table(gsea_res$CellType)

gsea_sign <- gsea_res %>%
  group_by(CellType) %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES)))
head(gsea_sign)


# run same DGE as above for Cd4 T cells ####
VlnPlot(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int@meta.data$KNN %in% c(1,4:7,14) & sc_data_no_int@assays$SCT@scale.data["CD8A",] < 0 )]),
        features = c("IFNG","GZMB", "PRF1", "ZAP70", "LAG3", "MKI67", "PDCD1", "CTLA4"),
        group.by = "KNN", 
        ncol = 2,
        cols = c25)
FeatureScatter(sc_data_no_int,
               slot = "scale.data",
               cells = colnames(sc_data_no_int)[which(sc_data_no_int@meta.data$KNN %in% c(5))],
               feature1 = "CD4",
               feature2 = "CD8A",
               group.by = "KNN"
)

# compare C5 vs. C4 ####
sc_data_no_int$KNN
sc_data_c5_c4 <- subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$KNN %in% c(5, 4))])
tt_c5_c4 <- as.data.frame.matrix(table(sc_data_c5_c4$KNN, sc_data_c5_c4$AHH)) %>%
  mutate(KNN = rownames(.)) %>% 
  melt()
VlnPlot(sc_data_c5_c4,
        features = c("ADT-CD278-TotalSeqC"),
        group.by = "KNN", 
        cols = c25[14:15])

VlnPlot(sc_data_c5_c4,
        features = c("IFNG","GZMB", "PRF1", "ZAP70", "LAG3", "MKI67","CD8A", "PDCD1", "CTLA4"),
        group.by = "KNN", 
        cols = c25[14:15])

diff_abundance_c5_c4 <- FindMarkers(sc_data_c5_c4, 
                                    
                                    min.cells.group = 3,
                                    ident.1 = "5", ident.2 = "4",
                                    logfc.threshold = 0.1, 
                                    min.pct = 0.1)
diff_abundance_c5_c4$Gene <- rownames(diff_abundance_c5_c4)
diff_abundance_c5_c4$Direction <- ifelse(diff_abundance_c5_c4$avg_log2FC > 0 & diff_abundance_c5_c4$p_val_adj < 0.05, "UP",
                                         ifelse(diff_abundance_c5_c4$avg_log2FC < 0 & diff_abundance_c5_c4$p_val_adj < 0.05, "DOWN", "NS"))
diff_abundance_c5_c4$Labels <- ifelse(diff_abundance_c5_c4$Direction != "NS" & 
                                        diff_abundance_c5_c4$pct.1 > 0.4 , 
                                      diff_abundance_c5_c4$Gene,"")
diff_abundance_c5_c4$Labels <- ifelse(diff_abundance_c5_c4$Labels == "ACTB" | grepl("^MT-",diff_abundance_c5_c4$Labels), "", diff_abundance_c5_c4$Labels)
head(diff_abundance_c5_c4)
ggplot(diff_abundance_c5_c4, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_abundance_c5_c4$Labels,
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

rownames(sc_data_c10_c13)[grep("cx3cr", rownames(sc_data_c10_c13), ignore.case = T)]
kable(diff_abundance_c5_c4[which(diff_abundance_c5_c4$Gene %in% c("BATF3", "IRF4", "BCL2L11", "CX3CR1")),], "pipe", align = "c")


# run re-clustering of AHH04 alone ####
sc_04 <- subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04")]) %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.25)

DefaultAssay(sc_04) <- "Protein"
FeaturePlot(sc_04,
            slot = "data",
            features = c("ADT-CD8-TotalSeqC", "ADT-CD4-TotalSeqC", "ADT-CD25-TotalSeqC","ADT-CD278-TotalSeqC"),
            min.cutoff = "q25",
            ncol=2)

# use proteins for 04 ####
sc_04 <- RunPCA(sc_04, reduction.name = 'pca_prot_04') %>%
  RunUMAP(features = VariableFeatures(sc_04), reduction.name = 'UMAP_prot_o4', reduction = "pca_prot_04") %>%
  FindNeighbors(reduction = "pca_prot_04", dims = 1:20) %>% 
  FindClusters(resolution = 0.25)
loadings_prot_04 <- Loadings(sc_04[["pca_prot_04"]])
loadings_prot_04[order(abs(loadings_prot_04[,1]), decreasing = T),]

DimPlot(sc_04, 
        reduction = "UMAP_prot_o4")
DimPlot(sc_04, group.by = "PatientID.x",
        reduction = "UMAP_prot_o4")

FeaturePlot(sc_04, 
            order = T,
            reduction = "UMAP_prot_o4",
            slot = "data",
            features = c("ADT-CD8-TotalSeqC", "ADT-CD4-TotalSeqC","ADT-CD278-TotalSeqC",  "ADT-CD25-TotalSeqC"),
            min.cutoff = "q25",
            ncol=2)

sc_04$KM2 <- kmeans(x = t(sc_04@assays$Protein@data), centers = 12)$cluster

DimPlot(sc_04, 
        reduction = "UMAP_prot_o4",
        group.by = "seurat_clusters",
        cols = c25)
DimPlot(sc_04, 
        reduction = "UMAP_prot_o4",
        group.by = "KM2",
        cols = c25)
DimPlot(sc_04, 
        reduction = "UMAP_prot_o4",
        group.by = "PatientID",
        cols = c25)
tt_knn_04 <- as.data.frame.matrix(table(sc_04$Protein_snn_res.0.25, sc_04$PatientID)) %>%
  mutate(Protein_clusters = rownames(.))
tt_knn_04 <- melt(tt_knn_04, id.vars="Protein_clusters")
head(tt_knn_04)
ggplot(tt_knn_04, aes(x=factor(Protein_clusters, levels = 0:6), y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  xlab("Protein_clusters") + ylab("Number of Cells") + labs(fill="PtID") +
  scale_fill_manual(values = c25)

ggplot(NULL, aes(x=factor(unique(sc_04$PatientID), levels = rev(paste0("Pt", c(4,5,64,89,9,93)))), 
                 y=as.numeric(table(sc_04$PatientID)))) +
  geom_bar(stat="identity", position="dodge", width = 0.75) +
  theme_classic() +
  xlab("Patient ID") + ylab("Number of Cells") +
  coord_flip()

ggplot(NULL, aes(x=factor(unique(sc_04$KNN), levels = 12:1), 
                 y=as.numeric(table(sc_04$KNN)))) +
  geom_bar(stat="identity", position="dodge", width = 0.75) +
  theme_classic() +
  xlab("KNN cluster") + ylab("Number of Cells") +
  coord_flip()

DefaultAssay(sc_04) <- "SCT"
VlnPlot(sc_04,
        slot = "data", 
        group.by = "KM2",
        features = c("CD8A","CD8B", "CD4", 
                     "IFNG","GZMB", "PRF1", "ZAP70", 
                     "MKI67", "CD28",  "FAS", "JUNB", "BATF3",
                     "LAG3", "PDCD1", "CTLA4", "FOXP3",
                     "CCR7", "SELL", "IL7R"),
        adjust = 1.5,
        stack = T)

FeatureScatter(sc_04,
               feature1 = "ADT-CD8-TotalSeqC", 
               feature2 = "ADT-CD4-TotalSeqC",
               group.by = "seurat_clusters",
               cols = c25, 
               plot.cor = F) +
  facet_wrap(~ sc_04$seurat_clusters, ncol=7) +
  labs(col="Seurat cluster")

VlnPlot(sc_04,
        features = c("ADT-CD8-TotalSeqC", "ADT-CD4-TotalSeqC"),
        group.by = "seurat_clusters",
        stack = T, flip = T)

Idents(sc_04) <- sc_04$seurat_clusters
top_markers_04_by_seurat <- FindAllMarkers(sc_04,
                                           only.pos = T,
                                           min.pct = 0.4,
                                           logfc.threshold = 0.4)

top_markers_04_by_knn %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(20)) %>%
  View

top20_04 <- top_markers_04_by_seurat %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(5))


# CD183 = CXCR3
# CD223 = LAG3
# CD197 = CCR7
# CD366 = TIM3
# CD127 = IL7R
# CD163 > only monocytes?
# CD137L = 41BB
# CD196 = CCR6 - effector/memory T cells
# CD134 = OX40 - activates NFKB
# CD278 = ICOS
# CD319 = NK cells activation
# CD244 = NK cells receptor
# CD352 = promotes T differentiation into Th17

heat <- data.frame(sc_04@assays$Protein@data[which(rownames(sc_04@assays$Protein@data) %in% c("ADT-CD8-TotalSeqC","ADT-CD4-TotalSeqC", "ADT-CD3-TotalSeqC", "ADT-CD56-TotalSeqC",
                                                                                              "ADT-CD137L-TotalSeqC","ADT-CD319-TotalSeqC","ADT-CD25-TotalSeqC",
                                                                                              "ADT-CD279-TotalSeqC", "ADT-CD223-TotalSeqC", "ADT-CD366-TotalSeqC",
                                                                                              "ADT-CD197-TotalSeqC", "ADT-CD278-TotalSeqC", "ADT-CD127-TotalSeqC",
                                                                                              "ADT-CD196-TotalSeqC", "ADT-CD244-TotalSeqC")),])

colnames(heat) <- gsub("\\.","-",colnames(heat))
cls <- colnames(heat)
heat <- t(apply(heat, 1, scale))
colnames(heat) <- cls
rownames(heat) <- gsub("ADT-","",gsub("-TotalSeqC","",rownames(heat)))
heat[1:4,1:4]
df_knn <- data.frame(SC = as.numeric(as.character(sc_04$seurat_clusters)))

df_knn$SC <- factor(df_knn$SC, levels = as.character(0:max(df_knn$SC)))
ind_order <- order(df_knn$SC, decreasing = F)
df_knn <- data.frame(SC=df_knn[ind_order,])
rownames(df_knn) <- colnames(sc_04)[ind_order]
head(df_knn)

Heatmap(heat[,match(rownames(df_knn), colnames(heat))],
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        
        clustering_method_rows = "ward.D",
        
        top_annotation = HeatmapAnnotation(df=df_knn,
                                           col = list(SC=c("0" = c25[1],"1" = c25[2],"2" = c25[3],"3" = c25[4],"4" = c25[5],"5" = c25[6],
                                                           "6" = c25[7],"8" = c25[8],"9" = c25[9],"10" = c25[10],"11" = c25[11],
                                                           "12" = c25[12]))),
        col = colorRamp2(breaks = c(-2,0,2),
                         colors = c("#1E90FF", "#FFF8DC", "#FF7F50")),
        name = "Expression",
        row_names_gp = gpar(fontsize=10),
        use_raster = F)

Heatmap(as.data.frame.matrix(table(sc_04$Manual_compiled, sc_04$seurat_clusters)),
        cluster_rows = F,
        cluster_columns = F,
        col=colorRamp2(breaks = c(0,50,400),
                       colors = c("white", "lightgray", "black")),
        name = "NumCells",
        border = T)


# add classification for AHH04 based on previous plots ####
sc_04$Manual_compiled_v2 <- ifelse(sc_04$seurat_clusters == 0 & sc_04@assays$Protein@scale.data["ADT-CD278-TotalSeqC",] > 2 &
                                     sc_04@assays$Protein@scale.data["ADT-CD352-TotalSeqC",] > 2 &
                                     sc_04@assays$SCT@data["GZMB",] > 2, "Tregs_Effector", "Tregs")
sc_04$Manual_compiled_v2 <- ifelse(sc_04$seurat_clusters == 1, "CD4_Effector_like", 
                                   ifelse(sc_04$seurat_clusters == 2, "CD8_Exhausted", 
                                          ifelse(sc_04$seurat_clusters == 3, "CD4_Naive", 
                                                 ifelse(sc_04$seurat_clusters == 4, "NK", 
                                                        ifelse(sc_04$seurat_clusters == 5, "CD8_Effector_like", 
                                                               ifelse(sc_04$seurat_clusters == 6, "Mixed", sc_04$Manual_compiled_v2))))))
table(sc_04$Manual_compiled_v2)
table(sc_04$Manual_compiled_v2, 
      sc_04$seurat_clusters)
ggplot(NULL, aes(x=names(table(sc_04$Manual_compiled_v2)),
                 y=round((as.numeric(table(sc_04$Manual_compiled_v2)) / sum(as.numeric(table(sc_04$Manual_compiled_v2))))*100,2))) +
  geom_bar(stat="identity", position = "dodge", fill="lightgray", col="black")+
  theme_classic() +
  xlab("CellType") + ylab("PercentageCells") +
  theme(axis.text.x = element_text(angle = 90,
                                   size=12)) +
  coord_flip()

# address Allart's questions 12.21 ####
rownames(sc_data_no_int)[grep("TRBV5", rownames(sc_data_no_int), ignore.case = T)]
VlnPlot(sc_data_no_int,
        slot = "data", 
        group.by = "seurat_clusters",
        log = T,
        features = c("TRBV6-1", "TRBV6-2", "TRBV6-4", "TRBV6-5", "TRBV6-6", "TRBV6-7", "TRBV6-8",
                     "TRBV5-1", "TRBV5-3" ,"TRBV5-4", "TRBV5-5", "TRBV5-6", "TRBV5-7"),
        stack = T)


df_knn <- data.frame(SC = sc_data_no_int$KNN,
                     AHH = sc_data_no_int$AHH)

df_knn <- df_knn[order(df_knn$SC, decreasing = F),]
df_knn$SC <- factor(df_knn$SC, levels = as.character(1:12))
head(df_knn)
table(df_knn$SC)

heat <- data.frame(sc_data_no_int@assays$SCT@data) %>%
  subset(rownames(.) %in% rownames(sc_data_no_int)[grep("TRBV", rownames(sc_data_no_int), ignore.case = T)])
colnames(heat) <- gsub("\\.","-",colnames(heat))
heat <- heat[,which(colnames(heat) %in% rownames(df_knn))]

heat <- heat[,-which(apply(heat,2,sum) == 0)]

cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
heat[1:10,1:4]
df_knn <- df_knn[which(rownames(df_knn) %in% colnames(heat)),]
all(rownames(df_knn) %in% colnames(heat))

for(ii in 1:12){
  temp_heat <- heat[,which(colnames(heat) %in% rownames(df_knn)[which(df_knn$SC == ii)])]
  temp_heat <- temp_heat[,order(colSums(temp_heat))]
  if(ii == 1){
    heat2 <- temp_heat
  } else {
    heat2 <- data.frame(cbind(heat2, temp_heat))
  }
}
all(rownames(df_knn) %in% colnames(heat))

Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        clustering_distance_rows = "pearson",
        top_annotation = HeatmapAnnotation(df = df_knn,
                                           col = list(SC=c("0" = c25[1],"1" = c25[2],"2" = c25[3],"3" = c25[4],"4" = c25[5],"5" = c25[6],
                                                           "6" = c25[7],"8" = c25[8],"9" = c25[9],"10" = c25[10],"11" = c25[11],
                                                           "12" = c25[12],"7"=c25[20]))),
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c("#AFEEEE", "white", "#B22222")),
        row_names_gp = gpar(fontsize=10),
        border = T,
        name="Expression")

# cd4 and CD8 split in each cluster - email question ####
FeatureScatter(sc_data_no_int,
               feature1 = "ADT-CD8-TotalSeqC",
               feature2 = "ADT-CD4-TotalSeqC",
               group.by = "KNN",
               cols = c25,
               plot.cor = F) +
  facet_wrap(~ sc_data_no_int$KNN)
ggplot(NULL, aes(x=as.numeric(sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",]),
                 y=as.numeric(sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",]))) +
  geom_point() +
  theme_classic() +
  labs(col="Cluster") + xlab("CD8") + ylab("CD4") +
  facet_wrap(~ sc_data_no_int$KNN)


# address Allart's questions 12.28 ####
trbv_genes <- rownames(sc_04)[grep("^TRBV", rownames(sc_04), ignore.case = T)]
ind_random <- sample(1:ncol(sc_04@assays$SCT), 3000, replace = F)
heat <- as.data.frame.matrix(sc_04@assays$SCT[which(rownames(sc_04@assays$SCT) %in% trbv_genes), ind_random])
heat <- heat[-which(apply(heat,1,sum) < 1), ]
heat <- heat[, -which(apply(heat,2,sum) < 1)]
rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls

# heat <- heat[-which(colVars(heat) < 0.1), ]
heat[1:10,1:5]
ind_heat <- which(rownames(sc_04@meta.data) %in% colnames(heat))
df_heat <- data.frame(sc_04@meta.data)[ind_heat,c("Manual_compiled_v2","seurat_clusters")]
df_heat <- df_heat[order(as.numeric(df_heat$seurat_clusters), decreasing = F),]
order_rows <- order(rowSums(heat), decreasing = T)
order_cols <- order(colSums(heat), decreasing = T)
head(df_heat)
Heatmap(heat[order_rows,match(rownames(df_heat), colnames(heat))],
        show_column_names = F,
        cluster_columns = F,
        cluster_rows = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list("seurat_clusters"=c("0"=c25[1],"1"=c25[2],"2"=c25[3],"3"=c25[4],"4"=c25[5],"5"=c25[6], "6"=c25[7]))),
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c("lightblue", "lightblue", "red")),
        name=paste("Scaled", "Expression", sep="\n"))


VlnPlot(sc_04, 
        slot = "data",
        features = ,
        group.by = "seurat_clusters",
        stack = T)

rownames(sc_04)[grep("^TRBV6", rownames(sc_04), ignore.case = T)]
ggplot(NULL, aes(x=as.numeric(sc_04@assays$SCT@data["TRBV6-1",]),
                 y=as.numeric(sc_04@assays$SCT@data["CD8A",]))) +
  geom_point() +
  theme_classic() +
  labs(col="Cluster") + 
  facet_wrap(~ sc_04$seurat_clusters, scales = "free")


# add VDJ data ####
library(djvdj)
library(scRepertoire)

paths <- c(
  Pt4 = "VDJ_4",
  Pt5 = "VDJ_5",
  Pt64 = "VDJ_64",
  Pt89 = "VDJ_89",
  Pt9 = "VDJ_9",
  Pt93 = "VDJ_93"
)

table(sc_04$PatientID.x)

for(ii in 1:length(paths)){
  
  pt_num <- strsplit(paths[ii], "_")[[1]][2]
  temp_add <- paste0("Pt", pt_num)
  
  temp_df <- read.table(paste0(paths[ii],"/filtered_contig_annotations.csv"),header = T, sep = ",")
  head(temp_df)
  
  temp_treatment <- paste0(pt_num, "-Categories-HashtagAssignment.csv")
  temp_treatment <- read.table(temp_treatment, sep=",", header = T)
  temp_treatment <- temp_treatment[-which(temp_treatment$AHH %in% c("Doublet", "Negative")),]
  head(temp_treatment)
  
  temp_df <- merge(temp_df, temp_treatment, by="barcode")
  
  
  if(ii == 1){
    vdj_data <- list(temp_df)
  } else {
    vdj_data[[ii]] <- temp_df
  }
  names(vdj_data)[ii] <- temp_add
}
vdj_data[[1]] %>%
  
  filter(grepl("TRBV6",v_gene)) %>%
  View

# extract beta chain data only ####
vdj_data_TRB <- rbindlist(lapply(vdj_data, function(xx) xx %>%
                                   filter(productive == "true" & chain == "TRB")),
                          idcol = T)
colnames(vdj_data_TRB)[1] <- "PatientID"
vdj_data_TRB$TRBV6 <- ifelse(grepl("TRBV6",vdj_data_TRB$v_gene), "YES", "NO")
vdj_data_TRB$TRBV6.1 <- ifelse(grepl("TRBV6-1",vdj_data_TRB$v_gene), "YES", "NO")
vdj_data_TRB$TRBV6.2 <- ifelse(grepl("TRBV6-2",vdj_data_TRB$v_gene), "YES", "NO")
vdj_data_TRB$TRBV6.5 <- ifelse(grepl("TRBV6-5",vdj_data_TRB$v_gene), "YES", "NO")
vdj_data_TRB$TRBV10.3 <- ifelse(grepl("TRBV10-3",vdj_data_TRB$v_gene), "YES", "NO")
vdj_data_TRB$TRBV20 <- ifelse(grepl("TRBV20",vdj_data_TRB$v_gene), "YES", "NO")

vdj_data_TRB$PtID_barcode <- paste0(vdj_data_TRB$barcode, "_", vdj_data_TRB$PatientID)

table(vdj_data_TRB$PatientID,
      vdj_data_TRB$AHH)
table(vdj_data_TRB$TRBV20,
      vdj_data_TRB$TRBV6.5)

# annotate main objects for TRVB usage ####
sc_data_no_int@meta.data[1:10,1:4]
sc_data_no_int$TRBV6 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV6 == "YES")], "YES", "NO")
sc_data_no_int$TRBV6.1 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV6.1 == "YES")], "YES", "NO")
sc_data_no_int$TRBV6.2 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV6.2 == "YES")], "YES", "NO")
sc_data_no_int$TRBV6.5 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV6.5 == "YES")], "YES", "NO")
sc_data_no_int$TRBV10.3 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV10.3 == "YES")], "YES", "NO")
sc_data_no_int$TRBV6 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV6 == "YES")], "YES", "NO")
sc_data_no_int$TRBV20 <- ifelse(sc_data_no_int@meta.data$barcode_ptID %in% vdj_data_TRB$PtID_barcode[which(vdj_data_TRB$TRBV20 == "YES")], "YES", "NO")
table(sc_data_no_int$AHH,
      sc_data_no_int$TRBV20)

df_TRBV_cells_AHH <- data.frame(AHH = c("AHH01", "AHH03", "AHH04"),
                                TRBV6 = as.data.frame.matrix(table(sc_data_no_int$AHH,
                                                                   sc_data_no_int$TRBV6))$YES,
                                TRBV6.1 = as.data.frame.matrix(table(sc_data_no_int$AHH,
                                                                     sc_data_no_int$TRBV6.1))$YES,
                                TRBV6.5 = as.data.frame.matrix(table(sc_data_no_int$AHH,
                                                                     sc_data_no_int$TRBV6.5))$YES,
                                TRBV20 = as.data.frame.matrix(table(sc_data_no_int$AHH,
                                                                    sc_data_no_int$TRBV20))$YES)
df_TRBV_cells_AHH <- melt(df_TRBV_cells_AHH)
ggplot(df_TRBV_cells_AHH, 
       aes(x = AHH,
           y=value,
           fill=variable)) +
  geom_bar(stat="identity", position = "dodge", col="black", width = 0.75) +
  theme_classic() +
  labs(fill = "Chain") + ylab("Number of cells")

deg_trbv6_pos_ned_cd8 <- subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which()])

sc_data_no_int@meta.data %>%
  filter(TRBV6 == "YES" & TRBV20 == "YES") %>%
  View

DimPlot(sc_data_no_int, 
        cells.highlight = colnames(sc_data_no_int)[which(sc_data_no_int$TRBV6.5 == "YES")], 
        reduction = "prot_UMAP",
        split.by = "AHH",
        pt.size = 0.1)
VlnPlot(sc_data_no_int, 
        group.by = "TRBV20",
        features = c("ADT-CD278-TotalSeqC","ADT-CD25-TotalSeqC"),
        split.by = "AHH",
        pt.size = 0)

# label CD8 and CD4 positive cells for VDJ analysis ####
VlnPlot(sc_data_no_int,
        assay = "Protein",
        group.by = "seurat_clusters",
        features = c("ADT-CD8-TotalSeqC", "ADT-CD4-TotalSeqC"),
        stack = T)

ggplot(data.frame(sc_data_no_int@meta.data),
       aes(x=sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",],
           y=sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",])) +
  geom_point(size=0.25, alpha=0.5) +
  theme_bw() +
  xlab("CD8 protein") + ylab("CD4 protein") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3))) +
  geom_vline(xintercept = 1.75) +
  geom_hline(yintercept = 1.5)

))
sc_data_no_int$CD8_CD4 <- ifelse(sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] > 1.75 & sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] < 1.5, "CD8",
                                 ifelse(sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] < 1.75 & sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.5, "CD4", 
                                        ifelse(sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",] > 1.75 & sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.5, "DP", "DN")))

ggplot(data.frame(sc_data_no_int@meta.data),
       aes(x=sc_data_no_int@assays$Protein@data["ADT-CD8-TotalSeqC",],
           y=sc_data_no_int@assays$Protein@data["ADT-CD4-TotalSeqC",],
           col=CD8_CD4)) +
  geom_point(size=0.25, alpha=0.5) +
  theme_bw() +
  xlab("CD8 protein") + ylab("CD4 protein") +
  facet_wrap(~ Protein_snn_res.0.25, ncol=5) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))


colnames(sc_data_no_int@meta.data)
sc_data_no_int$TRBV_used.v3 <- apply(sc_data_no_int@meta.data[,c("TRBV6.1", "TRBV6.2", "TRBV6.5","TRBV20", "TRBV10.3")], 1, function(xx) {
  
  
  xx <- as.character(xx)
  if(xx[1] == "YES" & sum(xx == "YES") == 1){
    chain <- "TRBV6.1"
  } else  if(xx[2] == "YES" & sum(xx == "YES") == 1){
    chain <- "TRBV6.2"
  } else  if(xx[3] == "YES" & sum(xx == "YES") == 1){
    chain <- "TRBV6.5"
  } else  if(xx[5] == "YES" & sum(xx == "YES") == 1){
    chain <- "TRBV10.3"
  }
  
  # test for multiple chains
  if(sum(xx == "YES") > 1){
    chain <- "mixed"
  }
  
  # test for TRBV20
  if(xx[1] == "NO" & xx[2] == "NO" & xx[3] == "NO" & xx[5] == "NO" & xx[4] == "YES"){
    chain <- "TRBV20"
  }
  
  # test for all negative
  if(xx[1] == "NO" & xx[2] == "NO" & xx[3] == "NO" & xx[5] == "NO" & xx[4] == "NO"){
    chain <- "other"
  }
  
  
  return(chain)
})

table(sc_data_no_int$TRBV_used.v3)
# TRBV6.1 vs. TRBV20 in AHH04
Idents(sc_data_no_int) <- sc_data_no_int$TRBV_used
res_dge_cd8_6.1_20_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 == "CD8"  & sc_data_no_int$AHH == "AHH04")]), 
                                     # slot = "data",
                                     min.cells.group = 3,
                                     ident.1 = "TRBV6.1", ident.2 = "TRBV20",
                                     logfc.threshold = 0.1, 
                                     min.pct = 0)
head(res_dge_cd8_6.1_20_04)
res_dge_cd4_6.1_20_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 == "CD4"  & sc_data_no_int$AHH == "AHH04")]), 
                                     # slot = "data",
                                     min.cells.group = 3,
                                     ident.1 = "TRBV6.1", ident.2 = "TRBV20",
                                     logfc.threshold = 0.1, 
                                     min.pct = 0)
head(res_dge_cd4_6.1_20_04)
res_dge_cd8cd4_6.1_20_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & sc_data_no_int$AHH == "AHH04")]), 
                                        # slot = "data",
                                        min.cells.group = 3,
                                        ident.1 = "TRBV6.1", ident.2 = "TRBV20",
                                        logfc.threshold = 0.1, 
                                        min.pct = 0)
head(res_dge_cd8cd4_6.1_20_04)

# TRBV6.5 vs. TRBV20 in AHH04
res_dge_cd8_6.5_20_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 == "CD8"  & sc_data_no_int$AHH == "AHH04")]), 
                                     # slot = "data",
                                     min.cells.group = 3,
                                     ident.1 = "TRBV6.5", ident.2 = "TRBV20",
                                     logfc.threshold = 0.1, 
                                     min.pct = 0)
head(res_dge_cd8_6.5_20_04)
res_dge_cd4_6.5_20_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 == "CD4"  & sc_data_no_int$AHH == "AHH04")]), 
                                     # slot = "data",
                                     min.cells.group = 3,
                                     ident.1 = "TRBV6.5", ident.2 = "TRBV20",
                                     logfc.threshold = 0.1, 
                                     min.pct = 0)
head(res_dge_cd4_6.5_20_04)
res_dge_cd8cd4_6.5_20_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & sc_data_no_int$AHH == "AHH04")]), 
                                        # slot = "data",
                                        min.cells.group = 3,
                                        ident.1 = "TRBV6.5", ident.2 = "TRBV20",
                                        logfc.threshold = 0.1, 
                                        min.pct = 0)
head(res_dge_cd8cd4_6.5_20_04)


# TRBV6.1 vs. TRBV6.5 in AHH04
res_dge_cd8_6.1vs6.5_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & sc_data_no_int$AHH == "AHH04")]), 
                                       # slot = "data",
                                       min.cells.group = 3,
                                       ident.1 = "TRBV6.1", ident.2 = "TRBV6.5",
                                       logfc.threshold = 0.1, 
                                       min.pct = 0)
head(res_dge_cd8_6.1vs6.5_04)
res_dge_cd4_6.1vs6.5_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & sc_data_no_int$AHH == "AHH04")]), 
                                       # slot = "data",
                                       min.cells.group = 3,
                                       ident.1 = "TRBV6.1", ident.2 = "TRBV6.5",
                                       logfc.threshold = 0.1, 
                                       min.pct = 0)
head(res_dge_cd4_6.1vs6.5_04)
res_dge_cd8cd4_6.1vs6.5_04 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & sc_data_no_int$AHH == "AHH04")]), 
                                          # slot = "data",
                                          min.cells.group = 3,
                                          ident.1 = "TRBV6.1", ident.2 = "TRBV6.5",
                                          logfc.threshold = 0.1, 
                                          min.pct = 0)
head(res_dge_cd8cd4_6.1vs6.5_04)


# TRBV6.1 AHH04 vs. AHH01
Idents(sc_data_no_int) <- sc_data_no_int$AHH
res_dge_cd8_6.1_04vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH04") &
                                                                                                    sc_data_no_int$TRBV6.1 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd8_6.1_04vs01)
res_dge_cd4_6.1_04vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH04") &
                                                                                                    sc_data_no_int$TRBV6.1 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd4_6.1_04vs01)
res_dge_cd8cd4_6.1_04vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & 
                                                                                                       sc_data_no_int$AHH %in% c("AHH01","AHH04") &
                                                                                                       sc_data_no_int$TRBV6.1 == "YES")]), 
                                         # slot = "data",
                                         min.cells.group = 3,
                                         ident.1 = "AHH04", ident.2 = "AHH01",
                                         logfc.threshold = 0.1, 
                                         min.pct = 0)
head(res_dge_cd8cd4_6.1_04vs01)


# TRBV6.1 AHH03 vs. AHH01
res_dge_cd8_6.1_03vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.1 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH03", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd8_6.1_03vs01)
write.table(res_dge_cd8_6.1_03vs01, "res_dge_cd8_6.1_03vs01.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd4_6.1_03vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.1 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH03", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd4_6.1_03vs01)
write.table(res_dge_cd4_6.1_03vs01, "res_dge_cd4_6.1_03vs01.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd8cd4_6.1_03vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & 
                                                                                                       sc_data_no_int$AHH %in% c("AHH01","AHH03") &
                                                                                                       sc_data_no_int$TRBV6.1 == "YES")]), 
                                         # slot = "data",
                                         min.cells.group = 3,
                                         ident.1 = "AHH03", ident.2 = "AHH01",
                                         logfc.threshold = 0.1, 
                                         min.pct = 0)
head(res_dge_cd8cd4_6.1_03vs01)
# TRBV6.1 AHH04 vs. AHH03
res_dge_cd8_6.1_04vs03 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH04","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.1 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH03",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd8_6.1_04vs03)
write.table(res_dge_cd8_6.1_04vs03, "res_dge_cd8_6.1_04vs03.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd4_6.1_04vs03 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH04","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.1 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH03",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)    
head(res_dge_cd4_6.1_04vs03)
write.table(res_dge_cd4_6.1_04vs03, "res_dge_cd4_6.1_04vs03.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd8cd4_6.1_04vs03 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & 
                                                                                                       sc_data_no_int$AHH %in% c("AHH04","AHH03") &
                                                                                                       sc_data_no_int$TRBV6.1 == "YES")]), 
                                         # slot = "data",
                                         min.cells.group = 3,
                                         ident.1 = "AHH04", ident.2 = "AHH03",
                                         logfc.threshold = 0.1, 
                                         min.pct = 0)    
head(res_dge_cd8cd4_6.1_04vs03)


# TRBV6.5 AHH04 vs. AHH01
res_dge_cd8_6.5_04vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH04") &
                                                                                                    sc_data_no_int$TRBV6.5 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd8_6.5_04vs01)
write.table(res_dge_cd8_6.5_04vs01, "res_dge_cd8_6.5_04vs01.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd4_6.5_04vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH04") &
                                                                                                    sc_data_no_int$TRBV6.5 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)    
head(res_dge_cd4_6.5_04vs01)
write.table(res_dge_cd4_6.5_04vs01, "res_dge_cd4_6.5_04vs01.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd8cd4_6.5_04vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & 
                                                                                                       sc_data_no_int$AHH %in% c("AHH01","AHH04") &
                                                                                                       sc_data_no_int$TRBV6.5 == "YES")]), 
                                         # slot = "data",
                                         min.cells.group = 3,
                                         ident.1 = "AHH04", ident.2 = "AHH01",
                                         logfc.threshold = 0.1, 
                                         min.pct = 0)    
head(res_dge_cd8cd4_6.5_04vs01)
write.table(res_dge_cd8cd4_6.5_04vs01, "res_dge_cd8cd4_6.5_04vs01.txt", sep="\t", quote = F, col.names = NA)
# TRBV6.5 AHH03 vs. AHH01
res_dge_cd8_6.5_03vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.5 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH03", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd8_6.5_03vs01)
write.table(res_dge_cd8_6.5_03vs01, "res_dge_cd8_6.5_03vs01.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd4_6.5_03vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH01","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.5 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH03", ident.2 = "AHH01",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)    
head(res_dge_cd4_6.5_03vs01)
write.table(res_dge_cd4_6.5_03vs01, "res_dge_cd4_6.5_03vs01.txt", sep="\t", quote = F, col.names = NA)
res_dge_cd8cd4_6.5_03vs01 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & 
                                                                                                       sc_data_no_int$AHH %in% c("AHH01","AHH03") &
                                                                                                       sc_data_no_int$TRBV6.5 == "YES")]), 
                                         # slot = "data",
                                         min.cells.group = 3,
                                         ident.1 = "AHH03", ident.2 = "AHH01",
                                         logfc.threshold = 0.1, 
                                         min.pct = 0)    


head(res_dge_cd8cd4_6.5_04vs03)

# TRBV6.5 AHH04 vs. AHH03
res_dge_cd8cd4_6.5_04vs03 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8","CD4")  & 
                                                                                                       sc_data_no_int$AHH %in% c("AHH04","AHH03") &
                                                                                                       sc_data_no_int$TRBV6.5 == "YES")]), 
                                         # slot = "data",
                                         min.cells.group = 3,
                                         ident.1 = "AHH04", ident.2 = "AHH03",
                                         logfc.threshold = 0.1, 
                                         min.pct = 0)    
head(res_dge_cd8cd4_6.5_04vs03)
res_dge_cd8_6.5_04vs03 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD8")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH04","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.5 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH03",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)
head(res_dge_cd8_6.5_04vs03)
res_dge_cd4_6.5_04vs03 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 %in% c("CD4")  & 
                                                                                                    sc_data_no_int$AHH %in% c("AHH04","AHH03") &
                                                                                                    sc_data_no_int$TRBV6.5 == "YES")]), 
                                      # slot = "data",
                                      min.cells.group = 3,
                                      ident.1 = "AHH04", ident.2 = "AHH03",
                                      logfc.threshold = 0.1, 
                                      min.pct = 0)    
head(res_dge_cd4_6.5_04vs03)

plot_violin <- function(results, name_plot){
  
  results$Significant <- ifelse(results$p_val_adj < 0.05, "YES", "NO")
  results$Direction <- ifelse(results$Significant == "YES" & results$avg_log2FC > 0, "UP",
                              ifelse(results$Significant == "YES" & results$avg_log2FC < 0, "DOWN", "NS"))
  results$size <- ifelse(results$Direction == "UP", 0.5,
                         ifelse(results$Direction == "DOWN", 0.5, 0.2))
  results$Labels <- ifelse(!grepl("^MT-",rownames(results)) & results$Direction != "NS", rownames(results), "")
  head(results)
  
  gg <- ggplot(results,
               aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
    geom_point() + 
    theme_pubclean() +
    scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                  "DOWN" = alpha("lightblue", 0.5),
                                  "NS" = "lightgray")) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_text_repel(label=results$Labels,
                    box.padding = 0.3, 
                    force = 0.15,
                    min.segment.length = 0,
                    max.overlaps = 20,
                    fontface = 'bold', 
                    color = 'black',
                    fill = "white") +
    guides(color=guide_legend(ncol=1, 
                              override.aes = list(size=3))) +
    ggtitle(name_plot)
  print(gg)
}

all_dge_results <- ls()[grep("res_dge_cd", ls())]
pdf("Violin_plots_dge_trbvs_ahhs.pdf", width = 5, height = 5)
for(ii in 1:length(all_dge_results)){
  plot_violin(get(all_dge_results[ii]),
              all_dge_results[ii])
  
  message(paste("done for", ii, all_dge_results[ii]))
}
dev.off()

plot_violin(res_dge_cd4_6.1_03vs01,
            deparse(substitute(res_dge_cd4_6.1_03vs01)))


DimPlot(sc_04,
        reduction = "UMAP_prot_o4",
        group.by = "PatientID.x")
VlnPlot(subset(sc_04, cells=colnames(sc_04)[which(sc_04$PatientID.x == "Pt64")]),
        features = "ADT-CD25-TotalSeqC",
        group.by = "Protein_snn_res.0.25",
        pt.size = 0.15)


# run GSEA for TRBV 27 DGE data ####
for(ii in 1:length(all_dge_results)){
  temp_res <- get(all_dge_results[ii])
  ranked_genes <- temp_res$avg_log2FC
  names(ranked_genes) <- rownames(temp_res)
  ranked_genes_ind <- order(ranked_genes, decreasing = T)
  ranked_genes <- ranked_genes[ranked_genes_ind]
  
  sig_temp <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.reactome.v7.2.symbols.gmt")
  fgsea_temp <- tryCatch(fgsea::fgsea(pathways = sig_temp,
                                      stats = ranked_genes),
                         error = function(xx){
                           message(xx)
                           message("\nadding NAs\n")
                           dummy_list <- NA
                           return(dummy_list)
                         })
  
  if(ii == 1){
    res_fsgea_reactome <- list(fgsea_temp)
    names(res_fsgea_reactome) <- all_dge_results[ii]
  } else {
    res_fsgea_reactome[[ii]] <- fgsea_temp
    names(res_fsgea_reactome)[ii] <- all_dge_results[ii]
  }
  
  message(paste("done for", ii, all_dge_results[ii]))
}
head(res_fsgea_reactome[[1]])
lapply(res_fsgea_reactome, function(xx) nrow(xx %>% filter(padj < 0.05)))
res_fsgea_reactome_sig <- rbindlist(lapply(res_fsgea_reactome, function(xx) xx %>% filter(padj < 0.05)), idcol = T)
colnames(res_fsgea_reactome_sig)[1] <- "Comparison"
res_fsgea_reactome_sig$leadingEdge <- as.character(res_fsgea_reactome_sig$leadingEdge)
res_fsgea_reactome_sig$Labels <- gsub("REACTOME_", "", res_fsgea_reactome_sig$pathway)
head(res_fsgea_reactome_sig)

table(res_fsgea_reactome_sig$Comparison)
ggplot(res_fsgea_reactome_sig %>% 
         filter(Comparison == "res_dge_cd4_6.5_04vs03") %>% 
         filter(NES > 0) %>%
         arrange(desc(NES)) %>%
         slice(seq_len(20)),
       aes(x = factor(Labels, levels = rev(as.character(Labels))),
           y=NES)) +
  geom_bar(stat="identity", position = "dodge", col="black", fill="#F0F8FF",width = 0.8) +
  theme_classic() +
  coord_flip() +
  ggtitle("res_dge_cd4_6.5_04vs03") +
  
  xlab("Pathways")

ggplot(as.data.frame.matrix(table(sc_04$PatientID.x,
                                  sc_04$Protein_snn_res.0.25)) %>%
         mutate(PatientID = rownames(.)) %>%
         melt(),
       aes(x=variable, y=value, fill=PatientID)) +
  geom_bar(stat="identity", position="dodge", width = 0.75, col="darkgray") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12)) +
  ylab("Number of cells") +
  scale_fill_manual(values = c25)

### continue from here !!!!!!!!!!!! ##########

sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c7.all.v7.4.symbols.gmt")
all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all,
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      })


# combine other VDJ AB with seurat object ####
combined <- combineTCR(lapply(vdj_data, function(xx) xx %>%
                                filter(chain == "TRB")), 
                       samples = rep("Pt",6),
                       removeNA = T,
                       removeMulti = T,
                       ID = gsub("Pt","",names(vdj_data)),
                       cells ="T-AB")


# remove trailing end -1
combined <- lapply(combined, function(xx){ 
  temp <- xx
  temp$barcode <- sapply(temp$barcode, function(yy) strsplit(yy, "-")[[1]][1])
  xx$barcode <- temp$barcode
  return(xx)
})

# add treatment information
for(ii in 1:length(combined)){
  
  pt_num <- gsub("Pt_","",names(combined)[ii])
  temp_add <- paste0("Pt_", pt_num, "_")
  
  temp_df <- combined[[ii]]
  
  temp_treatment <- paste0(pt_num, "-Categories-HashtagAssignment.csv")
  temp_treatment <- read.table(temp_treatment, sep=",", header = T)
  temp_treatment <- temp_treatment[-which(temp_treatment$AHH %in% c("Doublet", "Negative")),]
  temp_treatment$barcode <- gsub("-1","",paste0(temp_add, temp_treatment$barcode))
  head(temp_treatment)
  
  temp_df <- merge(temp_df, temp_treatment, by="barcode")

  combined[[ii]] <- temp_df
  
}


new_combined <- rbindlist(combined, use.names = T, idcol = "Pt_ID")
new_combined <- split(new_combined, f=new_combined$AHH)
lapply(new_combined, function(xx) table(xx$Pt_ID))

quantContig(new_combined, cloneCall="gene+nt", scale = TRUE)
quantContig_output <- quantContig(new_combined, cloneCall="gene+nt", 
                                  scale = TRUE, exportTable = TRUE)
kable(quantContig_output, format = "rst", align = "c")
quantContig(new_combined, cloneCall="aa", group = "AHH", scale = TRUE, exportTable = T)
abundanceContig(new_combined, cloneCall = "gene", scale = FALSE)

lengthContig(new_combined, cloneCall="aa", chains = "combined") 
lengthContig(new_combined, cloneCall="nt", chains = "single") 


clonalHomeostasis(new_combined, cloneCall = "aa")
clonalProportion(new_combined, cloneCall = "nt") 

clonalOverlap(new_combined, cloneCall = "gene+nt", method = "morisita")

clonalDiversity(new_combined, cloneCall = "aa", group = "samples", n.boots = 100)

compareClonotypes(new_combined, numbers = 20, samples = c("AHH01", "AHH03", "AHH04"), 
                  cloneCall="gene+nt", graph = "alluvial") +
  scale_fill_manual(values = c25)

clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")

sc_vdj <- sc_04

rownames(sc_vdj@meta.data) <- paste0("Pt_", gsub("Pt","",sc_vdj$PatientID),"_", rownames(sc_vdj@meta.data))
head(sc_vdj@meta.data)
rownames(sc_vdj@meta.data) <- sapply(rownames(sc_vdj@meta.data), function(xx) strsplit(xx, "-")[[1]][1])
table(sapply(rownames(sc_vdj@meta.data), function(xx) paste0(strsplit(xx, "_")[[1]][1],
                                                             strsplit(xx, "_")[[1]][2])))
head(rbindlist(new_combined), n=2)
table(sapply(rbindlist(combined)$barcode, function(xx) paste0(strsplit(xx, "_")[[1]][1],
                                                              strsplit(xx, "_")[[1]][2])))

sum(rownames(sc_vdj@meta.data) %in% combined[[1]]$barcode)
sc_vdj <- combineExpression(df = combined, 
                            sc = sc_vdj, 
                            cloneCall="aa", 
                            proportion = T)
head(sc_vdj)
table(sapply(rownames(sc_vdj@meta.data), function(xx) paste0(strsplit(xx, "_")[[1]][1],
                                                             strsplit(xx, "_")[[1]][2])))
sort(table(sc_vdj$CTaa), decreasing = T)[1:10]

if(all(all_of(sc_vdj$barcode_ptID %in% sc_data_no_int$barcode_ptID))){
  rownames(sc_vdj@meta.data) <- rownames(sc_04@meta.data)
 
}
sc_vdj$TRAV <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][1])),
  function(yy) strsplit(yy, "\\.")[[1]][1])
sc_vdj$TRAJ <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][1])),
  function(yy) strsplit(yy, "\\.")[[1]][2])
sc_vdj$TRAC <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][1])),
  function(yy) strsplit(yy, "\\.")[[1]][3])

sc_vdj$TRBV <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][3])),
  function(yy) strsplit(yy, "\\.")[[1]][1])
sc_vdj$TRBJ <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][3])),
  function(yy) strsplit(yy, "\\.")[[1]][2])
sc_vdj$TRBD <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][3])),
  function(yy) strsplit(yy, "\\.")[[1]][3])
sc_vdj$TRBC <- sapply(
  as.character(sapply(sc_vdj$CTstrict, function(xx) strsplit(xx, "_")[[1]][3])),
  function(yy) strsplit(yy, "\\.")[[1]][4])


sc_vdj@meta.data %>%
  View

colnames(sc_vdj@meta.data)
head(sc_vdj@meta.data)
DimPlot(sc_vdj %>%
          subset(cells=colnames(sc_vdj)[grepl("TRBV6", sc_vdj@meta.data$TRBV)]),
        reduction = "aUMAP2",
        cols = c25,
        group.by = "seurat_clusters") +
  xlab("UMAP1") + ylab("UMAP2")


sc_vdj$length_CDR3A <- nchar(sapply(sc_vdj$CTaa, function(xx) strsplit(xx, "_")[[1]][1]))
sc_vdj$length_CDR3B <- nchar(sapply(sc_vdj$CTaa, function(xx) strsplit(xx, "_")[[1]][2]))

VlnPlot(sc_vdj %>%
          subset(cells = colnames(sc_vdj)[-which(is.na(sc_vdj@meta.data$length_CDR3B))]),
        features = c("length_CDR3A", "length_CDR3B"),
        group.by = "Manual_compiled_v2")
FeatureScatter(sc_vdj %>%
                 subset(cells = colnames(sc_vdj)[-which(is.na(sc_vdj@meta.data$length_CDR3B))]),
               feature1 = "ADT-CD8-TotalSeqC",
               feature2 = "length_CDR3A",
               group.by = "Manual_compiled_v2")

Heatmap(as.data.frame.matrix(table(sc_vdj$TRBV[which(sc_vdj$AHH == "AHH04")], 
                                   sc_vdj$TRBJ[which(sc_vdj$AHH == "AHH04")])),
        col=colorRamp2(breaks = c(0,3,10,30),
                       colors = c("white","#FFF0F5","#6495ED", "#191970")),
        border = T,
        name = paste("Cell","Number",sep="\n"))


table(sc_vdj$length_CDR3A,
      sc_vdj$Manual_compiled_v2)


# label cells that use TRBV6 and compare to those who don't ####
rownames(sc_vdj@assays$Protein@data)
VlnPlot(sc_vdj,
        group.by = "seurat_clusters",
        features = c("ADT-CD8-TotalSeqC",
                     "ADT-CD4-TotalSeqC"),
        pt.size = 0,
        stack = T)

sc_vdj$TRBV6_use <- ifelse(grepl("TRBV6", sc_vdj$TRBV), "YES", "NO")
as.data.frame.matrix(table(sc_vdj$Protein_snn_res.0.25, sc_vdj$TRBV6_use)) %>%
  mutate(Cluster = rownames(.)) %>%
  mutate(NO = NULL) %>%
  melt(id.vars="Cluster") %>%
  ggplot(., aes(x=factor(Cluster, levels = as.character(as.numeric(0:6))),
                y = value)) +
  geom_bar(stat="identity", position="dodge", width=0.75, col="black", fill="orange") +
  theme_classic() +
  labs(fill=paste("TRBV6","recombination",sep="\n"))+
  xlab("Cluster") + ylab("Number of Cells")

ggplot(data.frame(sc_vdj@meta.data), 
       aes(x=sc_vdj@assays$Protein@data["ADT-CD8-TotalSeqC",],
           y=sc_vdj@assays$Protein@data["ADT-CD4-TotalSeqC",])) +
  geom_vline(xintercept = 1.75, col="red") +
  geom_hline(yintercept = 1.75, col="red") +
  geom_point(size=0.5, alpha=0.5) +
  theme_bw() +
  xlab("CD8 protein") + ylab("CD4 protein")

sc_vdj$CD8CD4 <- ifelse(sc_vdj@assays$Protein@data["ADT-CD8-TotalSeqC",] > 1.75 & sc_vdj@assays$Protein@data["ADT-CD4-TotalSeqC",] < 1.75, "CD8",
                        ifelse(sc_vdj@assays$Protein@data["ADT-CD8-TotalSeqC",] < 1.75 & sc_vdj@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.75, "CD4",
                               ifelse(sc_vdj@assays$Protein@data["ADT-CD8-TotalSeqC",] > 1.75 & sc_vdj@assays$Protein@data["ADT-CD4-TotalSeqC",] > 1.75, "DP_CD4CD8", "DN")))
table(sc_vdj$CD8CD4)

table(sc_vdj$Protein_snn_res.0.25, sc_vdj$TRBV6_use)



# run DGE for TRVB6 rearrangment ####

Idents(sc_vdj) <- sc_vdj$TRBV6_use
diff_trbv6_cd4 <- FindMarkers(subset(sc_vdj, cells=colnames(sc_vdj)[which(sc_vdj$CD8CD4 == "CD4"  & sc_vdj$AHH == "AHH04")]), 
                              
                              min.cells.group = 3,
                              ident.1 = "YES", ident.2 = "NO",
                              logfc.threshold = 0.1, 
                              min.pct = 0.1)
head(diff_trbv6_cd4)
diff_trbv6_cd4$Gene <- rownames(diff_trbv6_cd4)
diff_trbv6_cd4$Direction <- ifelse(diff_trbv6_cd4$avg_log2FC > 0 & diff_trbv6_cd4$p_val_adj < 0.05, "UP",
                                   ifelse(diff_trbv6_cd4$avg_log2FC < 0 & diff_trbv6_cd4$p_val_adj < 0.05, "DOWN", "NS"))
diff_trbv6_cd4$Labels <- ifelse(diff_trbv6_cd4$Direction != "NS",
                                
                                diff_trbv6_cd4$Gene,"")
diff_trbv6_cd4$Labels <- ifelse(diff_trbv6_cd4$Labels == "ACTB" | grepl("^MT-",diff_trbv6_cd4$Labels), "", diff_trbv6_cd4$Labels)
head(diff_trbv6_cd4)
ggplot(diff_trbv6_cd4, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_trbv6_cd4$Labels,
                  box.padding = 0.3, 
                  force = 0.5,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3))) +
  xlim(-1.75, 1.75)

diff_trbv6_cd8 <- FindMarkers(subset(sc_vdj, cells=colnames(sc_vdj)[which(sc_vdj$CD8CD4 == "CD8" & sc_vdj$AHH == "AHH04")]), 
                              
                              min.cells.group = 3,
                              ident.1 = "YES", ident.2 = "NO",
                              logfc.threshold = 0.1, 
                              min.pct = 0.1)
head(diff_trbv6_cd8)
diff_trbv6_cd8$Gene <- rownames(diff_trbv6_cd8)
diff_trbv6_cd8$Direction <- ifelse(diff_trbv6_cd8$avg_log2FC > 0 & diff_trbv6_cd8$p_val_adj < 0.05, "UP",
                                   ifelse(diff_trbv6_cd8$avg_log2FC < 0 & diff_trbv6_cd8$p_val_adj < 0.05, "DOWN", "NS"))
diff_trbv6_cd8$Labels <- ifelse(diff_trbv6_cd8$Direction != "NS",
                                
                                diff_trbv6_cd8$Gene,"")
diff_trbv6_cd8$Labels <- ifelse(diff_trbv6_cd8$Labels == "ACTB" | grepl("^MT-",diff_trbv6_cd8$Labels), "", diff_trbv6_cd8$Labels)
head(diff_trbv6_cd8)
ggplot(diff_trbv6_cd8, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_trbv6_cd8$Labels,
                  box.padding = 0.3, 
                  force = 0.5,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3))) +
  xlim(-1.75, 1.75)


# run GSEA for DGE above for TRBV6 ####
ranked_genes <- diff_trbv6_cd4$avg_log2FC
ranked_genes_ind <- order(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[ranked_genes_ind]
names(ranked_genes) <- rownames(diff_trbv6_cd4)[ranked_genes_ind]

sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c7.all.v7.4.symbols.gmt")
all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all,
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      }) %>%
  filter(padj <= 0.05)
View(all_fgsea[unique(c(grep("_TCELL", all_fgsea$pathway, ignore.case = T),
                        grep("_T_CELL", all_fgsea$pathway, ignore.case = T),
                        grep("_LYMPHOCYTE", all_fgsea$pathway, ignore.case = T))),] %>%
       filter(grepl("_UP", .$pathway, ignore.case = T)) %>%
       filter(grepl("CD4", .$pathway, ignore.case = T)))


fgsea_trbv6_cd8_hall <- all_fgsea
fgsea_trbv6_cd8_hall <- fgsea_trbv6_cd8_hall[order(fgsea_trbv6_cd8_hall$NES, decreasing = T),] %>%
  mutate(Labels = gsub("HALLMARK_", "", pathway))
fgsea_trbv6_cd8_reactome <- all_fgsea
fgsea_trbv6_cd8_reactome <- fgsea_trbv6_cd8_reactome[order(fgsea_trbv6_cd8_reactome$NES, decreasing = T),] %>%
  mutate(Labels = gsub("REACTOME_", "", pathway))

fgsea_trbv6_cd4_hall <- all_fgsea
fgsea_trbv6_cd4_hall <- fgsea_trbv6_cd4_hall[order(fgsea_trbv6_cd4_hall$NES, decreasing = T),] %>%
  mutate(Labels = gsub("HALLMARK_", "", pathway))
fgsea_trbv6_cd4_reactome <- all_fgsea
fgsea_trbv6_cd4_reactome <- fgsea_trbv6_cd4_reactome[order(fgsea_trbv6_cd4_reactome$NES, decreasing = T),] %>%
  mutate(Labels = gsub("REACTOME_", "", pathway))

ggplot(fgsea_trbv6_cd4_hall,
       aes(x = factor(Labels, levels = rev(as.character(Labels))),
           y=NES)) +
  geom_bar(stat="identity", position = "dodge", col="black", fill="#F0F8FF",width = 0.8) +
  theme_classic() +
  coord_flip() +
  geom_text_repel(label=c(fgsea_trbv6_cd4_hall$Labels),
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways")

# check intersections
ggVennDiagram::ggVennDiagram(x = list(CD8_UP = diff_trbv6_cd8$Gene[which(diff_trbv6_cd8$Direction == "UP")],
                                      CD4_UP = diff_trbv6_cd4$Gene[which(diff_trbv6_cd4$Direction == "UP")]))
ggVennDiagram::ggVennDiagram(x = list(CD8_DOWN = diff_trbv6_cd8$Gene[which(diff_trbv6_cd8$Direction == "DOWN")],
                                      CD4_DOWN = diff_trbv6_cd4$Gene[which(diff_trbv6_cd4$Direction == "DOWN")]))
degs_trbv6_common_up <- intersect(diff_trbv6_cd8$Gene[which(diff_trbv6_cd8$Direction == "UP")],
                                  diff_trbv6_cd4$Gene[which(diff_trbv6_cd4$Direction == "UP")])
degs_trbv6_common_down <- intersect(diff_trbv6_cd8$Gene[which(diff_trbv6_cd8$Direction == "DOWN")],
                                    diff_trbv6_cd4$Gene[which(diff_trbv6_cd4$Direction == "DOWN")])

# check TF in TRBV6 ####
degs_TF_trbv6_common_up <- degs_trbv6_common_up[which(degs_trbv6_common_up %in% all_TFs)]
degs_TF_trbv6_common_down <- degs_trbv6_common_down[which(degs_trbv6_common_down %in% all_TFs)]

# plot common TF logFC
TF_common_for_plot <- diff_trbv6_cd8 %>%
  filter(p_val_adj < 0.05) %>%
  
  filter(Gene %in% degs_TF_trbv6_common_down) %>%
  mutate(CD8CD4 = "CD8") %>%
  arrange(desc(avg_log2FC))
order_genes <- TF_common_for_plot$Gene
TF_common_for_plot <- diff_trbv6_cd4 %>%
  filter(p_val_adj < 0.05) %>%
  
  filter(Gene %in% degs_TF_trbv6_common_down) %>%
  mutate(CD8CD4 = "CD4") %>%
  rbind(TF_common_for_plot) %>%
  arrange(desc(avg_log2FC))
ggplot(TF_common_for_plot, 
       aes(x=factor(Gene, levels = rev(order_genes)),
           y=avg_log2FC,
           fill=CD8CD4)) +
  geom_bar(stat="identity", position = "dodge", col="gray", width = 0.9) +
  theme_classic() +
  
  scale_fill_manual(values = c("CD4" = "#1E90FF", "CD8" = "#000080")) +
  coord_flip() +
  xlab("Gene Name") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=10))


# run DGE analysis analysis C2 and C5 in sc_04 ####
# run CD8 C5 vs. C2
Idents(sc_04) <- sc_04$seurat_clusters
diff_sc04_c2_c5 <- FindMarkers(subset(sc_04, cells=colnames(sc_04)[which(sc_04$seurat_clusters %in% c(2,5))]), 
                               min.cells.group = 3,
                               ident.1 = "5", ident.2 = "2",
                               logfc.threshold = 0.1, 
                               min.pct = 0.1)
head(diff_sc04_c2_c5)
diff_sc04_c2_c5$Gene <- rownames(diff_sc04_c2_c5)
diff_sc04_c2_c5$Direction <- ifelse(diff_sc04_c2_c5$avg_log2FC > 0 & diff_sc04_c2_c5$p_val_adj < 0.05, "UP",
                                    ifelse(diff_sc04_c2_c5$avg_log2FC < 0 & diff_sc04_c2_c5$p_val_adj < 0.05, "DOWN", "NS"))
diff_sc04_c2_c5$Labels <- ifelse(diff_sc04_c2_c5$Direction != "NS" &
                                   diff_trbv6_cd4$pct.1 > 0.4 ,
                                 diff_sc04_c2_c5$Gene,"")
diff_sc04_c2_c5$Labels <- ifelse(diff_sc04_c2_c5$Labels == "ACTB" | grepl("^MT-",diff_sc04_c2_c5$Labels), "", diff_sc04_c2_c5$Labels)
head(diff_sc04_c2_c5)
ggplot(diff_sc04_c2_c5, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_sc04_c2_c5$Labels,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
# check for TRBV6
diff_sc04_c2_c5 %>% 
  filter(grepl("TRBV6", Gene)) %>%
  kable(format = "rst", align = "c")
write.table(diff_sc04_c2_c5, "DEG_12312021_AHH04_C5vsC2.txt", sep="\t", quote = F)

# run CD4 C vs. C1
# Idents(sc_04) <- sc_04$seurat_clusters
diff_sc04_c0_c1 <- FindMarkers(subset(sc_04, cells=colnames(sc_04)[which(sc_04$seurat_clusters %in% c(0,1))]), 
                               # slot = "data",
                               min.cells.group = 3,
                               ident.1 = "1", ident.2 = "0",
                               logfc.threshold = 0.1, 
                               min.pct = 0.1)
head(diff_sc04_c0_c1)
diff_sc04_c0_c1$Gene <- rownames(diff_sc04_c0_c1)
diff_sc04_c0_c1$Direction <- ifelse(diff_sc04_c0_c1$avg_log2FC > 0 & diff_sc04_c0_c1$p_val_adj < 0.05, "UP",
                                    ifelse(diff_sc04_c0_c1$avg_log2FC < 0 & diff_sc04_c0_c1$p_val_adj < 0.05, "DOWN", "NS"))
diff_sc04_c0_c1$Labels <- ifelse(diff_sc04_c0_c1$Direction != "NS",
                                 diff_sc04_c0_c1$Gene,"")
diff_sc04_c0_c1$Labels <- ifelse(diff_sc04_c0_c1$Labels == "ACTB" | grepl("^MT-",diff_sc04_c0_c1$Labels), "", diff_sc04_c0_c1$Labels)
head(diff_sc04_c0_c1)
ggplot(diff_sc04_c0_c1, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_sc04_c0_c1$Labels,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
# check for TRBV6
diff_sc04_c0_c1 %>% 
  filter(grepl("TRBV6", Gene)) %>%
  kable(format = "rst", align = "c")
write.table(diff_sc04_c0_c1, "DEG_12312021_AHH04_CD4_C1vsC0.txt", sep="\t", quote = F)

table(sc_04$Manual_compiled_v2, sc_04$seurat_clusters)

# re-cluster AHH04 based on rna data ####
sc_04 <-  RunPCA(sc_04, reduction.name = 'pca_SCT_04') %>%
  RunUMAP(features = VariableFeatures(sc_04), reduction.name = 'UMAP_SCT_o4', reduction = "pca_SCT_04") %>%
  FindNeighbors(reduction = "pca_SCT_04", dims = 1:20) %>% 
  FindClusters(resolution = 0.25)

DimPlot(sc_04, 
        reduction = "UMAP_SCT_o4",
        cols = c25,
        group.by = "seurat_clusters") +
  xlab("UMAP1_RNA") + ylab("UMAP2_RNA")


Heatmap(as.data.frame.matrix(table(sc_04$Protein_snn_res.0.25, sc_04$seurat_clusters)),
        cluster_rows = F,
        cluster_columns = F,
        col=colorRamp2(breaks = c(0,100,300),
                       colors = c("white", "#7FFFD4", "#0000CD")),
        name = "NumCells", 
        row_title = "Protein Clusters", column_title = "RNA Clusters",
        border = T)

FeaturePlot(sc_04, 
            reduction = "UMAP_SCT_o4",
            
            features = c("ADT-CD8-TotalSeqC"),
            min.cutoff = "q50") +
  xlab("UMAP1_RNA") + ylab("UMAP2_RNA")


# check contribution of CD278 and CD25 to protein PCA for clusters 0,1,2,5 ####
pca_04_loadings <- RunPCA(sc_04, assay = "Protein", reduction.name = 'PCA_temp_protein')
pca_04_loadings <- data.frame(pca_04_loadings[["PCA_temp_protein"]]@feature.loadings)[,1:5] %>%
  arrange(desc(abs(PC_1)))
pca_04_loadings <- data.frame(Feature = rownames(pca_04_loadings),
                              Contribution = abs(pca_04_loadings$PC_1))
tt_ncells_for_pca <- table(sc_04$Protein_snn_res.0.25)[c(1,2,3,6)]
total_cells <- 4754
pca_04_loadings <- pca_04_loadings %>%
  mutate(C0 = (Contribution * tt_ncells_for_pca[[1]]) / total_cells) %>%
  mutate(C1 = (Contribution * tt_ncells_for_pca[[2]]) / total_cells) %>%
  mutate(C2 = (Contribution * tt_ncells_for_pca[[3]]) / total_cells) %>%
  mutate(C5 = (Contribution * tt_ncells_for_pca[[4]]) / total_cells)



head(pca_04_loadings)
ggplot(pca_04_loadings %>%
         filter(Feature %in% c("ADT-CD25-TotalSeqC", "ADT-CD278-TotalSeqC")),
       aes(x=variable, y=value, fill=Feature)) +
  geom_bar(stat="identity", position="dodge", width = 0.75, col="black") +
  theme_classic() +
  scale_fill_manual(values = c("ADT-CD25-TotalSeqC" = "#696969", 
                               "ADT-CD278-TotalSeqC" = "#F5F5F5")) +
  ylab("Contribution Score")

ggplot(pca_04_loadings %>%
         mutate(Feature = gsub("ADT-","",gsub("-TotalSeqC","", Feature))),
       aes(x=factor(Feature, levels = unique(as.character(Feature))), y=value, fill=Feature)) +
  geom_bar(stat="identity", position="dodge", width = 0.75, col="black") +
  theme_classic() +
  scale_fill_manual(values = c("ADT-CD25-TotalSeqC" = "#696969", 
                               "ADT-CD278-TotalSeqC" = "#F5F5F5")) +
  ylab("Contribution Score") +
  facet_wrap(~ variable, scales = "free_x") +
  theme(axis.text.x = element_text(angle=90, size=12))



ind_clusters <- which(as.numeric(as.character(sc_04$Protein_snn_res.0.25)) %in% c(0,1,2,5))
ggplot(data.frame(sc_04@meta.data[ind_clusters,]),
       aes(x=sc_04@assays$SCT@scale.data["IL2RA", ind_clusters],
           y=sc_04@assays$Protein@data["ADT-CD25-TotalSeqC",ind_clusters])) +
  geom_point(size=0.5, alpha=0.5) +
  theme_bw() +
  xlab("CD25/IL2RA RNA") + ylab("CD25/IL2RA Protein") +
  facet_wrap(~ Protein_snn_res.0.25, scales = "free") +
  stat_cor(method = "spearman",
           aes(label = ..r.label..), 
           label.x=5, label.y=1)

ggplot(data.frame(sc_04@meta.data[ind_clusters,]),
       aes(x=sc_04@assays$SCT@scale.data["ICOS", ind_clusters],
           y=sc_04@assays$Protein@data["ADT-CD278-TotalSeqC",ind_clusters])) +
  geom_point(size=0.5, alpha=0.5) + 
  theme_bw() +
  xlab("CD278/ICOS RNA") + ylab("CD278/ICOS Protein") +
  facet_wrap(~ Protein_snn_res.0.25, scales = "free") +
  stat_cor(method = "spearman",
           aes(label = ..r.label..), 
           label.x=5, label.y=1)

diff_sc04_c0_c1 %>%
  filter(Gene %in% c("ICOS","IL2RA"))
diff_sc04_c2_c5 %>%
  filter(Gene %in% c("ICOS","IL2RA"))

score_CD25_C0C1 <- -log2((0.09045160 * 0.08344835) * (2^-7.147514e-79) * (0.49 + 0.37))
score_CD278_C0C1 <- -log2((0.07024597 * 0.0648071) * (2^-1) * (0.22 + 0.1))
score_CD25_C2C5 <- -log2((0.06631992 * 0.04649145) * (2^-1.100577e-94) * (0.57 + 0.37))
score_CD278_C2C5 <- -log2((0.05150498 * 0.03610591) * (2^-1) * (0.12 + 0.11))

df_scores <- data.frame(Scores = c(score_CD25_C0C1,
                                   score_CD278_C0C1,
                                   score_CD25_C2C5,
                                   score_CD278_C2C5),
                        Gene = rep(c("CD25/IL2RA","CD278/ICOS"), 2),
                        Clusters = rep(c("C0_C1","C2_C5"), each=2))
ggplot(df_scores, aes(x=Clusters, y=Scores, fill=Gene)) +
  geom_bar(stat="identity", position="dodge", col="black") + 
  theme_classic() +
  scale_fill_manual(values = c("CD25/IL2RA" = "#F5F5F5","CD278/ICOS" = "#808080"))

pca_04_loadings <- melt(pca_04_loadings, id.vars=c("Feature", "Contribution"))

## otehre studd ####
DimPlot(subset(sc_vdj, cells=colnames(sc_vdj)[which(sc_vdj$TRBV6_use == "YES")]), 
        reduction = "aUMAP2",
        group.by = "CD8CD4",
        split.by = "CD8CD4")

genes2check <- c("GZMB", "IFNG", "CD69", "ZAP70", "LCK", "IL2RA","KLF2", "PRF1", "IRF3", "BATF3", "IRF4", "MAF", "JUNB", "EOMES", "TCF7")
VlnPlot(subset(sc_vdj, cells=colnames(sc_vdj)[which(sc_vdj$TRBV6_use == "YES")]),
        features = genes2check) +
  facet_wrap(~ factor(sc_vdj$KNN))

ggplot(data.frame(sc_vdj@assays$SCT@data)[genes2check,] %>%
         t() %>%
         data.frame(cbind(sc_vdj@meta.data[,c("TRBV6_use", "TILPRED_compiled", "TIL.cluster")])) %>%
         melt(id.vars = c("TRBV6_use", "TILPRED_compiled", "TIL.cluster")), 
       aes(x=TIL.cluster,
           y=log2(value),
           fill=TRBV6_use)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90, size=11)) +
  facet_wrap(~ variable, scales = "free_y", ncol=5)


View(sc_vdj@meta.data)
table(sc_vdj@meta.data$Manual_compiled_v2[which(grepl("Large", sc_vdj$cloneType))] %>%
        subset(grepl("TRBV6", .$CTgene)))

rownames(sc_04@assays$Protein@data)
rownames(sc_04@assays$SCT@data)[grep("cd69", rownames(sc_04@assays$SCT@data), ignore.case = T)]

FeatureScatter(subset(sc_vdj, cells=colnames(sc_vdj)[which(grepl("TRBV6", sc_vdj$CTgene))]),
               feature1 = "ADT-CD4-TotalSeqC",
               feature2 = "ADT-CD8-TotalSeqC",
               group.by = "Manual_compiled_v2") +
  facet_wrap(~ sc_04$Manual_compiled_v2)

VlnPlot(sc_vdj,
        features = c("GZMB", "IFNG", "CD69", "GZMA", "IL2RA", "PRF1"),
        group.by = "CTaa",
        pt.size = 0) 



occupiedscRepertoire(sc_vdj, x.axis = "cluster")

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
head(sc_vdj@meta.data)

table(sc_04$Protein_snn_res.0.25)
sc_vdj$KMeans <- sc_vdj$KNN

alluvialClonotypes(subset(sc_vdj, cells=colnames(sc_vdj)[which(sc_vdj$CD8CD4 %in% c("CD8", "CD4"))]),
                   cloneCall = "aa",
                   y.axes = c("PatientID", "Protein_snn_res.0.25", "CD8CD4"),
                   color = "TRBV6_use") +
  geom_stratum() +
  geom_text(stat = ggalluvial::StatStratum, 
            infer.label = F, reverse = TRUE, size = 4, colour="black") +
  theme(axis.text.x = element_text(angle=90, size=12),
        axis.text.y = element_text(size=10))



# grossly divide CD4 vs. CD8 clusters, run DGE with TRBV groups ###

sc_vdj$TRBV6_use
alluvialClonotypes(sc_vdj , 
                   cloneCall = "aa", 
                   y.axes = c("PatientID", "KMeans", "AHH"),
                   color = ifelse(grepl("TRBV6",sc_vdj$TRBV), "red", "lightgray")) +
  
  circles <- getCirclize(sc_vdj, cloneCall = "aa", groupBy = "seurat_clusters")

#Just assigning the normal colors to each cluster
Idents(sc_vdj) <- sc_vdj$seurat_clusters
grid.cols <- hue_pal()(length(unique(Idents(sc_vdj))))
names(grid.cols) <- levels(slot(sc_vdj, "active.ident"))


#Graphing the chord diagram
chordDiagram(circles, self.link = 2, grid.col = grid.cols)


# make various QC plots of all cells ####
head(sc_data_no_int)
kable(as.data.frame.matrix(table(sc_data_no_int$KNN,
                                 sc_data_no_int$PatientID)) %>%
        mutate(ClusterID = rownames(.)), format = "rst", align = "c")
VlnPlot(sc_data_no_int,
        group.by = "KNN",
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        split.by = "PatientID",
        pt.size = 0, ncol=1)
apply(sc_data_no_int@meta.data[,c("nCount_RNA", "nFeature_RNA", "nCount_Protein")], 2, median)

ggplot(as.data.frame.matrix(table(sc_data_no_int$PatientID,
                                  sc_data_no_int$AHH)) %>%
         mutate(PatientID = rownames(.)) %>%
         melt(),
       aes(x=variable, y=value, fill=PatientID)) +
  geom_bar(stat="identity", position="dodge", col="black") +
  theme_classic() +
  xlab("Treatment Group") +
  ylab("num cells") +
  scale_fill_manual(values = c25)


DimPlot(sc_vdj)

kable(as.data.frame.matrix(table(sc_04$Protein_snn_res.0.25,
                                 sc_04$PatientID)) %>%
        mutate(ClusterID = rownames(.)), format = "rst", align = "c")


# run GSEA for C0-C1 and C2-C5 ####
head(diff_sc04_c0_c1)
head(diff_sc04_c2_c5)
VlnPlot(sc_04, assay = "SCT", features = "LGALS1", group.by = "Protein_snn_res.0.25")

sc_04
# run GSEA for DGE above for TRBV6 ####
ranked_genes <- diff_sc04_c0_c1$avg_log2FC
ranked_genes_ind <- order(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[ranked_genes_ind]
names(ranked_genes) <- rownames(diff_sc04_c0_c1)[ranked_genes_ind]

sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c7.all.v7.2.symbols.gmt")
all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all,
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      }) %>%
  filter(padj <= 0.05)
head(all_fgsea)
all_fgsea %>%
  filter(padj < 0.05 & NES > 0) %>%
  filter(grepl("_UP", pathway)) %>%
  filter(grepl("_CD4_", pathway)) %>%
  View

# up NES = high in C1 (Effector) vs C0 (Exhausted)
# up NES = high in C5 (Effector) vs C2 (Exhausted)

sign_all_reactome_c0c1 <- all_fgsea
sign_all_reactome_c0c1 <- sign_all_reactome_c0c1 %>%
  arrange(desc(NES)) %>%
  mutate(Labels = gsub("REACTOME_", "", pathway))


sign_all_immune_c0c1 <- all_fgsea %>%
  filter(padj < 0.05 & NES > 0) %>%
  filter(grepl("_UP", pathway)) %>%
  filter(grepl("_CD4_", pathway))
sign_all_immune_c0c1 <- sign_all_immune_c0c1[-unique(c(grep("_ACT_", sign_all_immune_c0c1$pathway, ignore.case = T),
                                                       grep("_[0-9]H_", sign_all_immune_c0c1$pathway, ignore.case = T),
                                                       grep("_[0-9][0-9]H_", sign_all_immune_c0c1$pathway, ignore.case = T),
                                                       grep("_DAY[0-9]_", sign_all_immune_c0c1$pathway, ignore.case = T)))]
sign_all_immune_c0c1 <- sign_all_immune_c0c1 %>%
  arrange(desc(NES)) %>%
  mutate(Labels = gsub("GSE[0-9]*_", "", pathway))

# make  plot with leading edge genes from the first few pathways

ind_paths <- 33:38
genes4heat <- unique(unlist(sign_all_reactome_c0c1$leadingEdge[ind_paths]))

ind_cells <- which(as.character(sc_04$Protein_snn_res.0.25) %in% c("0", "1"))
heat <- sc_04@assays$SCT@data[genes4heat,ind_cells]


rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
heat[1:10,1:10]
colnames(heat) <- cls
rownames(heat) <- rws
df_anno_cols <- data.frame(Clusters = sc_04$Protein_snn_res.0.25[ind_cells],
                           AHH = sc_04$AHH[ind_cells])

df_anno_rows <- data.frame(Pathway = rep(sign_all_reactome_c0c1$pathway[ind_paths], lapply(sign_all_reactome_c0c1$leadingEdge[ind_paths], length)),
                           Genes = unlist(sign_all_reactome_c0c1$leadingEdge[ind_paths]))
df_anno_rows <- df_anno_rows[-which(duplicated(df_anno_rows$Genes)),]
rownames(df_anno_rows) <- df_anno_rows$Genes

df_anno_rows$Pathway <- gsub("REACTOME_","",df_anno_rows$Pathway)
df_anno_rows$Pathway <- gsub("_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS_","",df_anno_rows$Pathway)
df_anno_rows$Pathway <- gsub("THE_CITRIC_ACID_","",df_anno_rows$Pathway)
df_anno_rows <- data.frame(df_anno_rows[-2])

heat <- heat[match(rownames(df_anno_rows), rownames(heat)), match(rownames(df_anno_cols),colnames(heat))]


Heatmap(heat,
        show_row_names = T,
        show_column_names = F,
        cluster_columns = T,
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c("#87CEEB", scales::alpha("snow",0.3), "#DC143C")),
        top_annotation = HeatmapAnnotation(df = df_anno_cols,
                                           col = list(Clusters = c("0" = "green",
                                                                   "1" = "gray"),
                                                      AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        left_annotation = rowAnnotation(df = df_anno_rows),
        use_raster = F)



npaths <- nrow(sign_all_immune_c0c1)
ggplot(sign_all_immune_c0c1 %>%
         arrange(desc(NES)) %>%
         slice(seq_len(npaths)),
       aes(x = factor(Labels, levels = rev(as.character(Labels))),
           y=NES)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           fill=ifelse(sign_all_immune_c0c1$NES > 0, "#FFA07A", "#E0FFFF"), 
           alpha=0.8,
           width = 0.8) +
  theme_classic() +
  coord_flip() +
  geom_text_repel(label=sign_all_immune_c0c1$Labels[1:npaths],
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways")



sign_all_reactome_c2c5 <- all_fgsea
sign_all_reactome_c2c5 <- sign_all_reactome_c2c5 %>%
  mutate(Labels = sign_all_reactome_c2c5$pathway) %>%
  arrange(desc(NES)) %>%
  mutate(Labels = gsub("REACTOME_","", Labels))

sign_all_immune_c2c5 <- all_fgsea %>%
  filter(padj < 0.05 & NES > 0) %>%
  filter(grepl("_UP", pathway)) %>%
  filter(grepl("_CD8_", pathway))
sign_all_immune_c2c5 <- sign_all_immune_c2c5[-unique(c(grep("_ACT_", sign_all_immune_c2c5$pathway, ignore.case = T),
                                                       grep("_[0-9]H_", sign_all_immune_c2c5$pathway, ignore.case = T),
                                                       grep("_[0-9][0-9]H_", sign_all_immune_c2c5$pathway, ignore.case = T),
                                                       grep("_DAY[0-9]_", sign_all_immune_c2c5$pathway, ignore.case = T)))]


sign_all_immune_c2c5 <- sign_all_immune_c2c5 %>%
  arrange(desc(NES)) %>%
  mutate(Labels = gsub("GSE[0-9]*_", "", pathway))

# make  plot with leading edge genes from the first few pathways
# ind_paths <- c(28,31,32,36:38)
ind_paths <- c(1:3, 9, 13)
genes4heat <- unique(unlist(sign_all_immune_c2c5$leadingEdge[ind_paths]))
# heat <- sample(1:ncol(sc_04), 5000, replace = F)
ind_cells <- which(as.character(sc_04$Protein_snn_res.0.25) %in% c("5", "2"))
heat <- sc_04@assays$SCT@data[genes4heat,ind_cells]


rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
heat[1:10,1:10]
colnames(heat) <- cls
rownames(heat) <- rws
df_anno_cols <- data.frame(Clusters = sc_04$Protein_snn_res.0.25[ind_cells],
                           AHH = sc_04$AHH[ind_cells])

df_anno_rows <- data.frame(Pathway = rep(sign_all_immune_c2c5$pathway[ind_paths], lapply(sign_all_immune_c2c5$leadingEdge[ind_paths], length)),
                           Genes = unlist(sign_all_immune_c2c5$leadingEdge[ind_paths]))
df_anno_rows <- df_anno_rows[-which(duplicated(df_anno_rows$Genes)),]
rownames(df_anno_rows) <- df_anno_rows$Genes
df_anno_rows$Pathway <- gsub("GSE[0-9]*_","",df_anno_rows$Pathway)

df_anno_rows <- data.frame(df_anno_rows[-2])

heat <- heat[match(rownames(df_anno_rows), rownames(heat)), match(rownames(df_anno_cols),colnames(heat))]


Heatmap(heat,
        show_row_names = T,
        show_column_names = F,
        cluster_columns = T,
        col = colorRamp2(breaks = c(-1,0,1),
                         colors = c("#87CEEB", scales::alpha("snow",0.3), "#DC143C")),
        top_annotation = HeatmapAnnotation(df = df_anno_cols,
                                           col = list(Clusters = c("5" = "orange",
                                                                   "2" = "lightgray"),
                                                      AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"))),
        left_annotation = rowAnnotation(df = df_anno_rows),
        row_names_gp = gpar(fontsize=10),
        use_raster = F)

npaths <- nrow(sign_all_immune_c2c5)
ggplot(sign_all_immune_c2c5 %>%
         slice(seq_len(npaths)),
       aes(x = factor(Labels, levels = rev(as.character(Labels))),
           y=NES)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           fill=ifelse(sign_all_immune_c2c5$NES > 0, "#FFA07A", "#E0FFFF"), 
           alpha=0.8,
           width = 0.8) +
  theme_classic() +
  coord_flip() +
  geom_text_repel(label=sign_all_immune_c2c5$Labels[1:npaths],
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways")


# >> addrss points from meeting with KCL ####
# check patients contribution to PCA
sc_data_no_int
ll <- data.frame(Loadings(sc_data_no_int[["apca2"]]))
ll <- ll[order(abs(ll$PC_1), decreasing = T),]
ll
DimPlot(sc_data_no_int, 
        reduction = "apca2",
        group.by = "PatientID", 
        cols = c25) +
  xlab("Protein PCA1") + ylab("Protein PCA2")

FeaturePlot(sc_data_no_int,
            reduction = "apca2",
            order = T,
            features = rownames(ll)[1:6],
            cols = c("#FFFACD", "#DC143C"),
            pt.size = 0.5)


sc_data_no_int@assays$Protein@data[1:10,1:4]
pp <- prcomp(sc_data_no_int@assays$Protein@data[-c(1:4),])
summary(pp)
pca_loadings <- data.frame(pp$rotation)
pca_loadings <- pca_loadings[,1:3] %>%
  mutate(ID = rownames(pca_loadings)) %>%
  mutate(AbsLoading = abs(PC1))
pca_loadings <- merge(pca_loadings, sc_data_no_int@meta.data[,c("AHH", "PatientID", "Protein_snn_res.0.25")], by="row.names")
pca_loadings <- pca_loadings[order(abs(pca_loadings$PC1), decreasing = T),]
head(pca_loadings)
ggplot(pca_loadings, 
       aes(x=PatientID,
           y=AbsLoading)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

df_x <- data.frame(pp$x[,1:2])
df_x$Labels <- gsub("ADT-","",gsub("-TotalSeqC","",rownames(df_x)))
df_x$CorCoef <- sapply(rownames(df_x), function(xx) cor(sc_data_no_int@assays$Protein@data[xx,],
                                                        pca_loadings$PC1))
ggplot(df_x, 
       aes(x=PC1,
           y=PC2,
           col=CorCoef)) +
  geom_point(size=4) +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) +
  geom_text_repel(label=df_x$Labels,
                  box.padding = 0.1, 
                  force = 2.75,
                  min.segment.length = 3.5,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red")

ggplot(df_x, 
       aes(x=PC1,
           y=CorCoef)) +
  geom_point(size=4) +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) +
  geom_text_repel(label=df_x$Labels,
                  box.padding = 0.1, 
                  force = 2.75,
                  min.segment.length = 3.5,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black')


VlnPlot(sc_data_no_int,
        assay = "Protein",
        features = c("ADT-CD8-TotalSeqC",
                     "ADT-CD4-TotalSeqC",
                     "ADT-CD278-TotalSeqC",
                     "ADT-CD25-TotalSeqC",
                     "ADT-CD163-TotalSeqC",
                     "ADT-CD196-TotalSeqC"),
        group.by = "PatientID",
        pt.size = 0.001)

DimPlot(sc_data_no_int, 
        reduction = "pca",
        group.by = "PatientID", 
        cols = c25) +
  xlab("RNA PCA1") + ylab("RNA PCA2")


# >> address Allat's email 01.25.2022 << ####
# MAIN Analysis 
# #1: Restricted to the AHH04 treatment group and to determine differences between non-responding (CD25neg) and responding (CD25pos) cells. 
# Perform analysis separately for CD4pos and CD8pos cells, 
# i.e. CD4pos CD25pos vs CD4pos CD25neg and CD8pos CD25pos versus CD8pos CD25neg

head(sc_04)
DimPlot(sc_04, reduction = "UMAP_prot_o4")

DefaultAssay(sc_04) <- "SCT"
FeaturePlot(sc_04,
            features = c("IL2RA", "ADT-CD25-TotalSeqC", "ADT-CD8-TotalSeqC","ADT-CD4-TotalSeqC"),
            reduction = "UMAP_prot_o4",
            ncol=2,
            min.cutoff = "q25",
            order = T)

# select CD4pos CD25pos vs CD4pos CD25neg
# CD4pos_CD25pos 
ggplot(sc_04@meta.data,
       aes(x = sc_04@assays$Protein@data["ADT-CD8-TotalSeqC",],
           y = sc_04@assays$Protein@data["ADT-CD25-TotalSeqC",],
           col=Protein_snn_res.0.25)) +
  geom_point(size=0.5) +
  theme_bw() +
  xlab("CD8") + ylab("CD25") +
  labs(col="Cluster") +
  scale_color_manual(values = c25) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5)))



# run DGE analysis analysis C0 and C1 == CD4-CD25 ####
Idents(sc_04) <- sc_04$seurat_clusters
diff_sc04_c0_c1 <- FindMarkers(subset(sc_04, cells=colnames(sc_04)[which(sc_04$seurat_clusters %in% c(0,1))]), 
                               assay = "SCT",
                               min.cells.group = 3,
                               ident.1 = "0", ident.2 = "1",
                               logfc.threshold = 0.25, 
                               min.pct = 0.25)
head(diff_sc04_c0_c1)
diff_sc04_c0_c1$Gene <- rownames(diff_sc04_c0_c1)
diff_sc04_c0_c1$Direction <- ifelse(diff_sc04_c0_c1$avg_log2FC > 0 & diff_sc04_c0_c1$p_val_adj < 0.05, "UP",
                                    ifelse(diff_sc04_c0_c1$avg_log2FC < 0 & diff_sc04_c0_c1$p_val_adj < 0.05, "DOWN", "NS"))
diff_sc04_c0_c1$Labels <- ifelse(diff_sc04_c0_c1$Direction != "NS",
                                 diff_sc04_c0_c1$Gene,"")
diff_sc04_c0_c1$Labels <- ifelse(diff_sc04_c0_c1$Labels == "ACTB" | grepl("^MT-",diff_sc04_c0_c1$Labels), "", diff_sc04_c0_c1$Labels)
head(diff_sc04_c0_c1)
ggplot(diff_sc04_c0_c1, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_sc04_c0_c1$Labels,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
write.table(diff_sc04_c0_c1, "DGE_CD4_CD25_PosVsNeg_01272022.txt", sep="\t", quote = F, col.names = NA)

# run DGE analysis analysis C4 and C4 == CD8-CD25 ####
Idents(sc_04) <- sc_04$seurat_clusters
diff_sc04_c4_c2 <- FindMarkers(subset(sc_04, cells=colnames(sc_04)[which(sc_04$seurat_clusters %in% c(2,4))]), 
                               assay = "SCT",
                               # slot = "data",
                               min.cells.group = 3,
                               ident.1 = "2", ident.2 = "4",
                               logfc.threshold = 0.25, 
                               min.pct = 0.25)
head(diff_sc04_c4_c2)
diff_sc04_c4_c2$Gene <- rownames(diff_sc04_c4_c2)
diff_sc04_c4_c2$Direction <- ifelse(diff_sc04_c4_c2$avg_log2FC > 0 & diff_sc04_c4_c2$p_val_adj < 0.05, "UP",
                                    ifelse(diff_sc04_c4_c2$avg_log2FC < 0 & diff_sc04_c4_c2$p_val_adj < 0.05, "DOWN", "NS"))
diff_sc04_c4_c2$Labels <- ifelse(diff_sc04_c4_c2$Direction != "NS",
                                 diff_sc04_c4_c2$Gene,"")
diff_sc04_c4_c2$Labels <- ifelse(diff_sc04_c4_c2$Labels == "ACTB" | grepl("^MT-",diff_sc04_c4_c2$Labels), "", diff_sc04_c4_c2$Labels)
head(diff_sc04_c4_c2)
ggplot(diff_sc04_c4_c2, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_sc04_c4_c2$Labels,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
write.table(diff_sc04_c4_c2, "DGE_CD8_CD25_PosVsNeg_01272022.txt", sep="\t", quote = F, col.names = NA)


# run fgsea analysis for the above DGE results ####
# run GSEA for C0-C1 and C2-C5 ####
head(diff_sc04_c0_c1)
head(diff_sc04_c2_c5)
VlnPlot(sc_04, assay = "SCT", features = "LGALS1", group.by = "Protein_snn_res.0.25")

sc_04


# run GSEA for DGE above ####
ranked_genes <- diff_sc04_c4_c2$avg_log2FC
ranked_genes_ind <- order(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[ranked_genes_ind]
names(ranked_genes) <- rownames(diff_sc04_c4_c2)[ranked_genes_ind]

sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c3.tft.v7.2.symbols.gmt")
all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all,
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      }) %>%
  filter(padj <= 0.05)
head(all_fgsea)
all_fgsea %>%
  filter(padj < 0.05 & NES > 0) %>%
  # filter(grepl("_UP", pathway)) %>%
  # filter(grepl("_CD4_", pathway)) %>%
  View

# save CD4 results
sign_all_reactome_c0c1 <- all_fgsea
sign_all_gobp_c0c1 <- all_fgsea
sign_all_hall_c0c1 <- all_fgsea
sign_all_immune_c0c1 <- all_fgsea %>%
  filter(padj < 0.05 & NES > 0) %>%
  filter(grepl("_UP", pathway)) %>%
  filter(grepl("_CD4_", pathway))
sign_all_immune_c0c1 <- sign_all_immune_c0c1[-unique(c(grep("_ACT_", sign_all_immune_c0c1$pathway, ignore.case = T),
                                                       grep("_[0-9]H_", sign_all_immune_c0c1$pathway, ignore.case = T),
                                                       grep("_[0-9][0-9]H_", sign_all_immune_c0c1$pathway, ignore.case = T),
                                                       grep("_DAY[0-9]_", sign_all_immune_c0c1$pathway, ignore.case = T)))]
sign_all_immune_c0c1 <- sign_all_immune_c0c1 %>%
  arrange(desc(NES)) %>%
  mutate(Labels = NULL)


# save CD8 results
sign_all_reactome_c2c4 <- all_fgsea
sign_all_gobp_c2c4 <- all_fgsea
sign_all_hall_c2c4 <- all_fgsea
sign_all_immune_c2c4 <- all_fgsea %>%
  filter(padj < 0.05 & NES > 0) %>%
  filter(grepl("_UP", pathway)) %>%
  filter(grepl("_CD8_", pathway))
sign_all_immune_c2c4 <- sign_all_immune_c2c4[-unique(c(grep("_ACT_", sign_all_immune_c2c4$pathway, ignore.case = T),
                                                       grep("_[0-9]H_", sign_all_immune_c2c4$pathway, ignore.case = T),
                                                       grep("_[0-9][0-9]H_", sign_all_immune_c2c4$pathway, ignore.case = T),
                                                       grep("_DAY[0-9]_", sign_all_immune_c2c4$pathway, ignore.case = T)))]
sign_all_immune_c2c4 <- sign_all_immune_c2c4 %>%
  arrange(desc(NES)) %>%
  mutate(pathway = gsub("GSE[0-9]*_", "", pathway))

# compile all results 
CD8_pathways_04 <- rbindlist(list(sign_all_reactome_c2c4, sign_all_gobp_c2c4, sign_all_hall_c2c4, sign_all_immune_c2c4), idcol = F) %>%
  mutate(CellType = "CD8") %>%
  mutate(leadingEdge = as.character(leadingEdge))

CD4_pathways_04 <- rbindlist(list(sign_all_reactome_c0c1, sign_all_gobp_c0c1, sign_all_hall_c0c1, sign_all_immune_c0c1), idcol = F) %>%
  mutate(CellType = "CD4") %>%
  mutate(leadingEdge = as.character(leadingEdge))

pathways_all_04 <- data.frame(rbind(CD8_pathways_04,
                                    CD4_pathways_04))
pathways_all_04 %>%
  filter(grepl("TCR", pathway)) %>%
  View

pathways_all_04 <- pathways_all_04 %>%
  arrange((NES))
write.table(pathways_all_04, "pathways_all_AHH04_CD4-CD8_CD25.txt", sep = "\t", row.names = F, quote = F)

npaths <- 40
ggplot(pathways_all_04 %>%
         slice(seq_len(npaths)),
       aes(x = factor(P_G, levels = rev(as.character(P_G))),
           y=NES,
           fill=CellType)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.8,
           width = 0.8) +
  theme_classic() +
  coord_flip() +
  geom_text_repel(label=pathways_all_04$pathway[1:npaths],
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways")


# cehck TF ####
diff_sc04_c0_c1 %>%
  filter(Gene %in% all_TF) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")

diff_sc04_c4_c2 %>%
  filter(Gene %in% all_TF) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")


# TCR signaling ####
rownames(sc_04@assays$RNA@counts)[grep("akt", rownames(sc_04@assays$RNA@counts), ignore.case = T)]
diff_sc04_c4_c2 %>%
  filter(Gene %in% c("LCK", "ZAP70", "CD28", "CD28", "AKT3", "AKT1", "AKT2", "LAT", "NFKBIA", "RELA", "RELB"))

diff_sc04_c0_c1 %>%
  filter(Gene %in% c("LCK", "ZAP70", "CD28", "CD28", "AKT3", "AKT1", "AKT2", "LAT", "NFKBIA", "RELA", "RELB"))


# Gene expression patterns for CD25 positive CD4 and CD8 ?? ####
rownames(sc_04@assays$RNA@counts)[grep("eomes", rownames(sc_04@assays$RNA@counts), ignore.case = T)]

diff_sc04_c0_c1 %>%
  filter(Gene %in% c("CCR7", "SELL", "IL7R", "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "BTLA", "CD160", "TGFB1", "TGFBR1", "CD69", "CX3CR1", "TBX21", "CD103", "CD44", "EOMES"))

diff_sc04_c4_c2 %>%
  filter(Gene %in% c("CCR7", "SELL", "IL7R", "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "BTLA", "CD160", "TGFB1", "TGFBR1", "CD69", "CX3CR1", "TBX21", "CD103", "CD44", "EOMES"))

# map slide 7 by patient + percentage ####
DimPlot(sc_data_no_int, 
        reduction = "prot_UMAP",
        group.by = "PatientID.x")
tt <- as.data.frame.matrix(table(sc_data_no_int@meta.data$PatientID.x,
                                 sc_data_no_int@meta.data$seurat_clusters))

# perc relative to all cells
data.frame(round(tt / sum(tt), 5)*100) %>% kable(format = "rst", align = "c")

# perc relative to each cluster
data.frame(round(apply(tt, 2, function(xx) xx/sum(xx)), 5)*100) %>% kable(format = "rst", align = "c")

DimPlot(sc_04, 
        reduction = "UMAP_prot_o4",
        group.by = "PatientID.x")
tt <- as.data.frame.matrix(table(sc_04@meta.data$PatientID.x,
                                 sc_04@meta.data$seurat_clusters))


# perc relative to all cells
data.frame(round(tt / sum(tt), 4)*100) %>% kable(format = "rst", align = "c")


# perc relative to each cluster
data.frame(round(apply(tt, 2, function(xx) xx/sum(xx)), 4)*100) %>% kable(format = "rst", align = "c")
rownames(sc_04@assays$Protein@data)

plot(sc_04@assays$Protein@data["ADT-CD3-TotalSeqC",],
     sc_04@assays$Protein@data["ADT-CD8-TotalSeqC",])

ggplot(data.frame(cbind(sc_04@meta.data,
                        CD8 = log10(sc_04@assays$Protein@counts["ADT-CD8-TotalSeqC",]),
                        IL2RA = log10(sc_04@assays$Protein@counts["ADT-CD25-TotalSeqC",]))),
       aes(x = CD8,
           y = IL2RA)) +
  geom_point(size=0.125) +
  theme_bw() +
  geom_vline(xintercept = 1.5, col="red") +
  geom_hline(yintercept = 1.6, col="red") +
  facet_wrap(~ PatientID.x) +
  stat_quadrant_counts(xintercept = 1.5, yintercept = 1.6)

# make table of patients and numbers:
CD25_table <- data.frame(CD8_CD25 = c(512, 55, 54, 12, 42, 202),
                         CD8_total = c(512+111, 55+70, 54+70, 12+17, 42+24, 202+76),
                         CD4_CD25 = c(640, 73, 132, 15, 90, 479),
                         CD4_total = c(640+288, 73+163, 132+118, 15+63, 90+92, 479+273))
CD25_table %>%
  mutate(perc_pos_over_neg = ((CD8_CD25 + CD4_CD25) / (CD8_total + CD4_total))*100) %>%
  kable(format = "rst", align = "c")


# Address remaining questions compiled in Allart's email from 01.25.2022 using notes from 01.27.2022 ####
extra_genes <- readxl::read_xlsx("20220120 List with genes of interest.xlsx",sheet = 1, skip = 1)


# ---------
# for the 10% just run UMAP plot CD25 vs. mitochondrial DNA to see if there is any association between CD25 expression, metabolic genes/oxphos and percentage of mitochondrial gene expression
# YES, if an association exists in certain donors where increasing the tolerance to a higher % of mitochondrial gene expression 
# also increases the number of CD25+, viable cells, then we should rerun the analysis including these additional cells. 
# The bottom-line is that we try to understand why we have so few cells for certain donors whereas these cells seem to have been there in earlier Loupe browser analyzed files.
# ---------


sc_data
head(sc_data)
max(sc_data$percent.mt)
table(sc_data$PatientID.x,
      sc_data$AHH)

ggplot(sc_data@meta.data[-which(sc_data$AHH %in% c("Doublet", "Negative")),],
       aes(x=PatientID.y,
           y=percent.mt,
           col=AHH)) +
  geom_violin(scale="width", alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12))

rownames(sc_data@assays$Protein@data)
ggplot(sc_data@meta.data[-which(sc_data$AHH %in% c("Doublet", "Negative")),],
       aes(x=PatientID.y,
           y=log10(sc_data@assays$Protein@data["ADT-CD25-TotalSeqC",-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           col=AHH)) +
  ylab("Log10(CD25 Protein)") +
  geom_violin(scale="width", alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12))

ggplot(sc_data@meta.data[-which(sc_data$AHH %in% c("Doublet", "Negative")),],
       aes(x=(sc_data@meta.data$percent.mt[-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           y=(sc_data@assays$Protein@data["ADT-CD25-TotalSeqC",-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           col=AHH)) +
  ylab("Log10(CD25 Protein)") +
  xlab("Percentage.Mt") +
  geom_point(size=0.125, alpha=0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12)) +
  geom_vline(aes(xintercept = 10), col="orange", linetype="dashed") +
  geom_vline(aes(xintercept = 5), col="purple", linetype="dashed") +
  facet_wrap(sc_data$PatientID.y[-which(sc_data$AHH %in% c("Doublet", "Negative"))] ~ sc_data$AHH[-which(sc_data$AHH %in% c("Doublet", "Negative"))],
             ncol=3)

ggplot(sc_data@meta.data[-which(sc_data$AHH %in% c("Doublet", "Negative")),],
       aes(x=(sc_data@meta.data$percent.mt[-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           y=log10(sc_data@assays$Protein@data["ADT-CD25-TotalSeqC",-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           col=AHH)) +
  ylab("Log10(CD25 Protein)") +
  xlab("Percentage.Mt") +
  geom_point(size=0.125, alpha=0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12)) +
  geom_vline(aes(xintercept = 10), col="orange", linetype="dashed") +
  geom_vline(aes(xintercept = 5), col="purple", linetype="dashed") +
  facet_wrap(sc_data$PatientID.y[-which(sc_data$AHH %in% c("Doublet", "Negative"))] ~ sc_data$AHH[-which(sc_data$AHH %in% c("Doublet", "Negative"))],
             ncol=3) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=5)))

sc_data$percent.mt_5 <- ifelse(sc_data$percent.mt < 5, 1, 0)
sc_data$percent.mt_5_10 <- ifelse(sc_data$percent.mt < 10 & sc_data$percent.mt_5 == 0, "5 < x < 10",
                                  ifelse(sc_data$percent.mt_5 == 1, "< 5", "> 10"))

ggplot(sc_data@meta.data[-which(sc_data$AHH %in% c("Doublet", "Negative")),],
       aes(x=(sc_data@meta.data$nFeature_RNA[-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           y=log10(sc_data@assays$Protein@data["ADT-CD25-TotalSeqC",-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
           col=percent.mt_5_10)) +
  ylab("Log10(CD25 Protein)") +
  xlab("nFeature_RNA") +
  geom_point(size=0.25, alpha=0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12)) +
  geom_vline(aes(xintercept = 300), col="brown", linetype="dashed") +
  geom_vline(aes(xintercept = 4000), col="brown", linetype="dashed") +
  facet_wrap(sc_data$PatientID.y[-which(sc_data$AHH %in% c("Doublet", "Negative"))] ~ sc_data$AHH[-which(sc_data$AHH %in% c("Doublet", "Negative"))],
             ncol=3) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=5)))

as.data.frame.matrix(table(paste0(sc_data@meta.data$AHH[-which(sc_data$AHH %in% c("Doublet", "Negative"))],
                                  "_",
                                  sc_data@meta.data$PatientID.x[-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
                           sc_data@meta.data$percent.mt_5_10[-which(sc_data$AHH %in% c("Doublet", "Negative"))])) %>%
  kbl(caption = "") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  row_spec(1:6, bold = F, color = "black", background = "#EEE8AA") %>%
  row_spec(7:12, bold = F, color = "black", background = "#E0FFFF") %>%
  row_spec(13:18, bold = F, color = "black", background = "#FFE4C4")

sc_data$nFeature_filtered <- ifelse(sc_data$nFeature_RNA > 300 & sc_data$nFeature_RNA < 4000, "300 < x < 4000", "S or D")

as.data.frame.matrix(table(paste0(sc_data@meta.data$AHH[-which(sc_data$AHH %in% c("Doublet", "Negative"))],
                                  "_",
                                  sc_data@meta.data$PatientID.x[-which(sc_data$AHH %in% c("Doublet", "Negative"))]),
                           paste0(sc_data@meta.data$percent.mt_5_10[-which(sc_data$AHH %in% c("Doublet", "Negative"))],
                                  "_",
                                  sc_data@meta.data$nFeature_filtered[-which(sc_data$AHH %in% c("Doublet", "Negative"))]))) %>%
  kbl(caption = "") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  row_spec(1:6, bold = F, color = "black", background = "#EEE8AA") %>%
  row_spec(7:12, bold = F, color = "black", background = "#E0FFFF") %>%
  row_spec(13:18, bold = F, color = "black", background = "#FFE4C4") %>%
  column_spec(column = 2, bold = T, color = "green", border_left = T, border_right = T, background = "#FFFAFA") %>%
  column_spec(column = 6, bold = T, color = "orange", border_left = T, border_right = T, background = "#FFFAFA") %>%
  row_spec(6, hline_after = T)


# -------
# Comparing differences in gene expression and pathways induced between AHH04 and AHH03 treatment. 
# First analysis to be top level and therefore use the Seurat clusters which included all treatment groups. 
# For CD4: compare cluster 0 (AHH03) vs Cluster 5 (AHH04); 
# for CD8: compare cluster 4 (AHH03) vs cluster 7 (AHH04)
# -------

sc_data_no_int
DimPlot(sc_data_no_int, reduction = "prot_UMAP")
FeaturePlot(sc_data_no_int, 
            reduction = "prot_UMAP", 
            features = "ADT-CD25-TotalSeqC", 
            order=T, 
            min.cutoff = "q50")


# run DGE analysis analysis For CD4: compare cluster 0 (AHH03) vs Cluster 5 (AHH04) ####
Idents(sc_data_no_int) <- sc_data_no_int$seurat_clusters
diff_04vs03_cd4_c5_c0 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$seurat_clusters %in% c(0,5))]), 
                                     assay = "SCT",
                                     # slot = "data",
                                     min.cells.group = 3,
                                     ident.1 = "5", ident.2 = "0",
                                     logfc.threshold = 0.25, 
                                     min.pct = 0.25)
head(diff_04vs03_cd4_c5_c0)
diff_04vs03_cd4_c5_c0$Gene <- rownames(diff_04vs03_cd4_c5_c0)
diff_04vs03_cd4_c5_c0$Direction <- ifelse(diff_04vs03_cd4_c5_c0$avg_log2FC > 0 & diff_04vs03_cd4_c5_c0$p_val_adj < 0.05, "UP",
                                          ifelse(diff_04vs03_cd4_c5_c0$avg_log2FC < 0 & diff_04vs03_cd4_c5_c0$p_val_adj < 0.05, "DOWN", "NS"))

diff_04vs03_cd4_c5_c0$Labels <- ifelse(diff_04vs03_cd4_c5_c0$Labels == "ACTB" | grepl("^MT-",diff_04vs03_cd4_c5_c0$Labels), "", diff_04vs03_cd4_c5_c0$Labels)
head(diff_04vs03_cd4_c5_c0)
ggplot(diff_04vs03_cd4_c5_c0, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_04vs03_cd4_c5_c0$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
diff_04vs03_cd4_c5_c0$InList <- ifelse(diff_04vs03_cd4_c5_c0$Gene %in% extra_genes$FeatureName, "YES","NO")

diff_04vs03_cd4_c5_c0 %>%
  filter(Gene == "RUNX3")

table(diff_04vs03_cd4_c5_c0$InList)
table(diff_04vs03_cd4_c5_c0$InList,
      diff_04vs03_cd4_c5_c0$Direction)


# for CD8: compare cluster 4 (AHH03) vs cluster 7 (AHH04)
diff_04vs03_cd8_c7_c4 <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$seurat_clusters %in% c(7,4))]), 
                                     assay = "SCT",
                                     min.cells.group = 3,
                                     ident.1 = "7", ident.2 = "4",
                                     logfc.threshold = 0.25, 
                                     min.pct = 0.25)
head(diff_04vs03_cd8_c7_c4)
diff_04vs03_cd8_c7_c4$Gene <- rownames(diff_04vs03_cd8_c7_c4)
diff_04vs03_cd8_c7_c4$Direction <- ifelse(diff_04vs03_cd8_c7_c4$avg_log2FC > 0 & diff_04vs03_cd8_c7_c4$p_val_adj < 0.05, "UP",
                                          ifelse(diff_04vs03_cd8_c7_c4$avg_log2FC < 0 & diff_04vs03_cd8_c7_c4$p_val_adj < 0.05, "DOWN", "NS"))

diff_04vs03_cd8_c7_c4$Labels <- ifelse(diff_04vs03_cd8_c7_c4$Labels == "ACTB" | grepl("^MT-",diff_04vs03_cd8_c7_c4$Labels), "", diff_04vs03_cd8_c7_c4$Labels)
head(diff_04vs03_cd8_c7_c4)

diff_04vs03_cd8_c7_c4 %>%
  filter(Gene == "RUNX3")

ggplot(diff_04vs03_cd8_c7_c4, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_04vs03_cd8_c7_c4$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
diff_04vs03_cd8_c7_c4$InList <- ifelse(diff_04vs03_cd8_c7_c4$Gene %in% extra_genes$FeatureName, "YES","NO")
table(diff_04vs03_cd8_c7_c4$InList)
table(diff_04vs03_cd8_c7_c4$InList,
      diff_04vs03_cd8_c7_c4$Direction)


# run GSEA for DGE above ####
ranked_genes <- diff_04vs03_cd8_c7_c4$avg_log2FC
ranked_genes_ind <- order(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[ranked_genes_ind]
names(ranked_genes) <- rownames(diff_04vs03_cd8_c7_c4)[ranked_genes_ind]

sig_all <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.reactome.v7.2.symbols.gmt")
all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all,
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      }) %>%
  filter(padj <= 0.05)
head(all_fgsea)
all_fgsea %>%
  filter(padj < 0.25 & NES > 0) %>%

  View

all_fgsea <- all_fgsea[order(all_fgsea$NES, decreasing = T),] %>%
  filter(padj <= 0.05) %>%
  mutate(P_G = gsub("REACTOME_", "", pathway))
npaths <- nrow(all_fgsea)
ggplot(all_fgsea %>%
         # arrange(desc(NES)),
         slice(seq_len(npaths)),
       aes(x = factor(P_G, levels = rev(as.character(P_G))),
           # aes(x = pathway,
           y=NES)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", fill="lightgray",
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  coord_flip() +
  geom_text_repel(label=all_fgsea$P_G[1:npaths],
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways")

# do leading edge analysis ####
all_fgsea_cd4 <- all_fgsea
all_fgsea_cd4$leadingEdge <- as.character(all_fgsea_cd4$leadingEdge)
all_fgsea_cd8 <- all_fgsea
all_fgsea_cd8$leadingEdge <- as.character(all_fgsea_cd8$leadingEdge)
write.table(all_fgsea_cd4, "CD4_pathways_02032022.txt", sep='\t', quote = F, row.names = F)
write.table(all_fgsea_cd8, "CD8_pathways_02032022.txt", sep='\t', quote = F, row.names = F)

cd4_lead_edge <- unique(unlist(all_fgsea_cd4$leadingEdge))
cd8_lead_edge <- unique(unlist(all_fgsea_cd8$leadingEdge))


cd8_lead_edge[-which(cd8_lead_edge %in% cd4_lead_edge)]

length(intersect(cd8_lead_edge, cd4_lead_edge))
ggVennDiagram::ggVennDiagram(x = list(CD8 = cd8_lead_edge, 
                                      CD4 = cd4_lead_edge))
all_leading_edge <- unique(c(cd8_lead_edge, cd4_lead_edge))

all_leading_edge <- unique(c(unique(unlist(all_fgsea_04vs03_cd8_c7_c4$leadingEdge)), 
                             unique(unlist(all_fgsea_04vs03_cd4_c5_c0$leadingEdge))))  # this is from the analysis on 02.06.2022




# make heatmap for leading edge #### 
ind_rows <- which(rownames(sc_data_no_int@assays$SCT@scale.data) %in% all_leading_edge)
ind_cols <- which(colnames(sc_data_no_int) %in% rownames(sc_data_no_int@meta.data)[
  which(as.character(sc_data_no_int$seurat_clusters) %in% c("0", "5", "4", "7"))
])
heat <- sc_data_no_int@assays$SCT@scale.data[ind_rows,
                                             ind_cols]

ncol(heat)
rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
rownames(heat) <- rws
heat[90:100,1:10]
df_anno_cols <- data.frame(Clusters = sc_data_no_int$seurat_clusters[ind_cols],
                           AHH = sc_data_no_int$AHH[ind_cols]) %>%
  mutate(Clusters= factor(Clusters, levels = as.character(c(0,4,5,7))))

table(df_anno_cols$Clusters, df_anno_cols$AHH)


for(ii in unique(as.numeric(as.character(df_anno_cols$Clusters)))){
  temp_clust <- ii
  temp_df <- df_anno_cols[which(df_anno_cols$Clusters == ii),]
  temp_df <- temp_df[order(temp_df$AHH),]
  if(ii == 0){
    df_anno_cols_2 <- temp_df
    print(ii)
  } else {
    df_anno_cols_2 <- data.frame(rbind(df_anno_cols_2,
                                       temp_df))
    print(ii)
  }
}
head(df_anno_cols_2)
rws_df <- rownames(df_anno_cols_2)

df_anno_cols_2 <- data.frame(df_anno_cols_2[order(as.numeric(as.character(df_anno_cols_2$Clusters)), decreasing = F),])
rownames(df_anno_cols_2) <- rws_df
all(rownames(df_anno_cols_2) == rws_df)
head(df_anno_cols_2)


df_anno_cols_2[1:10,]
Heatmap(heat[1:20,match(rownames(df_anno_cols_2), colnames(heat))],
        show_row_names = F,
        show_column_names = F,
        
        cluster_columns = F,
        col = colorRamp2(breaks = c(-2,0,2),
                         colors = c(scales::alpha(muted("blue"), 0.5), scales::alpha("lightgray",0.3), scales::alpha("tomato", 0.8))),
        top_annotation = HeatmapAnnotation(df = df_anno_cols_2,
                                           col = list(AHH=c("AHH01" = "gray",
                                                            "AHH03" = "green",
                                                            "AHH04" = "blue"),
                                                      Clusters=c("0"= c25[4],
                                                                 "5"= c25[14],
                                                                 "4"= c25[8],
                                                                 "7"= c25[24]))),
        
        use_raster = F,
        name = "Abundance")
table(df_anno_cols_2)

# run DGE for AHH01 ####
Idents(sc_data_no_int) <- sc_data_no_int$PatientID.x
Idents(sc_data_no_int) <- ifelse(sc_data_no_int$PatientID.x %in% c("Pt89", "Pt64"), "Pt89Pt64", "others")
diff_01_pt <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$AHH %in%"AHH01")]), 
                          assay = "SCT",
                          min.cells.group = 3,
                          ident.1 = "Pt89Pt64", ident.2 = "others",
                          logfc.threshold = 0.1, 
                          min.pct = 0.1)
head(diff_01_pt)
diff_01_pt$Gene <- rownames(diff_01_pt)
diff_01_pt$Direction <- ifelse(diff_01_pt$avg_log2FC > 0 & diff_01_pt$p_val_adj < 0.05, "UP",
                               ifelse(diff_01_pt$avg_log2FC < 0 & diff_01_pt$p_val_adj < 0.05, "DOWN", "NS"))

diff_01_pt$Labels <- ifelse(diff_01_pt$Labels == "ACTB" | grepl("^MT-",diff_01_pt$Labels), "", diff_01_pt$Labels)
head(diff_01_pt)
ggplot(diff_01_pt, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_01_pt$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

table(diff_01_pt$Direction)
VlnPlot(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$AHH %in%"AHH01")]),
        group.by = "PatientID.y",
        features = c("CSF2", "CCL4", "GZMB", "CCL3"),
        pt.size = 0.0125,
        ncol=2)

# check Allart's genes in pathways - email 02.04.2022 ####
extra_genes <- readxl::read_xlsx("20220120 List with genes of interest.xlsx",sheet = 1, skip = 1)
sig_all <- c(qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c5.go.bp.v7.2.symbols.gmt"),
             qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.kegg.v7.2.symbols.gmt"),
             qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/h.all.v7.2.symbols.gmt"),
             qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.reactome.v7.2.symbols.gmt"))
sig_all_df <- rbindlist(
  lapply(seq_along(sig_all), function(xx) xx <- data.frame(term = names(sig_all)[[xx]],
                                                           gene = sig_all[[xx]])),
  use.names = T)
head(sig_all_df);tail(sig_all_df)


res_extra_genes_pathways <- clusterProfiler::enricher(gene = extra_genes$FeatureName, 
                                                      TERM2GENE = sig_all_df, 
                                                      qvalueCutoff = 0.05)
dotplot(res_extra_genes_pathways, 
        showCategory =20)
res_extra_genes_pathways <- res_extra_genes_pathways@result[,c(2,7)] %>%
  filter(qvalue < 0.05)

# run GSEA for DGE in previous ppt using pathawys above ####
ranked_genes <- diff_04vs03_cd4_c5_c0$avg_log2FC

ranked_genes_ind <- order(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[ranked_genes_ind]
names(ranked_genes) <- rownames(diff_04vs03_cd4_c5_c0)[ranked_genes_ind]

all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all[which(names(sig_all) %in% res_extra_genes_pathways@result$Description[1:150])],
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      }) %>%
  filter(pval <= 0.1) %>%
  filter(!grepl("NEURO", pathway))
all_fgsea <- all_fgsea[order(all_fgsea$NES, decreasing = T),]
all_fgsea
fgsea::plotEnrichment(sig_all[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                      ranked_genes)
fgsea::plotGseaTable(sig_all[which(names(sig_all) %in% res_extra_genes_pathways@result$Description[which(res_extra_genes_pathways@result$Description %in% all_fgsea$pathway)])],
                     fgseaRes = all_fgsea,gseaParam = 0.5,
                     ranked_genes)

ggplot(all_fgsea %>%
         slice(seq_len(npaths)),
       aes(x = factor(pathway, levels = rev(as.character(pathway))),
           y=NES,
           fill=pval)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  coord_flip() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text_repel(label=all_fgsea$pathway[1:npaths],
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  labs(fill="p.value") +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways") +
  theme(axis.text.x = element_text(size=15))


all_fgsea_04vs03_cd8_c7_c4 <- all_fgsea
all_fgsea_04vs03_cd4_c5_c0 <- all_fgsea

# check cytokine chemokine receptors ####
lig_rec <- fread("~/Google Drive/Protocols_bfx/FANTOM5_ligand_receptor/Reactome_FIsInGene_122220_with_annotations.txt", header = T)
lig_rec <- lig_rec %>%
  filter(Annotation != "predicted")
head(lig_rec)
View(as.data.frame(table(lig_rec$Annotation)))
lig_rec %>%
  filter(Gene2 == "CCR5")
lig_rec <- lig_rec %>%
  filter(Gene1 %in% extra_genes$FeatureName | Gene2 %in% extra_genes$FeatureName)

lig_rec_cd4 <- lig_rec[,c(1,2)] %>%
  filter(Gene1 %in% diff_04vs03_cd4_c5_c0$Gene | Gene2 %in% diff_04vs03_cd4_c5_c0$Gene)
lig_rec_cd8 <- lig_rec[,c(1,2)] %>%
  filter(Gene1 %in% diff_04vs03_cd8_c7_c4$Gene | Gene2 %in% diff_04vs03_cd8_c7_c4$Gene)

lig_rec_cd4_plot <- lig_rec_cd4 %>%
  mutate(Gene1_FC = diff_04vs03_cd4_c5_c0$avg_log2FC[match(Gene1, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(test_Gene1 = diff_04vs03_cd4_c5_c0$Gene[match(Gene1, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(Gene1_sig = ifelse(Gene1 %in% diff_04vs03_cd4_c5_c0$Gene[which(diff_04vs03_cd4_c5_c0$p_val_adj < 0.05)], "YES", "NO")) %>%
  mutate(Gene2_FC = diff_04vs03_cd4_c5_c0$avg_log2FC[match(Gene2, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(test_Gene2 = diff_04vs03_cd4_c5_c0$Gene[match(Gene2, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(Gene2_sig = ifelse(Gene2 %in% diff_04vs03_cd4_c5_c0$Gene[which(diff_04vs03_cd4_c5_c0$p_val_adj < 0.05)], "YES", "NO")) %>%
  mutate(Gene1_Gene2 = paste0(Gene1, "_", Gene2)) %>%
  mutate(Gene1_Gene2_FC = Gene1_FC + Gene2_FC) %>%
  filter(Gene1_sig == "YES" & Gene2_sig == "YES") %>%
  arrange(desc(Gene1_Gene2_FC))


# plot ligand_receptor for CD4
ggplot(lig_rec_cd4_plot,
       aes(x=Gene1_FC,
           y=Gene2_FC,
           col=Gene1_Gene2_FC)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient2(low="#B0E0E6",mid = "lightgray", high="#B22222",midpoint = 0)


ggplot(lig_rec_cd4_plot %>%
         filter(abs(Gene1_Gene2_FC) > 1),
       aes(x = factor(Gene1_Gene2, levels = rev(as.character(Gene1_Gene2))),
           y=Gene1_FC + Gene2_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank())
ggplot(lig_rec_cd4_plot,
       aes(x = factor(paste0(Ligand_symbol, "_", Receptor_symbol), levels = rev(as.character(paste0(Ligand_symbol, "_", Receptor_symbol)))),
           y=Receptor_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  scale_x_discrete(limit = rev(as.character(paste0(lig_rec_cd4_plot$Ligand_symbol, "_", lig_rec_cd4_plot$Receptor_symbol))),
                   labels = rev(as.character(lig_rec_cd4_plot$Receptor_symbol))) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank()) 

ggplot(lig_rec_cd4_plot,
       aes(x = factor(paste0(Ligand_symbol, "_", Receptor_symbol), levels = rev(as.character(paste0(Ligand_symbol, "_", Receptor_symbol)))),
           y=Receptor_FC+Ligand_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           fill="#66CDAA",
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  scale_x_discrete(position="top") +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank()) 

# plot ligand_receptor for CD8
lig_rec_cd8_plot <- lig_rec_cd8 %>%
  mutate(Gene1_FC = diff_04vs03_cd8_c7_c4$avg_log2FC[match(Gene1, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(test_Gene1 = diff_04vs03_cd8_c7_c4$Gene[match(Gene1, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(Gene1_sig = ifelse(Gene1 %in% diff_04vs03_cd8_c7_c4$Gene[which(diff_04vs03_cd8_c7_c4$p_val_adj < 0.05)], "YES", "NO")) %>%
  mutate(Gene2_FC = diff_04vs03_cd8_c7_c4$avg_log2FC[match(Gene2, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(test_Gene2 = diff_04vs03_cd8_c7_c4$Gene[match(Gene2, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(Gene2_sig = ifelse(Gene2 %in% diff_04vs03_cd8_c7_c4$Gene[which(diff_04vs03_cd8_c7_c4$p_val_adj < 0.05)], "YES", "NO")) %>%
  mutate(Gene1_Gene2 = paste0(Gene1, "_", Gene2)) %>%
  mutate(Gene1_Gene2_FC = Gene1_FC + Gene2_FC) %>%
  filter(Gene1_sig == "YES" & Gene2_sig == "YES") %>%
  arrange(desc(Gene1_Gene2_FC))

# plot ligand_receptor for CD4
ggplot(lig_rec_cd8_plot,
       aes(x=Gene1_FC,
           y=Gene2_FC,
           col=Gene1_Gene2_FC)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient2(low="#B0E0E6",mid = "lightgray", high="#B22222",midpoint = 0)


ggplot(lig_rec_cd8_plot %>%
         filter(abs(Gene1_Gene2_FC) > 1),
       aes(x = factor(Gene1_Gene2, levels = rev(as.character(Gene1_Gene2))),
           y=Gene1_FC + Gene2_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank())

ggplot(lig_rec_cd8_plot,
       aes(x = factor(paste0(Ligand_symbol, "_", Receptor_symbol), levels = rev(as.character(paste0(Ligand_symbol, "_", Receptor_symbol)))),
           y=Ligand_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  scale_x_discrete(limit = rev(as.character(paste0(lig_rec_cd8_plot$Ligand_symbol, "_", lig_rec_cd8_plot$Receptor_symbol))),
                   labels = rev(as.character(lig_rec_cd8_plot$Ligand_symbol)),
                   position = "top") +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank())
ggplot(lig_rec_cd8_plot,
       aes(x = factor(paste0(Ligand_symbol, "_", Receptor_symbol), levels = rev(as.character(paste0(Ligand_symbol, "_", Receptor_symbol)))),
           y=Receptor_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  scale_x_discrete(limit = rev(as.character(paste0(lig_rec_cd8_plot$Ligand_symbol, "_", lig_rec_cd8_plot$Receptor_symbol))),
                   labels = rev(as.character(lig_rec_cd8_plot$Receptor_symbol))) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank()) 

ggplot(lig_rec_cd8_plot,
       aes(x = factor(paste0(Ligand_symbol, "_", Receptor_symbol), levels = rev(as.character(paste0(Ligand_symbol, "_", Receptor_symbol)))),
           y=Receptor_FC+Ligand_FC)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           fill="#FFD700",
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  scale_x_discrete(position="top") +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank()) 

# plot whether any ligand or receptor is significantly differentially expressed
lig_rec_cd4_plot <- lig_rec_cd4 %>%
  mutate(Ligand_FC = diff_04vs03_cd4_c5_c0$avg_log2FC[match(Ligand_symbol, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(test_order_ligand = diff_04vs03_cd4_c5_c0$Gene[match(Ligand_symbol, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(Ligand_sig = ifelse(Ligand_symbol %in% diff_04vs03_cd4_c5_c0$Gene[which(diff_04vs03_cd4_c5_c0$p_val_adj < 0.05)], "YES", "NO")) %>%
  mutate(Receptor_FC = diff_04vs03_cd4_c5_c0$avg_log2FC[match(Receptor_symbol, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(test_order_receptor = diff_04vs03_cd4_c5_c0$Gene[match(Receptor_symbol, diff_04vs03_cd4_c5_c0$Gene)]) %>%
  mutate(Receptor_sig = ifelse(Receptor_symbol %in% diff_04vs03_cd4_c5_c0$Gene[which(diff_04vs03_cd4_c5_c0$p_val_adj < 0.05)], "YES", "NO")) %>%
  filter(Ligand_sig == "YES" | Receptor_sig == "YES") %>%
  arrange(desc(Ligand_FC))
lig_rec_cd4_plot <- rbind(lig_rec_cd4_plot[,c(1,5)],
                          lig_rec_cd4_plot[,c(3,8)], use.names=F) %>%
  mutate(Type = ifelse(Ligand_symbol %in% lig_rec_cd4_plot$Ligand_symbol, "Ligand", "Receptor")) %>%
  arrange(desc(Ligand_FC)) %>%
  filter(!is.na(Ligand_FC)) %>%
  filter(!duplicated(Ligand_symbol)) %>%
  rename(GeneName = Ligand_symbol) %>%
  rename(logFC = Ligand_FC)

ggplot(lig_rec_cd4_plot,
       aes(x = factor(GeneName, levels = rev(as.character(GeneName))),
           y=logFC,
           fill=Type)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black", 
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  scale_fill_manual(values=c("Receptor" = "orange", "Ligand" = "gray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  labs(fill="p.value") +
  theme(axis.title.y = element_blank()) 
lig_rec %>%
  filter(Ligand_symbol == "IL7R")

lig_rec_cd8_plot <- lig_rec_cd8 %>%
  mutate(Ligand_FC = diff_04vs03_cd8_c7_c4$avg_log2FC[match(Ligand_symbol, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(test_order_ligand = diff_04vs03_cd8_c7_c4$Gene[match(Ligand_symbol, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(Ligand_sig = ifelse(Ligand_symbol %in% diff_04vs03_cd8_c7_c4$Gene[which(diff_04vs03_cd8_c7_c4$p_val_adj < 0.05)], "YES", "NO")) %>%
  mutate(Receptor_FC = diff_04vs03_cd8_c7_c4$avg_log2FC[match(Receptor_symbol, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(test_order_receptor = diff_04vs03_cd8_c7_c4$Gene[match(Receptor_symbol, diff_04vs03_cd8_c7_c4$Gene)]) %>%
  mutate(Receptor_sig = ifelse(Receptor_symbol %in% diff_04vs03_cd8_c7_c4$Gene[which(diff_04vs03_cd8_c7_c4$p_val_adj < 0.05)], "YES", "NO")) %>%
  filter(Ligand_sig == "YES" & Receptor_sig == "YES") %>%
  arrange(desc(Ligand_FC))


# Allart's email 02.07.2022 ####
# anti-CD3 comparisons. As these are very much determined by the single donor number 4, 
# if we would downsample this donor by 90%, to bring it more in line with the other donors, 
# would that change the analysis?
# generate indices to downsample cell numbers per patient
# analyze pt93 AHH03 633 cells vs.1012 cells
# have to downsample pt4 in ahh03 and ahh04 but also pt93 for ahh04
set.seed(1234)
ind_pt4 <- c(sample(which(sc_data_no_int$PatientID.x == "Pt4" & sc_data_no_int$AHH == "AHH03"), round(3848*0.9, 0), replace = F),
             sample(which(sc_data_no_int$PatientID.x == "Pt4" & sc_data_no_int$AHH == "AHH04"), round(1471*0.9, 0), replace = F))
ind_pt4c0 <- which(sc_data_no_int$AHH == "AHH03" & as.character(sc_data_no_int$seurat_clusters) == "0")
ind_pt4c0 <- ind_pt4c0[-which(ind_pt4c0 %in% ind_pt4)]
sample_pt4c0 <- sample(ind_pt4c0, length(ind_pt4c0)-40, replace = F) 
ind_pt4c4 <- which(sc_data_no_int$AHH == "AHH03" & as.character(sc_data_no_int$seurat_clusters) == "4")
ind_pt4c4 <- ind_pt4c4[-which(ind_pt4c4 %in% ind_pt4)]
sample_pt4c4 <- sample(ind_pt4c4, length(ind_pt4c4)-25, replace = F)
ind_pt4 <- c(ind_pt4, sample_pt4c0, sample_pt4c4)


ind_pt93 <- c(sample(which(sc_data_no_int$PatientID.x == "Pt93" & sc_data_no_int$AHH == "AHH03"), (633*0.75), replace = F),
              sample(which(sc_data_no_int$PatientID.x == "Pt93" & sc_data_no_int$AHH == "AHH04"), (1088*0.75), replace = F))
table(sc_data_no_int$PatientID.x, sc_data_no_int$AHH)
table(sc_data_no_int$PatientID.x[-c(ind_pt4, ind_pt93)], sc_data_no_int$AHH[-c(ind_pt4, ind_pt93)])

as.data.frame.matrix(table(paste0(sc_data_no_int$PatientID.x[-c(ind_pt4, ind_pt93)], "_", sc_data_no_int$AHH[-c(ind_pt4, ind_pt93)]),
                           sc_data_no_int$seurat_clusters[-c(ind_pt4, ind_pt93)])) %>% 
  mutate(PatientID_AHH = rownames(.)) %>%
  melt %>%
  mutate(PtID = sapply(PatientID_AHH, function(x) strsplit(x, "_")[[1]][1])) %>%
  mutate(AHH = sapply(PatientID_AHH, function(x) strsplit(x, "_")[[1]][2])) %>%
  filter(AHH != "AHH01") %>%
  ggplot(., aes(x=PtID,
                y=value,
                fill=AHH)) +
  geom_bar(stat="identity", position = "dodge", width = 0.85, col="black") +
  theme_bw() +
  facet_wrap(~ variable, scales = "free", ncol=5) +
  scale_fill_manual(values = c("AHH03" = muted("yellow"),
                               "AHH04" = muted("tomato")))

# new DGE with downsampled data ####
Idents(sc_data_no_int) <- sc_data_no_int$seurat_clusters
cd4_to_keep <- colnames(sc_data_no_int)[-c(ind_pt4, ind_pt93)]
cd4_to_keep <- cd4_to_keep[which(cd4_to_keep %in% colnames(sc_data_no_int)[which(as.character(sc_data_no_int$seurat_clusters) %in% as.character(c(0,5)))])]
diff_04vs03_cd4_c5_c0 <- FindMarkers(subset(sc_data_no_int, cells=cd4_to_keep), 
                                     assay = "SCT",
                                     # slot = "data",
                                     min.cells.group = 3,
                                     ident.1 = "5", ident.2 = "0",
                                     logfc.threshold = 0.25, 
                                     min.pct = 0.25)
head(diff_04vs03_cd4_c5_c0)
diff_04vs03_cd4_c5_c0$Gene <- rownames(diff_04vs03_cd4_c5_c0)
diff_04vs03_cd4_c5_c0$Direction <- ifelse(diff_04vs03_cd4_c5_c0$avg_log2FC > 0 & diff_04vs03_cd4_c5_c0$p_val_adj < 0.05, "UP",
                                          ifelse(diff_04vs03_cd4_c5_c0$avg_log2FC < 0 & diff_04vs03_cd4_c5_c0$p_val_adj < 0.05, "DOWN", "NS"))

diff_04vs03_cd4_c5_c0$Labels <- ifelse(diff_04vs03_cd4_c5_c0$Labels == "ACTB" | grepl("^MT-",diff_04vs03_cd4_c5_c0$Labels), "", diff_04vs03_cd4_c5_c0$Labels)
head(diff_04vs03_cd4_c5_c0)
ggplot(diff_04vs03_cd4_c5_c0, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_04vs03_cd4_c5_c0$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
diff_04vs03_cd4_c5_c0$InList <- ifelse(diff_04vs03_cd4_c5_c0$Gene %in% extra_genes$FeatureName, "YES","NO")

diff_04vs03_cd4_c5_c0 %>%
  filter(Gene == "IL2RA")

table(diff_04vs03_cd4_c5_c0$InList)
table(diff_04vs03_cd4_c5_c0$InList,
      diff_04vs03_cd4_c5_c0$Direction)

# for CD8: compare cluster 4 (AHH03) vs cluster 7 (AHH04)
cd8_to_keep <- colnames(sc_data_no_int)[-c(ind_pt4, ind_pt93)]
cd8_to_keep <- cd8_to_keep[which(cd8_to_keep %in% colnames(sc_data_no_int)[which(as.character(sc_data_no_int$seurat_clusters) %in% as.character(c(4,7)))])]
table(sc_data_no_int$AHH[cd8_to_keep],
      sc_data_no_int$seurat_clusters[cd8_to_keep])
diff_04vs03_cd8_c7_c4 <- FindMarkers(subset(sc_data_no_int, cells=cd8_to_keep), 
                                     assay = "SCT",
                                     min.cells.group = 2,
                                     ident.1 = "7", ident.2 = "4",
                                     logfc.threshold = 0.25, 
                                     min.pct = 0.25)
head(diff_04vs03_cd8_c7_c4)
diff_04vs03_cd8_c7_c4$Gene <- rownames(diff_04vs03_cd8_c7_c4)
diff_04vs03_cd8_c7_c4$Direction <- ifelse(diff_04vs03_cd8_c7_c4$avg_log2FC > 0 & diff_04vs03_cd8_c7_c4$p_val_adj < 0.05, "UP",
                                          ifelse(diff_04vs03_cd8_c7_c4$avg_log2FC < 0 & diff_04vs03_cd8_c7_c4$p_val_adj < 0.05, "DOWN", "NS"))

diff_04vs03_cd8_c7_c4$Labels <- ifelse(diff_04vs03_cd8_c7_c4$Labels == "ACTB" | grepl("^MT-",diff_04vs03_cd8_c7_c4$Labels), "", diff_04vs03_cd8_c7_c4$Labels)
head(diff_04vs03_cd8_c7_c4)

diff_04vs03_cd8_c7_c4 %>%
  filter(Gene == "GZMB")

ggplot(diff_04vs03_cd8_c7_c4, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=diff_04vs03_cd8_c7_c4$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))
diff_04vs03_cd8_c7_c4$InList <- ifelse(diff_04vs03_cd8_c7_c4$Gene %in% extra_genes$FeatureName, "YES","NO")
table(diff_04vs03_cd8_c7_c4$InList)
table(diff_04vs03_cd8_c7_c4$InList,
      diff_04vs03_cd8_c7_c4$Direction)


# run GSEA for DGE above ####
ranked_genes <- diff_04vs03_cd4_c5_c0$avg_log2FC
ranked_genes_ind <- order(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[ranked_genes_ind]
names(ranked_genes) <- rownames(diff_04vs03_cd4_c5_c0)[ranked_genes_ind]

sig_all <- c(qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c5.go.bp.v7.2.symbols.gmt"),
             qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.kegg.v7.2.symbols.gmt"),
             qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/h.all.v7.2.symbols.gmt"),
             qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.reactome.v7.2.symbols.gmt"))

all_fgsea <- tryCatch(fgsea::fgsea(pathways = sig_all[which(names(sig_all) %in% res_extra_genes_pathways@result$Description[1:150])],
                                   stats = ranked_genes),
                      error = function(xx){
                        message(xx)
                        message("\nadding NAs\n")
                        dummy_list <- NA
                        return(dummy_list)
                      }) %>%
  filter(pval <= 0.1)
all_fgsea %>%
  View

all_fgsea <- all_fgsea[order(all_fgsea$NES, decreasing = T),]

npaths <- nrow(all_fgsea)
ggplot(all_fgsea %>%
         slice(seq_len(npaths)),
       aes(x = factor(pathway, levels = rev(as.character(pathway))),
           y=NES,
           fill = pval)) +
  geom_bar(stat="identity", 
           position = "dodge", 
           col="black",
           alpha=0.5,
           width = 0.8) +
  theme_classic() +
  coord_flip() +
  scale_fill_gradient(low = "lightgray", high = "tomato") +
  geom_text_repel(label=all_fgsea$pathway[1:npaths],
                  box.padding = 0.1, 
                  force = 0,
                  min.segment.length = 0,
                  max.overlaps = 30,
                  fontface = 'bold', 
                  color = 'black') +
  theme(axis.text.y = element_blank()) +
  xlab("Pathways")


# change TRBV grouping ####
# DEG tables for the analyses that were shared in slides from Seb dated 20th Jan, 
# comparing specific TRBV subsets (6-1, 6-5, 20-1); 
# this will be helpful to cross-compare, for the validations.

table(sc_data_no_int@meta.data$TRBV_used_grouped[which(sc_data_no_int$CD8_CD4 == "CD4"  & sc_data_no_int$AHH == "AHH04")])
table(sc_data_no_int@meta.data$TRBV_used_grouped[which(sc_data_no_int$CD8_CD4 == "CD8"  & sc_data_no_int$AHH == "AHH04")])

table(sc_data_no_int$TRBV_used)
sc_data_no_int$TRBV_used_grouped <- ifelse(sc_data_no_int$TRBV_used.v2 == "TRBV20", "TRBV20",
                                           ifelse(sc_data_no_int$TRBV_used.v2 %in% c("TRBV6.1", "TRBV6.2", "TRBV6.5"), "TRBV6",
                                                  sc_data_no_int$TRBV_used.v2))
table(sc_data_no_int$TRBV_used_grouped)
Idents(sc_data_no_int) <- sc_data_no_int$TRBV_used_grouped
res_dge_cd8_6.1_20_04_grouped <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 == "CD8"  & sc_data_no_int$AHH == "AHH04")]), 
                                             min.cells.group = 3,
                                             ident.1 = "TRBV6", ident.2 = "TRBV20",
                                             logfc.threshold = 0.25, 
                                             min.pct = 0)
head(res_dge_cd8_6.1_20_04_grouped)
ggplot(res_dge_cd8_6.1_20_04_grouped %>%
         mutate(Direction = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                                   ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DOWN", "NS"))),
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
  geom_point() + 
  theme_pubclean() +
  scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                "DOWN" = alpha("lightblue", 0.5),
                                "NS" = "lightgray")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(label=ifelse(res_dge_cd8_6.1_20_04_grouped$p_val_adj < 0.05, rownames(res_dge_cd8_6.1_20_04_grouped), ""),
                  box.padding = 0.3, 
                  force = 0.15,
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3))) +
  ggtitle("")




res_dge_cd4_6.1_20_04_grouped <- FindMarkers(subset(sc_data_no_int, cells=colnames(sc_data_no_int)[which(sc_data_no_int$CD8_CD4 == "CD4"  & sc_data_no_int$AHH == "AHH04")]), 
                                             min.cells.group = 3,
                                             ident.1 = "TRBV6", ident.2 = "TRBV20",
                                             logfc.threshold = 0.25, 
                                             min.pct = 0)
head(res_dge_cd4_6.1_20_04_grouped)

ggplot(res_dge_cd4_6.1_20_04_grouped %>%
         mutate(Direction = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                                   ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DOWN", "NS"))),
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) + 
  geom_point() + 
  theme_pubclean() +
  scale_color_manual(values = c("UP" = alpha("red", 0.5), 
                                "DOWN" = alpha("lightblue", 0.5),
                                "NS" = "lightgray")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(label=ifelse(res_dge_cd4_6.1_20_04_grouped$p_val_adj < 0.05, rownames(res_dge_cd4_6.1_20_04_grouped), ""),
                  box.padding = 0.3, 
                  force = 0.15,
                  min.segment.length = 0,
                  max.overlaps = 20,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3))) +
  ylim(0,30) +
  ggtitle("") 


# create cloupe files for all and AHH04 seaprately ####
library(funcutils)
DimPlot(sc_04, reduction = "UMAP_prot_o4")
seurat2cloupe(sc_04,
              dims = c(1,2),
              reduction = "UMAP_prot_o4",
              keyword = "sc_AHH04",
              metadata = colnames(sc_04@meta.data),
              opdir = "."
)
DimPlot(sc_data_no_int, reduction = "prot_UMAP")
seurat2cloupe(sc_data_no_int,
              dims = c(1,2),
              reduction = "prot_UMAP",
              keyword = "sc_data_all_AHH",
              metadata = colnames(sc_data_no_int@meta.data),
              opdir = "."
)


dd <- read.table("res_dge_cd8_6.5_04vs01.txt", header = T, sep = "\t")
head(dd)
ggplot(dd, aes(x = pct.1 - pct.2,
               y = avg_log2FC,
               col=p_val_adj)) +
  geom_point() +
  theme_bw() +
  geom_text_repel(label=dd$X,
                  box.padding = 0.3, 
                  min.segment.length = 0,
                  max.overlaps = 10,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white")


# check KLRG1 as in email 03302022 ####
head(sc_data_no_int)
Idents(sc_data_no_int) <- sc_data_no_int$Protein_snn_res.0.25
table(sc_data_no_int$CD8_CD4)
VlnPlot(sc_data_no_int %>% subset(cells = which(sc_data_no_int$CD8_CD4 == "CD8")),
        idents = c(7,4),
        features = "KLRG1", 
        log = T,
        pt.size = 1) +
  stat_compare_means(method="t.test")

mean(sc_data_no_int@assays$SCT@data["KLRG1", which(sc_data_no_int$CD8_CD4 == "CD8" & sc_data_no_int$Protein_snn_res.0.25 == "7")]) - mean(sc_data_no_int@assays$SCT@data["KLRG1", which(sc_data_no_int$CD8_CD4 == "CD8" & sc_data_no_int$Protein_snn_res.0.25 == "4")])
mean(sc_data_no_int@assays$SCT@data["KLRG1", which(sc_data_no_int$CD8_CD4 == "CD4" & sc_data_no_int$Protein_snn_res.0.25 == "5")]) - mean(sc_data_no_int@assays$SCT@data["KLRG1", which(sc_data_no_int$CD8_CD4 == "CD4" & sc_data_no_int$Protein_snn_res.0.25 == "0")])


# DEG within AHH04 ####
# Using the AHH04 treatment group, 
# compare on the one hand TRBV6+ (including 6-1, 6-2, 6-5 and 10-3) 
# versus TRBV6- (any cell expressing a TRBV but not those in the TRBV6+ group), 
# separately for CD4 and CD8 populations.
table(sc_data_no_int$TRBV_used_grouped)
sc_data_no_int$TRBV_for_dge <- ifelse(sc_data_no_int$TRBV_used.v3 %in% c("TRBV10.3", "TRBV6.1", "TRBV6.2", "TRBV6.5"), "test", "control")
table(sc_data_no_int$TRBV_for_dge)
Idents(sc_data_no_int) <- sc_data_no_int$TRBV_for_dge
dge_trbv_6s10_vs_other <- FindMarkers(sc_data_no_int %>% subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04")]),
                                      # slot = "data",
                                      min.cells.group = 5,
                                      ident.1 = "test", ident.2 = "control",
                                      logfc.threshold = 0.2, 
                                      min.pct = 0.3)
dge_trbv_6s10_vs_other$Gene <- rownames(dge_trbv_6s10_vs_other)
write.table(dge_trbv_6s10_vs_other, "dge_trbv_6s10_vs_other_03312022.txt", sep="\t", quote = F, row.names = F)

dge_trbv_6s10_vs_other_cd8 <- FindMarkers(sc_data_no_int %>% subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD8")]),
                                          # slot = "data",
                                          min.cells.group = 5,
                                          ident.1 = "test", ident.2 = "control",
                                          logfc.threshold = 0.2, 
                                          min.pct = 0.3)
dge_trbv_6s10_vs_other_cd8$Gene <- rownames(dge_trbv_6s10_vs_other_cd8)
dge_trbv_6s10_vs_other_cd8$Direction <- ifelse(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05 & dge_trbv_6s10_vs_other_cd8$avg_log2FC > 0, "UP",
                                               ifelse(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05 & dge_trbv_6s10_vs_other_cd8$avg_log2FC < 0, "DOWN", "NS"))
write.table(dge_trbv_6s10_vs_other_cd8, "dge_trbv_6s10_vs_other_cd8_03312022.txt", sep="\t", quote = F, row.names = F)

dge_trbv_6s10_vs_other_cd4 <- FindMarkers(sc_data_no_int %>% subset(cells = colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD4")]),
                                          # slot = "data",
                                          min.cells.group = 5,
                                          ident.1 = "test", ident.2 = "control",
                                          logfc.threshold = 0.2, 
                                          min.pct = 0.3)
dge_trbv_6s10_vs_other_cd4$Gene <- rownames(dge_trbv_6s10_vs_other_cd4)
dge_trbv_6s10_vs_other_cd4$Direction <- ifelse(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05 & dge_trbv_6s10_vs_other_cd4$avg_log2FC > 0, "UP",
                                               ifelse(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05 & dge_trbv_6s10_vs_other_cd4$avg_log2FC < 0, "DOWN", "NS"))
write.table(dge_trbv_6s10_vs_other_cd4, "dge_trbv_6s10_vs_other_cd4_03312022.txt", sep="\t", quote = F, row.names = F)

all_TF <- read.table("~/Google Drive/Protocols_bfx/GSEA_signatures/all_TF.txt", sep="\t", header = F)[,1]

dge_trbv_6s10_vs_other_cd8 %>%
  filter(Gene %in% all_TF) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")

dge_trbv_6s10_vs_other_cd4 %>%
  filter(Gene %in% all_TF) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")

table(sc_data_no_int$TRBV_used.v3[which(sc_data_no_int$CD8_CD4 == "CD4")],
      sc_data_no_int$AHH[which(sc_data_no_int$CD8_CD4 == "CD4")])

# do heatmap as we discussed in call 03312022 ####
# take AHH04 and annotate TRBV+ and TRBV-
table(sc_data_no_int$TRBV_for_dge)
df_anno <- data.frame(TRBV = ifelse(sc_data_no_int$TRBV_for_dge == "test", "TRBVpos", "TRBVneg"),
                      AHH = sc_data_no_int$AHH)
df_anno <- df_anno %>% filter(AHH == "AHH04")
head(df_anno)

# heat 1
# 1.	whole list in xls file
table(df_anno$TRBV)
genes4heat <- readxl::read_xlsx("~/Google Drive/BridgeBfx_SB/MarengoTx/Genes of interest_03312022.xlsx", skip = 1, col_names = T)
genes4heat$FeatureName[-which(c(genes4heat$FeatureName, "IL8RA") %in% rownames(sc_data_no_int@assays$SCT@data))]
heat_1 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% c(genes4heat$FeatureName, "IL8RA")),
                                         which(sc_data_no_int$AHH == "AHH04")]
cls <- colnames(heat_1)
heat_1 <- data.matrix(t(apply(heat_1, 1, scale)))
colnames(heat_1) <- cls
heat_1 <- heat_1[-which(apply(heat_1, 1, function(xx) sum(is.na(xx))) == ncol(heat_1)),]
heat_1[1:5,1:5]
dim(heat_1)
Heatmap(heat_1, 
        show_column_names = F,
        show_row_names = T, 
        top_annotation = HeatmapAnnotation(df = df_anno, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "spearman",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")

# heat 2
# 2.	CD4 with DGE
dge_trbv_6s10_vs_other_cd4$inXlsFile <- ifelse(dge_trbv_6s10_vs_other_cd4$Gene %in% genes4heat$FeatureName, "YES", "NO")
write.table(dge_trbv_6s10_vs_other_cd4, "dge_trbv_6s10_vs_other_cd4_w_xls.txt", sep = "\t", quote=F, row.names = F)

genes4heat2 <- dge_trbv_6s10_vs_other_cd4$Gene[which(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05)]
heat_2 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat2),
                                         which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD4")]
cls <- colnames(heat_2)
heat_2 <- data.matrix(t(apply(heat_2, 1, scale)))
colnames(heat_2) <- cls
heat_2[1:5,1:5]
dim(heat_2)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD4")]),]
Heatmap(heat_2, 
        show_column_names = F,
        show_row_names = F, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "pearson",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")

# heat 3
# 3.	CD8 with DGE

dge_trbv_6s10_vs_other_cd8$inXlsFile <- ifelse(dge_trbv_6s10_vs_other_cd8$Gene %in% genes4heat$FeatureName, "YES", "NO")
write.table(dge_trbv_6s10_vs_other_cd8, "dge_trbv_6s10_vs_other_cd8_w_xls.txt", sep = "\t", quote=F, row.names = F)
genes4heat3 <- dge_trbv_6s10_vs_other_cd8$Gene[which(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05)]
heat_3 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat3),
                                         which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD8")]
cls <- colnames(heat_3)
heat_3 <- data.matrix(t(apply(heat_3, 1, scale)))
colnames(heat_3) <- cls
heat_3[1:5,1:5]
dim(heat_3)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD8")]),]
Heatmap(heat_3, 
        show_column_names = F,
        show_row_names = F, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "pearson",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")

# heat 4
# 4.	CD4 with DGE overlapping xls file
dge_trbv_6s10_vs_other_cd4$inXlsFile <- ifelse(dge_trbv_6s10_vs_other_cd4$Gene %in% genes4heat$FeatureName, "YES", "NO")
write.table(dge_trbv_6s10_vs_other_cd4, "dge_trbv_6s10_vs_other_cd4_w_xls.txt", sep = "\t", quote=F, row.names = F)
genes4heat4 <- intersect(dge_trbv_6s10_vs_other_cd4$Gene[which(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05)],
                         c(genes4heat$FeatureName, "IL8RA"))
heat_4 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat4),
                                         which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD4")]
cls <- colnames(heat_4)
heat_4 <- data.matrix(t(apply(heat_4, 1, scale)))
colnames(heat_4) <- cls
heat_4[1:5,1:5]
dim(heat_4)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD4")]),]
Heatmap(heat_4, 
        show_column_names = F,
        show_row_names = T, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "pearson",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")

# heat 5
# 5.	CD8 with DGE overlapping xls file
dge_trbv_6s10_vs_other_cd8$inXlsFile <- ifelse(dge_trbv_6s10_vs_other_cd8$Gene %in% genes4heat$FeatureName, "YES", "NO")
write.table(dge_trbv_6s10_vs_other_cd8, "dge_trbv_6s10_vs_other_cd8_w_xls.txt", sep = "\t", quote=F, row.names = F)
genes4heat5 <- intersect(c(genes4heat$FeatureName, "IL8RA"), 
                         dge_trbv_6s10_vs_other_cd8$Gene[which(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05)])
heat_5 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat5),
                                         which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD8")]
cls <- colnames(heat_5)
heat_5 <- data.matrix(t(apply(heat_5, 1, scale)))
colnames(heat_5) <- cls
heat_5[1:5,1:5]
dim(heat_5)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 == "CD8")]),]
Heatmap(heat_5, 
        show_column_names = F,
        show_row_names = T, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "pearson",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")

# heat 6
# 6.	CD4 DGE that overlap with CD8 DGE and then overlap with xls file
genes4heat6 <- intersect(c(genes4heat$FeatureName, "IL8RA"), 
                         intersect(dge_trbv_6s10_vs_other_cd4$Gene[which(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05)],
                                   dge_trbv_6s10_vs_other_cd8$Gene[which(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05)]))
heat_6 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat6),
                                         which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 %in% c("CD8", "CD4"))]
cls <- colnames(heat_6)
heat_6 <- data.matrix(t(apply(heat_6, 1, scale)))
colnames(heat_6) <- cls
heat_6[1:5,1:5]
dim(heat_6)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 %in% c("CD8", "CD4"))]),]
Heatmap(heat_6, 
        show_column_names = F,
        show_row_names = T, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "spearman",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")


# heat 7
# 7.	CD4 DGE overlap with CD8 DGE
genes4heat7 <- intersect(dge_trbv_6s10_vs_other_cd4$Gene[which(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05)],
                         dge_trbv_6s10_vs_other_cd8$Gene[which(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05)])
heat_7 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat7),
                                         which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4 %in% c("CD8", "CD4"))]
cls <- colnames(heat_7)
heat_7 <- data.matrix(t(apply(heat_7, 1, scale)))
colnames(heat_7) <- cls
heat_7[1:5,1:5]
dim(heat_7)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$CD8_CD4  %in% c("CD8", "CD4"))]),]
Heatmap(heat_7, 
        show_column_names = F,
        show_row_names = F, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen"))),
        raster_quality = 2, 
        clustering_distance_columns = "spearman",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")


# check the overall overlap between CD4 and CD8 DEGs 
par(mar=c(0,0,0,0))
gplots::venn(list(CD4 = dge_trbv_6s10_vs_other_cd4$Gene[which(dge_trbv_6s10_vs_other_cd4$p_val_adj < 0.05)],
                  CD8 = dge_trbv_6s10_vs_other_cd8$Gene[which(dge_trbv_6s10_vs_other_cd8$p_val_adj < 0.05)],
                  XlsList = genes4heat$FeatureName))

# make volcano
head(dge_trbv_6s10_vs_other_cd4)
ggplot(dge_trbv_6s10_vs_other_cd4, 
       aes(x=avg_log2FC, 
           y=-log10(p_val_adj), 
           col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=dge_trbv_6s10_vs_other_cd4$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  ggtitle("dge_trbv_6s10_vs_other_cd4") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

ggplot(dge_trbv_6s10_vs_other_cd4, 
       aes(x=avg_log2FC, 
           y=-log10(p_val_adj), 
           col=inXlsFile)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("YES" = "green", "NO"=muted("blue"))) +
  geom_text_repel(label=dge_trbv_6s10_vs_other_cd4$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  ggtitle("dge_trbv_6s10_vs_other_cd4") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))


ggplot(dge_trbv_6s10_vs_other_cd8, 
       aes(x=avg_log2FC, 
           y=-log10(p_val_adj), 
           col=Direction)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  geom_text_repel(label=dge_trbv_6s10_vs_other_cd8$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  ggtitle("dge_trbv_6s10_vs_other_cd8") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))

ggplot(dge_trbv_6s10_vs_other_cd8, 
       aes(x=avg_log2FC, 
           y=-log10(p_val_adj), 
           col=inXlsFile)) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values = c("YES" = "green", "NO"=muted("blue"))) +
  geom_text_repel(label=dge_trbv_6s10_vs_other_cd8$Gene,
                  box.padding = 0.3, 
                  force = 0.1,
                  min.segment.length = 0,
                  max.overlaps = 15,
                  fontface = 'bold', 
                  color = 'black',
                  fill = "white") +
  ggtitle("dge_trbv_6s10_vs_other_cd8") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=3)))


# email from Allart 04/11/2022 ####
# Some of these analyses have already been run, e.g. DGE_diff_04vs03_cd4_c5_c0_02032022 and associated cd8 file comparing c7 vs c4. From this file, 
# could you please plot the TFs, just as you did for the TRBV6s, and heatmaps of all significant genes (CD4, CD8 and overlapping genes

# TF plots
diff_04vs03_cd4_c5_c0[which(diff_04vs03_cd4_c5_c0$Gene %in% c("COPS2", "ENO1", "TCF7")),]
diff_04vs03_cd4_c5_c0 %>%
  filter(Gene %in% all_TF) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")

diff_04vs03_cd8_c7_c4 %>%
  filter(Gene %in% all_TF) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")


# make heatmaps
# take AHH04 and annotate TRBV+ and TRBV-
table(sc_data_no_int$TRBV_for_dge)
df_anno <- data.frame(TRBV = ifelse(sc_data_no_int$TRBV_for_dge == "test", "TRBVpos", "TRBVneg"),
                      AHH = sc_data_no_int$AHH)
df_anno <- df_anno %>% filter(AHH %in% c("AHH03","AHH04"))
head(df_anno)

# heat 1
# 1.	CD4 with DGE
genes4heat_cd4 <- diff_04vs03_cd4_c5_c0$Gene[which(diff_04vs03_cd4_c5_c0$p_val_adj < 0.05)]
heat_1 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat_cd4),
                                         which(sc_data_no_int$AHH %in% c("AHH03","AHH04") & sc_data_no_int$CD8_CD4 == "CD4")]
cls <- colnames(heat_1)
heat_1 <- data.matrix(t(apply(heat_1, 1, scale)))
colnames(heat_1) <- cls
heat_1[1:5,1:5]
dim(heat_1)
df_anno1 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH %in% c("AHH03","AHH04") & sc_data_no_int$CD8_CD4 == "CD4")]),]
Heatmap(heat_1, 
        show_column_names = F,
        show_row_names = F, 
        top_annotation = HeatmapAnnotation(df = df_anno1, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen",
                                                              "AHH03" = "lightblue"))),
        raster_quality = 2, 
        clustering_distance_columns = "pearson",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")


# add TF (per email 05.16.2022)
# check TF ####
diff_04vs03_cd4_c5_c0 %>%
  filter(Gene %in% all_TF) %>%
  filter(Gene %in% rownames(heat_1)) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")


# heat 2
# 2.	CD8 with DGE
dim(diff_04vs03_cd8_c7_c4)
diff_04vs03_cd8_c7_c4 <- read.table("DGE_diff_04vs03_cd8_c7_c4_02032022.txt", header = T, sep="\t")
genes4heat_cd8 <- diff_04vs03_cd8_c7_c4$Gene[which(diff_04vs03_cd8_c7_c4$p_val_adj < 0.05)]
heat_2 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat_cd8),
                                         which(sc_data_no_int$AHH %in% c("AHH03","AHH04") & sc_data_no_int$CD8_CD4 == "CD8")]
cls <- colnames(heat_2)
heat_2 <- data.matrix(t(apply(heat_2, 1, scale)))
colnames(heat_2) <- cls
heat_1[1:5,1:5]
dim(heat_2)
df_anno2 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH %in% c("AHH03","AHH04") & sc_data_no_int$CD8_CD4 == "CD8")]),]
Heatmap(heat_2, 
        show_column_names = F,
        show_row_names = F, 
        top_annotation = HeatmapAnnotation(df = df_anno2, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen",
                                                              "AHH03" = "lightblue"))),
        row_km = 4,
        raster_quality = 2, 
        clustering_distance_columns = "euclidean",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")


# cehck TF ####
diff_04vs03_cd8_c7_c4 %>%
  filter(Gene %in% all_TF) %>%
  filter(Gene %in% rownames(heat_2)) %>%
  arrange(desc(avg_log2FC)) %>%
  ggplot(., aes(x = factor(Gene, levels = rev(as.character(Gene))),
                y = avg_log2FC,
                fill=Direction)) +
  geom_bar(stat="identity", position = "dodge", width = 0.75, col="black") +
  theme_bw() +
  scale_fill_manual(values = c("UP" = muted("red"), "DOWN"=muted("cyan"), "NS"="lightgray")) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip() +
  xlab("")


# heat 3
# 3.	overlap
genes4heat_cd8_cd4 <- intersect(genes4heat_cd8, genes4heat_cd4)
heat_3 <- sc_data_no_int@assays$SCT@data[which(rownames(sc_data_no_int@assays$SCT@data) %in% genes4heat_cd8_cd4),
                                         which(sc_data_no_int$AHH %in% c("AHH03","AHH04") & sc_data_no_int$CD8_CD4 %in% c("CD8", "CD4"))]
cls <- colnames(heat_3)
heat_3 <- data.matrix(t(apply(heat_3, 1, scale)))
colnames(heat_3) <- cls
dim(heat_3)
df_anno3 <- df_anno[which(rownames(df_anno) %in% colnames(sc_data_no_int)[which(sc_data_no_int$AHH %in% c("AHH03","AHH04") & sc_data_no_int$CD8_CD4 %in% c("CD8", "CD4"))]),]
df_anno3$CD4_CD8 <- ifelse(rownames(df_anno3) %in% rownames(df_anno1), "CD4", "CD8")
head(df_anno3)
Heatmap(heat_3, 
        show_column_names = F,
        show_row_names = F, 
        top_annotation = HeatmapAnnotation(df = df_anno3, 
                                           col = list(TRBV = c("TRBVpos" = "orange",
                                                               "TRBVneg" = "lightgray"),
                                                      AHH = c("AHH04" = "lightgreen",
                                                              "AHH03" = "lightblue"),
                                                      CD4_CD8 = c("CD4" = "purple", 
                                                                  "CD8" = "#DAA520"))),
        raster_quality = 2, 
        clustering_distance_columns = "pearson",
        col= colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("#B0E0E6", "white", "#DC143C")),
        border = T,
        name = "Expression")

# Prep intersection between cd4 and cd8 with stats for each.
head(diff_04vs03_cd8_c7_c4)
head(diff_04vs03_cd4_c5_c0)
df_both_stat <- data.frame(cbind(diff_04vs03_cd8_c7_c4[which(diff_04vs03_cd8_c7_c4$Gene %in% rownames(heat_3)),],
                                 diff_04vs03_cd4_c5_c0[which(diff_04vs03_cd4_c5_c0$Gene %in% rownames(heat_3)),]))
df_both_stat[,10:ncol(df_both_stat)] <- df_both_stat[match(df_both_stat$Gene, df_both_stat$Gene.1),10:ncol(df_both_stat)]
df_both_stat$Gene <- NULL
df_both_stat$Gene.1 <- NULL
df_both_stat$InList <- NULL
colnames(df_both_stat)[1] <- "GeneName"
colnames(df_both_stat) <- c(colnames(df_both_stat)[1],
                            paste0(colnames(df_both_stat)[2:7], "_CD8"),
                            paste0(colnames(df_both_stat)[8:13], "_CD4"))
head(df_both_stat)
write.table(df_both_stat, "results_intersection_cd8_cd4_05192022.txt", row.names = F, sep = "\t", quote = F)

# email allart 05.17.2022 ####
# extract cluster 
ha_cluster <- Heatmap(heat_2, 
                      show_column_names = F,
                      show_row_names = F, 
                      top_annotation = HeatmapAnnotation(df = df_anno2, 
                                                         col = list(TRBV = c("TRBVpos" = "orange",
                                                                             "TRBVneg" = "lightgray"),
                                                                    AHH = c("AHH04" = "lightgreen",
                                                                            "AHH03" = "lightblue"))),
                      row_km = 4,
                      raster_quality = 2, 
                      clustering_distance_columns = "euclidean",
                      col= colorRamp2(breaks = c(-1.5,0,1.5),
                                      colors = c("#B0E0E6", "white", "#DC143C")),
                      border = T,
                      name = "Expression")

r.dend <- row_dend(ha_cluster)  #If needed, extract row dendrogram
rcl.list <- row_order(ha_cluster)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%>%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(heat_2[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

#check
clu_df
table(clu_df$Cluster)
genes2keep <- clu_df$GeneID[which(clu_df$Cluster == "cluster1")]

# other question:
# plot in something resembling example below, but where we have eight columns: 
# 1 -- TRBV6+ CD4; 
# 2 -- TRBV6+ CD8; 
# 3 -- TRBV6- CD4; 
# 4 -- TRBV6- CD8; 
# 5 -- AHH04 CD4 (cluster 5); 
# 6 -- AHH04 CD8 (cluster 7); 
# 7 -- AHH03 CD4 (cluster 0) and 
# 8 -- AHH03 CD8 (cluster 4)
head(sc_data_no_int)
table(sc_data_no_int$TRBV_used.v3)
table(sc_data_no_int$TRBV_used.v3,
      sc_data_no_int$TRBV_for_dge)
table(sc_data_no_int$TRBV_used.v3,
      sc_data_no_int$AHH)


groups_4_heat <- list(AntiVB_CD4 = which(sc_data_no_int$TRBV_for_dge == "test" & sc_data_no_int$CD8_CD4 == "CD4" & sc_data_no_int$AHH == "AHH04"),
                      AntiVB_CD8 = which(sc_data_no_int$TRBV_for_dge == "test" & sc_data_no_int$CD8_CD4 == "CD8" & sc_data_no_int$AHH == "AHH04"),
                      AntiVB_neg_CD4 = which(sc_data_no_int$TRBV_for_dge == "control" & sc_data_no_int$CD8_CD4 == "CD4" & sc_data_no_int$AHH == "AHH04"),
                      AntiVB_neg_CD8 = which(sc_data_no_int$TRBV_for_dge == "control" & sc_data_no_int$CD8_CD4 == "CD8" & sc_data_no_int$AHH == "AHH04"),
                      AntiVB_c5_CD4 = which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$seurat_clusters == "5" & sc_data_no_int$CD8_CD4 == "CD4"),
                      AntiVB_c7_CD8 = which(sc_data_no_int$AHH == "AHH04" & sc_data_no_int$seurat_clusters == "7" & sc_data_no_int$CD8_CD4 == "CD8"),
                      AntiCD3_c0_CD4 = which(sc_data_no_int$AHH == "AHH03" & sc_data_no_int$seurat_clusters == "0" & sc_data_no_int$CD8_CD4 == "CD4"),
                      AntiCD3_c4_CD8 = which(sc_data_no_int$AHH == "AHH03" & sc_data_no_int$seurat_clusters == "4" & sc_data_no_int$CD8_CD4 == "CD8")
)


lapply(groups_4_heat, length)


# extra genes 05.20.2022 ####
genes2use <- readxl::read_xlsx("20220520 Additional genes for plotting.xlsx", 
                               sheet = 1)$Gene
genes2usetest <- c(genes2use, "LEF1", "MKI67", "PRF1", "SCL6A3")
genes2usetest[-which(genes2usetest %in% rownames(sc_data_no_int))]

ind_genes <- which(rownames(sc_data_no_int) %in% genes2use)
for(ii in 1:length(groups_4_heat)){
  temp_heat <- data.matrix(sc_data_no_int@assays$SCT@counts[ind_genes,groups_4_heat[[ii]]])
  temp_val <- apply(temp_heat, 1, mean)
  if(ii == 1){
    df_heat <- matrix(nrow = length(genes2use),
                      ncol=length(groups_4_heat),dimnames = list(genes2use,
                                                                 names(groups_4_heat)))
  }
  df_heat[,ii] <- temp_val[match(rownames(df_heat), names(temp_val))]
}

head(df_heat)
par(mar=c(8,4,2,2))
barplot(df_heat["IL2RB",], names=colnames(df_heat), las=2)
heat <- t(apply(df_heat, 1, scale))
colnames(heat) <- names(groups_4_heat)
Heatmap(heat,
        cluster_rows = T,
        cluster_columns = T, 
        border = T, 
        column_km = 2, 
        row_km = 3,
        rect_gp = gpar(col = "darkgray", lwd = 1),
        name = "scaled_expression", 
        show_row_dend = F, 
        show_column_dend = F)


# do venn diagrams 
# anti-VB against anti-CD3 (clusters 0 vs 5 and 4 vs 7)
par(mar=c(1,1,1,1))
venn(list(CD4=diff_04vs03_cd4_c5_c0$Gene,
          CD8=diff_04vs03_cd8_c7_c4$Gene),
     main="04 vs. 03 - CD4 vs. CD8")


# address email 05.23.22 ####
# splitting the graphs up into figures only showing a comparison between three groups being either CD4 or CD8:
# Comparing: antiVB_CD4; antiCD3_c0_CD4 and antiVBneg_CD4
# Comparing: antiVB_CD8, antiCD3_c4_CD8 and antiVBneg_CD8
# Could you please do this for the genes in both gene lists

# extra genes 05.20.2022 ####
# genes2use <- readxl::read_xlsx("20220517 Selection cytokine chemokines IR DGEs.xlsx",
#                                sheet = 5)$Gene
genes2use <- readxl::read_xlsx("20220520 Additional genes for plotting.xlsx",
                               sheet = 1)$Gene
genes2usetest <- c(genes2use, "LEF1", "MKI67", "PRF1", "SLC2A3")
genes2usetest[-which(genes2usetest %in% rownames(sc_data_no_int))]

ind_genes <- which(rownames(sc_data_no_int) %in% genes2use)
for(ii in 1:length(groups_4_heat)){
  temp_heat <- data.matrix(sc_data_no_int@assays$SCT@counts[ind_genes,groups_4_heat[[ii]]])
  temp_val <- apply(temp_heat, 1, mean)
  if(ii == 1){
    df_heat <- matrix(nrow = length(genes2use),
                      ncol=length(groups_4_heat),dimnames = list(genes2use,
                                                                 names(groups_4_heat)))
  }
  df_heat[,ii] <- temp_val[match(rownames(df_heat), names(temp_val))]
}

head(df_heat)
par(mar=c(8,4,2,2))
barplot(df_heat["IL2RB",], names=colnames(df_heat), las=2)


# h1 = antiVB_CD4; antiCD3_c0_CD4 and antiVBneg_CD4
ind_cols_heat <- which(colnames(df_heat) %in% c("AntiVB_CD4", "AntiCD3_c0_CD4", "AntiVB_neg_CD4"))
heat <- t(apply(df_heat[,ind_cols_heat], 1, scale))

colnames(heat) <- names(groups_4_heat)[ind_cols_heat]
head(heat)
Heatmap(heat,
        cluster_rows = T,
        cluster_columns = T, 
        border = T, 
        column_km = 2, 
        row_km = 4,
        rect_gp = gpar(col = "darkgray", lwd = 1),
        name = "scaled_expression", 
        show_row_dend = F, 
        show_column_dend = F)


ind_cols_heat <- which(colnames(df_heat) %in% c("AntiVB_CD8", "AntiCD3_c4_CD8", "AntiVB_neg_CD8"))
heat <- t(apply(df_heat[,ind_cols_heat], 1, scale))

colnames(heat) <- names(groups_4_heat)[ind_cols_heat]
head(heat)
Heatmap(heat,
        cluster_rows = T,
        cluster_columns = T, 
        border = T, 
        column_km = 2, 
        row_km = 4,
        rect_gp = gpar(col = "darkgray", lwd = 1),
        name = "scaled_expression", 
        show_row_dend = F, 
        show_column_dend = F)

# Allart's email 05.30.2022 ####
# summarise the overall gene expressions for the attached gene list into a single violin plot comparing 
# anti-VB positive (= TRBV6+ cells), 
# anti-CD3 (cluster 0 or cluster 4) and 
# anti-VB negative cells for 
# CD4 and a separate one for CD8
# merge all genes together and make two violin plots, 3 violins each
genes2use <- as.character(data.frame(readxl::read_xlsx("20220520 Trm genes for plotting fig S6E.xlsx"))[,1])

rownames(sc_data_no_int@assays$RNA@data)[grep("LY6", rownames(sc_data_no_int@assays$RNA@data), ignore.case = T)]
genes2use[-which(genes2use %in% rownames(sc_data_no_int@assays$RNA@data))]
colnames(sc_data_no_int@meta.data)

# CD4
groups_4_violin_cd4 <- list(AntiVB_pos_CD4 = which(sc_data_no_int$TRBV_for_dge == "test" & sc_data_no_int$CD8_CD4 == "CD4" & sc_data_no_int$AHH == "AHH04"),
                            AntiVB_neg_CD4 = which(sc_data_no_int$TRBV_for_dge == "control" & sc_data_no_int$CD8_CD4 == "CD4" & sc_data_no_int$AHH == "AHH04"),
                            AntiCD3_c0_CD4 = which(sc_data_no_int$AHH == "AHH03" & sc_data_no_int$seurat_clusters == "0" & sc_data_no_int$CD8_CD4 == "CD4")
)
sc_data_no_int$CD4_new_list <- "NA"
sc_data_no_int$CD4_new_list <- ifelse(1:length(sc_data_no_int$CD4_new_list) %in% groups_4_violin_cd4[[3]], names(groups_4_violin_cd4)[3],
                                      ifelse(1:length(sc_data_no_int$CD4_new_list) %in% groups_4_violin_cd4[[1]], names(groups_4_violin_cd4)[1],
                                             ifelse(1:length(sc_data_no_int$CD4_new_list) %in% groups_4_violin_cd4[[2]], names(groups_4_violin_cd4)[2],"")))
table(sc_data_no_int$CD4_new_list)
temp_row_ind <- which(rownames(sc_data_no_int@assays$RNA@data) %in% genes2use)
temp_col_ind <- which(sc_data_no_int$CD4_new_list != "")
df_data4vln <- data.frame(t(sc_data_no_int@assays$RNA@data[temp_row_ind, temp_col_ind]),
                          Group = sc_data_no_int$CD4_new_list[temp_col_ind])
df_data4vln <- data.frame(Log2Counts = apply(log2(df_data4vln[,1:32]+1), 1, mean), 
                          Group = df_data4vln$Group)
head(df_data4vln)
temp_vln <- df_data4vln


ggplot(temp_vln,
       aes(x = factor(Group, levels = c("AntiVB_neg_CD4", "AntiCD3_c0_CD4", "AntiVB_pos_CD4")),
           y = Log2Counts,
           fill = Group)) +
  geom_violin(scale = "width") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=12)) +
  xlab("") 
stat_compare_means(method = "t.test", 
                   comparisons = list(c("AntiVB_neg_CD4", "AntiCD3_c0_CD4"),
                                      c("AntiCD3_c0_CD4", "AntiVB_pos_CD4"),
                                      c("AntiVB_neg_CD4", "AntiVB_pos_CD4")), 
                   label = "p.format",
                   vjust = 0.1,
                   step.increase = 0.15)

ggplot(temp_vln,
       aes(x=Log2Counts, 
           col=Group)) +
  stat_ecdf(geom = "step", size=1.25) +
  theme_bw() +
  labs(col = "Group") +
  ylab("cumulative frequency") +
  ggtitle("")

v_AntiVB_neg_CD4 <- ecdf(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_neg_CD4")])
f_v_AntiVB_neg_CD4 <- v_AntiVB_neg_CD4(sort(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_neg_CD4")]))
v_AntiCD3_c0_CD4 <- ecdf(temp_vln$Log2Counts[which(temp_vln$Group == "AntiCD3_c0_CD4")])
f_v_AntiCD3_c0_CD4 <- v_AntiCD3_c0_CD4(sort(temp_vln$Log2Counts[which(temp_vln$Group == "AntiCD3_c0_CD4")]))
v_AntiVB_pos_CD4 <- ecdf(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_pos_CD4")])
f_v_AntiVB_pos_CD4 <- v_AntiVB_pos_CD4(sort(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_pos_CD4")]))
ks.test(knots(v_AntiVB_neg_CD4), knots(v_AntiCD3_c0_CD4))
twosamples::cvm_test(knots(v_AntiVB_neg_CD4), knots(v_AntiCD3_c0_CD4))
ks.test(knots(v_AntiVB_neg_CD4), knots(v_AntiVB_pos_CD4))
twosamples::cvm_test(knots(v_AntiVB_neg_CD4), knots(v_AntiVB_pos_CD4))
ks.test(knots(v_AntiCD3_c0_CD4), knots(v_AntiVB_pos_CD4))
twosamples::cvm_test(knots(v_AntiCD3_c0_CD4), knots(v_AntiVB_pos_CD4))


# CD8
groups_4_violin_cd8 <- list(AntiVB_pos_CD8 = which(sc_data_no_int$TRBV_for_dge == "test" & sc_data_no_int$CD8_CD4 == "CD8" & sc_data_no_int$AHH == "AHH04"),
                            AntiVB_neg_CD8 = which(sc_data_no_int$TRBV_for_dge == "control" & sc_data_no_int$CD8_CD4 == "CD8" & sc_data_no_int$AHH == "AHH04"),
                            AntiCD3_c4_CD8 = which(sc_data_no_int$AHH == "AHH03" & sc_data_no_int$seurat_clusters == "4" & sc_data_no_int$CD8_CD4 == "CD8")
)
sc_data_no_int$CD8_new_list <- "NA"
sc_data_no_int$CD8_new_list <- ifelse(1:length(sc_data_no_int$CD8_new_list) %in% groups_4_violin_cd8[[3]], names(groups_4_violin_cd8)[3],
                                      ifelse(1:length(sc_data_no_int$CD8_new_list) %in% groups_4_violin_cd8[[1]], names(groups_4_violin_cd8)[1],
                                             ifelse(1:length(sc_data_no_int$CD8_new_list) %in% groups_4_violin_cd8[[2]], names(groups_4_violin_cd8)[2],"")))
table(sc_data_no_int$CD8_new_list)
temp_row_ind <- which(rownames(sc_data_no_int@assays$RNA@data) %in% genes2use)
temp_col_ind <- which(sc_data_no_int$CD8_new_list != "")
df_data4vln <- data.frame(t(sc_data_no_int@assays$RNA@data[temp_row_ind, temp_col_ind]),
                          Group = sc_data_no_int$CD8_new_list[temp_col_ind])
df_data4vln <- data.frame(Log2Counts = apply(log2(df_data4vln[,1:32]+1), 1, mean), 
                          Group = df_data4vln$Group)
head(df_data4vln)
temp_vln <- df_data4vln

ggplot(temp_vln,
       aes(x = factor(Group, levels = c("AntiVB_neg_CD8", "AntiCD3_c4_CD8", "AntiVB_pos_CD8")),
           y = Log2Counts,
           fill = Group)) +
  geom_violin(scale = "width") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=12)) +
  xlab("") 
stat_compare_means(method = "t.test", 
                   comparisons = list(c("AntiVB_neg_CD8", "AntiCD3_c4_CD8"),
                                      c("AntiCD3_c4_CD8", "AntiVB_pos_CD8"),
                                      c("AntiVB_neg_CD8", "AntiVB_pos_CD8")), 
                   label = "p.format",
                   vjust = 0.1,
                   step.increase = 0.15)

ggplot(temp_vln,
       aes(x=Log2Counts, 
           col=Group)) +
  stat_ecdf(geom = "step", size=1.25) +
  theme_bw() +
  labs(col = "Group") +
  ylab("cumulative frequency") +
  ggtitle("")

v_AntiVB_neg_CD8 <- ecdf(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_neg_CD8")])
f_v_AntiVB_neg_CD8 <- v_AntiVB_neg_CD8(sort(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_neg_CD8")]))
v_AntiCD3_c4_CD8 <- ecdf(temp_vln$Log2Counts[which(temp_vln$Group == "AntiCD3_c4_CD8")])
f_v_AntiCD3_c4_CD8 <- v_AntiCD3_c4_CD8(sort(temp_vln$Log2Counts[which(temp_vln$Group == "AntiCD3_c4_CD8")]))
v_AntiVB_pos_CD8 <- ecdf(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_pos_CD8")])
f_v_AntiVB_pos_CD8 <- v_AntiVB_pos_CD8(sort(temp_vln$Log2Counts[which(temp_vln$Group == "AntiVB_pos_CD8")]))
ks.test(knots(v_AntiVB_neg_CD8), knots(v_AntiCD3_c4_CD8))
twosamples::cvm_test(knots(v_AntiVB_neg_CD8), knots(v_AntiCD3_c4_CD8))
ks.test(knots(v_AntiVB_neg_CD8), knots(v_AntiVB_pos_CD8))
twosamples::cvm_test(knots(v_AntiVB_neg_CD8), knots(v_AntiVB_pos_CD8))
ks.test(knots(v_AntiCD3_c4_CD8), knots(v_AntiVB_pos_CD8))
twosamples::cvm_test(knots(v_AntiCD3_c4_CD8), knots(v_AntiVB_pos_CD8))


# prep STAT3 and 4 plots in cluster 4 and 7 as per meeting 09012022 ####
table(sc_data_no_int$Protein_snn_res.0.25)
DimPlot(sc_data_no_int,
        reduction = "prot_UMAP", 
        group.by = "Protein_snn_res.0.25")

ind_4_7 <- which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))

FeatureScatter(sc_data_no_int, 
               jitter = T,
               feature1 = "STAT3",
               feature2 = "STAT4",
               cells = ind_4_7)


# data need to be imputed
sc_data_no_int <- magic(sc_data_no_int, 
                        seed = 1234,
                        verbose = T,
                        npca = 40,
                        solver = "approximate")

DefaultAssay(sc_data_no_int) <- "MAGIC_SCT"


# calculate stats
df_4_cor <- data.frame(cbind(t(sc_data_no_int@assays$MAGIC_SCT@data[c("STAT3", "STAT4", "BCL6", "TCF7", "GZMB"),ind_4_7]),
                             t(sc_data_no_int@assays$SCT@data[c("STAT3", "STAT4", "BCL6", "TCF7", "GZMB"),ind_4_7])))
colnames(df_4_cor) <- c("STAT3_MAGIC", "STAT4_MAGIC","BCL6_MAGIC", "TCF7_MAGIC", "GZMB_MAGIC",
                        "STAT3_SCT","STAT4_SCT","BCL6_SCT", "TCF7_SCT", "GZMB_SCT")
df_4_cor <- df_4_cor %>%
  mutate(UID = rownames(.)) %>%
  left_join(sc_data_no_int@meta.data %>% mutate(UID = rownames(.)), by = "UID")

head(df_4_cor)


# plot STAT3 vs. STAT4
ggplot(df_4_cor,
       aes(x = STAT3_MAGIC,
           y = STAT4_MAGIC,
           col = AHH)) +
  geom_point(size = 0.5) +
  theme_classic2() +
  facet_wrap(Protein_snn_res.0.25 ~ AHH) +
  scale_color_manual(values = c("AHH01" = "#1c61a5",
                                "AHH03" = "#fc6812",
                                "AHH04" = "#269320")) +
  geom_smooth(method = "glm", se = T, show.legend = F) +
  stat_cor(label.sep = "\n",
           method = "spearman",
           p.accuracy = 0.001, 
           r.accuracy = 0.01,
           label.y = 0.44,
           label.x = 0.7,
           show.legend = F,
           cor.coef.name = "rho", 
           col = "black") +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5)))


# plot STAT3_MAGIC vs. STAT3_SCT
ggplot(df_4_cor,
       aes(x = STAT3_MAGIC,
           y = STAT3_SCT,
           col = AHH)) +
  geom_point(size = 0.5) +
  theme_classic2() +
  scale_color_manual(values = c("AHH01" = "#1c61a5",
                                "AHH03" = "#fc6812",
                                "AHH04" = "#269320")) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5)))


# plot STAT4_MAGIC vs. STAT4_SCT
ggplot(df_4_cor,
       aes(x = STAT4_MAGIC,
           y = STAT4_SCT,
           col = AHH)) +
  geom_point(size = 0.5) +
  theme_classic2() +
  scale_color_manual(values = c("AHH01" = "#1c61a5",
                                "AHH03" = "#fc6812",
                                "AHH04" = "#269320")) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5)))


# plot BCL6_MAGIC vs. BCL6_SCT
ggplot(df_4_cor,
       aes(x = GZMB_MAGIC,
           y = GZMB_SCT,
           col = AHH)) +
  geom_point(size = 0.5) +
  theme_classic2() +
  scale_color_manual(values = c("AHH01" = "#1c61a5",
                                "AHH03" = "#fc6812",
                                "AHH04" = "#269320")) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5)))


# address ALlart/Adrian questions 10042022 ####
# plot STAT3 vs. STAT4
ggplot(df_4_cor,
       aes(x = STAT3_MAGIC,
           y = STAT4_MAGIC,
           col = AHH)) +
  geom_point(size = 0.5) +
  theme_classic2() +
  facet_wrap(Protein_snn_res.0.25 ~ AHH) +
  scale_color_manual(values = c("AHH01" = "#1c61a5",
                                "AHH03" = "#fc6812",
                                "AHH04" = "#269320")) +
  guides(color=guide_legend(ncol=1, 
                            override.aes = list(size=5))) +
  geom_vline(xintercept = 0.6) +
  geom_hline(yintercept = 0.3)

ggplot(sc_data_no_int@assays$SCT@data["STAT3", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))],
       aes(x = STAT3)) +
  geom_density2d()

hist(as.numeric(sc_data_no_int@assays$SCT@data["STAT3", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]))
hist(as.numeric(sc_data_no_int@assays$MAGIC_SCT@data["GZMB", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]))

par(mfrow=c(4,1))
hist(as.numeric(sc_data_no_int@assays$MAGIC_SCT@data["STAT4", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]),
     xlim = c(0,7))
hist(as.numeric(sc_data_no_int@assays$MAGIC_SCT@data["STAT3", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]),
     xlim = c(0,7))
hist(as.numeric(sc_data_no_int@assays$MAGIC_SCT@data["GZMB", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]),
     xlim = c(0,7))
hist(as.numeric(sc_data_no_int@assays$MAGIC_SCT@data["BCL6", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]),
     xlim = c(0,7))
hist(as.numeric(sc_data_no_int@assays$MAGIC_SCT@data["GAPDH", which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))]),
     xlim = c(0,7))


# test genes for their overall expression, 10082022 ####
genes2test <- c("GZMB","STAT3", "STAT4", "KLRG1", "TCF7", "LEF1", "KLF2", "BCL6", "PRDM1", "MKI67", "LAG3", "TBX21", "IRF4", "IRF1", "TOX", "PDCD1", "RUNX3", "BATF3", "IFNG")
row_ind <- which(rownames(sc_data_no_int@assays$SCT@data) %in% genes2test)
col_meta_ind <- which(colnames(sc_data_no_int@meta.data) %in% c("Protein_snn_res.0.25", "AHH", "TRBV_for_dge"))
genes2test[-which(genes2test %in% rownames(sc_data_no_int@assays$SCT@data))]
df_genes <- data.frame(sc_data_no_int@assays$SCT@data[row_ind,]) %>%
  dplyr::mutate(GENE = rownames(.)) %>%
  tidyr::pivot_longer(!GENE,
                      names_to = "UID", 
                      values_to = "Expression") %>%
  dplyr::left_join(sc_data_no_int@meta.data[,col_meta_ind] %>% dplyr::mutate(UID = gsub("-","\\.",rownames(.))), by = "UID") %>%
  dplyr::rename(CLUSTER = Protein_snn_res.0.25) %>%
  dplyr::filter(CLUSTER %in% c(4,7))
head(df_genes)
ggplot(df_genes,
       aes(x = Expression,
           y = GENE)) +
  geom_density_ridges2(panel_scaling = T) + 
  theme_classic() +
  facet_wrap(~ CLUSTER, scales = "free") +
  xlab("Not Imputed Expression")

plot(sc_data_no_int@assays$MAGIC_SCT@data["TOX",],
     sc_data_no_int@assays$MAGIC_SCT@data["CD3E",])


# use new list of genes for correlation analysis 10112022 ####
##### also used for call with KCL on 02072022 ####
# [a] cluster 4 cells
# [b] cluster 7 cells
# [c] cluster 0 cells
# [d] cluster 5 cells
# [e] cluster 6 cells

gene.pairs <- list(
  c("TCF7", "LEF1"),
  c("TCF7","BCL6"),
  c("BCL6","PRDM1"),
  c("MKI67","LAG3"),
  c("LAG3","TBX21"),
  c("IRF4","IRF1"),
  c("IRF4","RUNX3"),
  c("IRF4","BATF3")
)
ind_clusters <- list("4" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 4),
                     "7" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 7),
                     "0" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 0),
                     "5" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 5),
                     "6" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 6))

# calcualte MAGIC ####
sc_data_no_int <- magic(sc_data_no_int, 
                        seed = 1234,
                        verbose = T,
                        npca = 40,
                        solver = "approximate")

DefaultAssay(sc_data_no_int) <- "MAGIC_SCT"


for(ii in 1:length(gene.pairs)){
  ind_genes <- c(which(rownames(sc_data_no_int@assays$MAGIC_SCT@data) %in% unlist(gene.pairs[[ii]][[1]])),
                 which(rownames(sc_data_no_int@assays$MAGIC_SCT@data) %in% unlist(gene.pairs[[ii]][[2]])))
  
  for(kk in 1:length(ind_clusters)){
    temp.sc <- data.frame(t(sc_data_no_int@assays$MAGIC_SCT@data[ind_genes, ind_clusters[[kk]]]))
    colnames(temp.sc) <- c("g1", "g2")
    
    if(kk == 1){
      xmax <- max(temp.sc$g1, na.rm = T)
      ymax <- max(temp.sc$g2, na.rm = T)
    }
    
    temp.gg <- ggplot(temp.sc,
                      aes(x = g1,
                          y = g2)) +
      geom_point(size = 0.75) +
      theme_classic() + 
      xlab(gene.pairs[[ii]][[1]]) + ylab(gene.pairs[[ii]][[2]]) +
      xlim(0,xmax*1.1) + ylim(0,ymax*1.1) +
      ggtitle(paste("Cluster", names(ind_clusters)[kk]))
    
    if(kk == 1){
      gg.to.plot <- list(ggplotify::as.grob(temp.gg))
    } else {
      gg.to.plot[[kk]] <- ggplotify::as.grob(temp.gg)
    }
  }
  gg <- ggarrange(plotlist = gg.to.plot, 
                  nrow = 1,
                  ncol = 5)
  print(gg)
}


# prep spider plots - 10242022 ####
library(fmsb)


# From the literature / databases
# Naïve CD8 T cells
# Naïve CD4 T cells
# CD8 TRM
# CD8 TCM
# CD4 TCM
# CD4 TEM
# CD8 TEM
# CD4 TEMRA
# CD8 TEMRA
# CD8 T SLEC
# CD8 T EX
# CD8 PEX (precursors of TEX)

# load data Science paper
# T cell subsets : Naïve, TEM, TCM, TRM, TEMRA, TEX 
meta.data.zeng <- fread("../Zeng_science_panCancer/GSE156728_metadata.txt.gz")

files <- list.files("../Zeng_science_panCancer", 
                    full.names = T, 
                    recursive = T,
                    pattern = "*.txt.gz")[-5]

genes4adar <- c("IL7R","IFNG", "CCR7", "CX3CR1", "PDCD1", "TCF7", "TOX", "RUNX3", "IRF4", "TBX21", "KLRG1", "KLF2", "IRF1", "IL7RA", "SELL", "IL2RA")

#' c0 - c5 CD4 cells
#' c4 - c7 CD8 cells
#' divide plot by cell type

for(ii in 1:length(files)){
  
  # skip the metadata file
  if(files[ii] == "../Zeng_science_panCancer/GSE156728_metadata.txt.gz")
    next
  
  # get tumor type
  message("")
  message(" -- -- -- ")
  message(paste("working with", files[ii]))
  
  if(grepl("GSE156728", files[ii])){
    temp.tumor.type <- gsub("\\.\\./Zeng_science_panCancer/GSE156728_","",files[ii])
  } 
  if(grepl("GSM4743199", files[ii])){
    temp.tumor.type <- gsub("\\.\\./Zeng_science_panCancer/GSM4743199_","",files[ii])
  } 
  if(grepl("GSM4743231", files[ii])){
    temp.tumor.type <- gsub("\\.\\./Zeng_science_panCancer/GSM4743231_","",files[ii])
  }
  temp.tumor.type <- unlist(strsplit(temp.tumor.type, "_")[[1]][1])
  
  temp.cd4.cd8 <- ifelse(grepl("CD4", files[ii]), "CD4","CD8")
  
  pdf(paste0("Radial_plots_IFNG_CD4CD8_",temp.tumor.type,"_",temp.cd4.cd8,".pdf"), width = 9.5, height = 6)
  
  # select the cell types that they need
  # here we assume Tm as TCM and Tn as Naive
  temp.cells.to.keep <- c("Tex","Temra","Trm","Tm","Tem","Tn")
  
  # get metadata
  # cancerType has the tumor type data
  # we assume that loc includes normal (N) vs. tumor (T)
  temp.meta <- meta.data.zeng %>%
    dplyr::filter(cancerType == temp.tumor.type,
                  grepl(paste(temp.cells.to.keep, collapse = ".|."), meta.cluster),
                  loc == "T") %>%
    dplyr::mutate(CellType = sapply(meta.cluster, function(xx) strsplit(xx, "\\.")[[1]][3]))
  
  # read gene expression data
  temp.data <- data.frame(fread(files[ii]))
  rownames(temp.data) <- temp.data$V1
  temp.data$V1 <- NULL
  temp.ind.cols <- c(1,which(colnames(temp.data) %in% temp.meta$cellID))
  temp.ind.rows <- which(rownames(temp.data) %in% genes4adar)
  # filter the data keeping only those annotated in meta.data
  temp.data <- data.frame(t(temp.data[temp.ind.rows,temp.ind.cols])) %>%
    dplyr::mutate(cellID = rownames(.)) %>%
    dplyr::inner_join(temp.meta, by = "cellID")
  
  # take the mean log2 value of the cell type
  temp.data.mean <- aggregate(temp.data[,c(1:15)], by=list(temp.data$CellType), FUN=mean)
  rownames(temp.data.mean) <- temp.data.mean$Group.1
  temp.data.mean$Group.1 <- NULL
  
  # merge the data from clusters 4,7,5,0
  if(temp.cd4.cd8 == "CD8"){
    temp.sc.mean <- t(AverageExpression(subset(sc_data_no_int,
                                               cells = which(sc_data_no_int$Protein_snn_res.0.25 %in% c(4,7))), 
                                        assays = "RNA",
                                        group.by = "Protein_snn_res.0.25",
                                        features = genes4adar, 
                                        slot = "counts")$RNA)
  } 
  if(temp.cd4.cd8 == "CD4"){
    temp.sc.mean <- t(AverageExpression(subset(sc_data_no_int,
                                               cells = which(sc_data_no_int$Protein_snn_res.0.25 %in% c(0,5))), 
                                        assays = "RNA",
                                        group.by = "Protein_snn_res.0.25",
                                        features = genes4adar, 
                                        slot = "counts")$RNA)
  } 
  
  # prep data for radial plot by cell type
 
  temp.avg.merged <- rbind(t(apply(log2(temp.sc.mean+1), 1, function(x) (x-min(x))/(max(x)-min(x)) )),
                           t(apply(log2(temp.data.mean+1), 1, function(x) (x-min(x))/(max(x)-min(x)) )))
  temp.max <- max(temp.avg.merged, na.rm = T)
  temp.avg.merged <- rbind(rep(temp.max, ncol(temp.avg.merged)),
                           rep(0, ncol(temp.avg.merged)),
                           temp.avg.merged)
  
  # plot radial plot by cell type
  for(kk in temp.cells.to.keep){
    message(paste(" >> plotting",kk))
    
    ind.cell.type <- which(rownames(temp.avg.merged) == kk)
    
    # if cell type is not there skip it
    if(length(ind.cell.type) == 0){
      next
    } else {
      
      colors_border <- c25
      colors_in <- scales::alpha(c25, alpha = 0.2)
      
      temp.radial.df <- temp.avg.merged[c(1:4,ind.cell.type),]
      
      # create title
      temp.title <- paste(paste(temp.cd4.cd8,kk), 
                          paste("Tumor =",temp.tumor.type), 
                          sep = "\n")
      
      # do the actual plot
      radarchart(data.frame(temp.radial.df),
                 axistype=1 , 
                 #custom polygon
                 pcol=colors_border , pfcol=colors_in , 
                 plwd=1 , plty=1,
                 #custom the grid
                 cglcol="grey", cglty=1, axislabcol="grey", 
                 caxislabels=seq(0,temp.max,ncol(temp.avg.merged)),
                 cglwd=0.8,
                 #custom labels
                 vlcex=1.4,
                 title = temp.title
      )
      legend(x = -2,
             y = 1, 
             legend = rownames(temp.radial.df[-c(1,2),]), 
             bty = "n", 
             pch=20 , 
             col=colors_in, 
             text.col = "black", 
             cex=1.2, 
             pt.cex=3)
    }
  }
  dev.off()
  message("")
  message("")
}



# load data - Andreatta paper
data <- readxl::read_xlsx("../Andreatta_Nature_scNASeq_MC38_mouse/41467_2021_23324_MOESM4_ESM.xlsx", 
                          sheet = 2,
                          skip = 2)


data <- data.frame(t(data %>%
                       dplyr::filter(Gene %in% genes4adar)))
colnames(data) <- data[1,]
data <- data[-1,]
rws <- rownames(data)
data <- apply(data, 2, as.numeric)
rownames(data) <- rws
max(data)


# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- data.frame(rbind(rep(2.25,ncol(data)) , rep(0,ncol(data)) , data))

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart( data[1:4,]  , axistype=1 , 
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,2.25,14), cglwd=0.8,
            vlcex=1.4
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),])[1:2], bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)

# using gene pairs above make MAGIC vs. SCT vs. counts ####
genes.scatter <- c(unique(unlist(gene.pairs)),"CD3D")
ind_clusters <- list("4" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 4),
                     "7" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 7),
                     "0" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 0),
                     "5" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 5),
                     "6" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 6))

for(ii in 12:length(genes.scatter)){
  for(kk in 1:length(ind_clusters)){
    
    # take data for magic
    x.magic <- sc_data_no_int@assays$MAGIC_SCT@data[genes.scatter[ii],
                                                    ind_clusters[[kk]]]
    # take data for SCT
    y.sct <- sc_data_no_int@assays$SCT@data[genes.scatter[ii],
                                            ind_clusters[[kk]]]
    
    # compile data
    if(kk == 1){
      df.res <- data.frame(MAGIC = x.magic,
                           SCT = y.sct,
                           CLUSTER = names(ind_clusters)[kk],
                           GENE = genes.scatter[ii])
    } else {
      df.res <- rbind(df.res,
                      data.frame(MAGIC = x.magic,
                                 SCT = y.sct,
                                 CLUSTER = names(ind_clusters)[kk],
                                 GENE = genes.scatter[ii]))
    }
  }
  
  gg <- ggplot(df.res,
               aes(x = MAGIC,
                   y = SCT)) +
    geom_point(size = 0.5) +
    theme_bw() +
    ggtitle(genes.scatter[ii]) +
    facet_wrap(~ CLUSTER, ncol = 5)
  print(gg)
}


# make violin plos as requested 02.17.2022 ####
genes4violin <- c("IL7R", "KLRG1", "TCF7", "IL2RA", "GZMB", "MKI67", "PDCD1", "HAVCR2", "LAG3", "TBX21", "KLF2", "PRDM1", "BCL6", "IRF4")
ind_clusters <- list("4" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 4),
                     "7" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 7),
                     "0" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 0),
                     "5" = which(sc_data_no_int@meta.data$Protein_snn_res.0.25 == 5))

temp.df <- sc_data_no_int@meta.data[as.numeric(unlist(ind_clusters)), ]
temp.genes.magic <- t(sc_data_no_int@assays$MAGIC_SCT@data[genes4violin, as.numeric(unlist(ind_clusters))])
temp.df.magic <- merge(temp.genes.magic, temp.df, by = "row.names")

temp.df.magic %>%
  dplyr::select(genes4violin, AHH, Protein_snn_res.0.25,TRBV_for_dge, CD8_CD4) %>%
  dplyr::rename(Cluster = Protein_snn_res.0.25) %>% 
  pivot_longer(-c(Cluster, AHH, TRBV_for_dge,CD8_CD4)) %>%
  dplyr::mutate(Cluster = as.character(Cluster)) -> temp.df.long.magic

temp.genes.sct <- t(sc_data_no_int@assays$SCT@data[genes4violin, as.numeric(unlist(ind_clusters))])
temp.df.sct <- merge(temp.genes.sct, temp.df, by = "row.names")

temp.df.sct %>%
  dplyr::select(genes4violin, AHH, Protein_snn_res.0.25,TRBV_for_dge,CD8_CD4) %>%
  dplyr::rename(Cluster = Protein_snn_res.0.25) %>% 
  pivot_longer(-c(Cluster, AHH, TRBV_for_dge,CD8_CD4)) %>%
  dplyr::mutate(Cluster = as.character(Cluster)) -> temp.df.long.sct


for(ii in 1:length(genes4violin)){
  temp.violin.magic <- temp.df.long.magic %>%
    dplyr::filter(name == genes4violin[ii])
  
  pdf(paste0("pdf_figures_02272023/violin_",genes4violin[ii],"_CD4_C5_C0.pdf"), width = 3, height = 2.7)
  temp.stats.0.5.gene <- diff_04vs03_cd4_c5_c0 %>%
    dplyr::filter(Gene == genes4violin[ii])
  adjp <- round(temp.stats.0.5.gene$p_val_adj, 4)
  temp.title <- paste(genes4violin[ii],
                      paste0("Log2FC = ", round(temp.stats.0.5.gene$avg_log2FC,2)),
                      ifelse(adjp == 0, "Adj.P.value < 1e-4", paste0("Adj.P.value = ", adjp)),
                      sep = "\n")
  
  temp.violin.0.5 <- temp.violin.magic %>%
    dplyr::filter(Cluster %in% c("0", "5"),
                  CD8_CD4 == "CD4")
  gg<- ggplot(temp.violin.0.5,
              aes(x = Cluster,
                  y = value,
                  fill = factor(Cluster, levels = c("0", "5")))) +
    geom_violin(scale = "width") +
    theme_classic() +
    labs(fill = "Cluster") + 
    ggtitle(temp.title) + ylab("Imputed Expression")
  print(gg)
  
  temp.violin.sct <- temp.df.long.sct %>%
    dplyr::filter(name == genes4violin[ii])
  temp.violin.0.5 <- temp.violin.sct %>%
    dplyr::filter(Cluster %in% c("0", "5"),
                  CD8_CD4 == "CD4")
  gg <- ggplot(temp.violin.0.5,
               aes(x = Cluster,
                   y = value,
                   fill = factor(Cluster, levels = c("0", "5")))) +
    geom_violin(scale = "width") +
    theme_classic() +
    labs(fill = "Cluster") + 
    ggtitle(temp.title) + ylab("SCT Expression")
  print(gg)
  dev.off()
  
  # same for 7 and 4
  pdf(paste0("pdf_figures_02272023/violin_",genes4violin[ii],"_CD8_C7_C4.pdf"), width = 3, height = 2.7)
  temp.stats.7.4.gene <- diff_04vs03_cd8_c7_c4 %>%
    dplyr::filter(Gene == genes4violin[ii])
  adjp <- round(temp.stats.7.4.gene$p_val_adj, 4)
  temp.title <- paste(genes4violin[ii],
                      paste0("Log2FC = ", round(temp.stats.7.4.gene$avg_log2FC,2)),
                      ifelse(adjp == 0, "Adj.P.value < 1e-4", paste0("Adj.P.value = ", adjp)),
                      sep = "\n")
  
  temp.violin.7.4 <- temp.violin.magic %>%
    dplyr::filter(Cluster %in% c("4", "7"),
                  CD8_CD4 == "CD8")
  gg<- ggplot(temp.violin.7.4,
              aes(x = Cluster,
                  y = value,
                  fill = factor(Cluster, levels = c("4", "7")))) +
    geom_violin(scale = "width") +
    theme_classic() +
    labs(fill = "Cluster") + 
    ggtitle(temp.title) + ylab("Imputed Expression")
  print(gg)
  
  temp.violin.sct <- temp.df.long.sct %>%
    dplyr::filter(name == genes4violin[ii])
  temp.violin.7.4 <- temp.violin.sct %>%
    dplyr::filter(Cluster %in% c("4", "7"),
                  CD8_CD4 == "CD8")
  gg <- ggplot(temp.violin.7.4,
               aes(x = Cluster,
                   y = value,
                   fill = factor(Cluster, levels = c("4", "7")))) +
    geom_violin(scale = "width") +
    theme_classic() +
    labs(fill = "Cluster") + 
    ggtitle(temp.title) + ylab("SCT Expression")
  print(gg)
  dev.off()
  
  message(paste("Done for", genes4violin[ii]))
  
}

#' make venn diagrams for TFs between:
#'  TRBV6/10 vs other  (dge_trbv_6s10_vs_other)
#'  CD25+ vs CD25-  (DGE_CD8_CD25_PosVsNeg_01272022.txt and DGE_CD4_CD25_PosVsNeg_01272022.txt)
#' for both CD4 and CD8 
dge_trbv_6s10_vs_other_cd8 <- read.table("dge_trbv_6s10_vs_other_cd8_03312022.txt", sep="\t",header = T) %>%
  filter(Gene %in% all_TF,
         p_val_adj <0.05)
dge_cd25_posneg_04_cd8 <- fread("DGE_CD8_CD25_PosVsNeg_01272022.txt")
venn(list(CD8_CD25pos_vs_CD25neg_UP = dge_cd25_posneg_04_cd8 %>%
            filter(Gene %in% all_TF,
                   p_val_adj <0.05,
                   Direction == "UP") %>%
            pull(Gene),
          CD8_TRBV6.10_vs_other_UP = dge_trbv_6s10_vs_other_cd8 %>%
            filter(Direction == "UP") %>%
            pull(Gene)))
phyper(q = 26,
       m = 32, 
       n = length(all_TF) - 32,
       k = 51,
       lower.tail=FALSE
)

venn(list(CD8_CD25pos_vs_CD25neg_DOWN = dge_cd25_posneg_04_cd8 %>%
            filter(Gene %in% all_TF,
                   p_val_adj <0.05,
                   Direction == "DOWN") %>%
            pull(Gene),
          CD8_TRBV6.10_vs_other_DOWN = dge_trbv_6s10_vs_other_cd8 %>%
            filter(Direction == "DOWN") %>%
            pull(Gene)))
phyper(q = 8,
       m = 12, 
       n = length(all_TF) - 12,
       k = 14,
       lower.tail=FALSE
)

venn(list(CD8_CD25pos_vs_CD25neg = dge_cd25_posneg_04_cd8 %>%
            filter(Gene %in% all_TF,
                   p_val_adj <0.05) %>%
            pull(Gene),
          CD8_TRBV6.10_vs_other = dge_trbv_6s10_vs_other_cd8 %>%
            pull(Gene)))
phyper(q = 34,
       m = 44, 
       n = length(all_TF) - 44,
       k = 41,
       lower.tail=FALSE
)

dge_cd25_posneg_04_cd4 <- fread("DGE_CD4_CD25_PosVsNeg_01272022.txt") %>%
  filter(Gene %in% all_TF,
         p_val_adj <0.05)
venn(list(CD4_CD25pos_vs_CD25neg_UP = dge_cd25_posneg_04_cd4 %>%
            filter(Gene %in% all_TF,
                   p_val_adj <0.05,
                   Direction == "UP") %>%
            pull(Gene),
          CD4_TRBV6.10_vs_other_UP = dge_trbv_6s10_vs_other_cd4 %>%
            filter(Direction == "UP") %>%
            pull(Gene)))
phyper(q = 33,
       m = 38, 
       n = length(all_TF) - 38,
       k = 46,
       lower.tail=FALSE
)

venn(list(CD4_CD25pos_vs_CD25neg_DOWN = dge_cd25_posneg_04_cd4 %>%
            filter(Gene %in% all_TF,
                   p_val_adj <0.05,
                   Direction == "DOWN") %>%
            pull(Gene),
          CD4_TRBV6.10_vs_other_DOWN = dge_trbv_6s10_vs_other_cd4 %>%
            filter(Direction == "DOWN") %>%
            pull(Gene)))
phyper(q = 13,
       m = 19, 
       n = length(all_TF) - 19,
       k = 15,
       lower.tail=FALSE
)

venn(list(CD4_CD25pos_vs_CD25neg = dge_cd25_posneg_04_cd4 %>%
            filter(Gene %in% all_TF,
                   p_val_adj <0.05) %>%
            pull(Gene),
          CD4_TRBV6.10_vs_other = dge_trbv_6s10_vs_other_cd4 %>%
            pull(Gene)))
phyper(q = 46,
       m = 57, 
       n = length(all_TF) - 57,
       k = 61,
       lower.tail=FALSE
)

dge_cd25_posneg_04_cd4



VlnPlot(subset(sc_data_no_int, cells = which(sc_data_no_int$Protein_snn_res.0.25 %in% c(0,4,5,7))), 
        features = "IL2RA",
        group.by = "Protein_snn_res.0.25",

        pt.size = 0) +
  facet_wrap(~ sc_data_no_int$CD8_CD4[which(sc_data_no_int$Protein_snn_res.0.25 %in% c(0,4,5,7))])

# make individual scatterplots
for(ii in 1:length(gene.pairs)){
  ind_genes <- c(which(rownames(sc_data_no_int@assays$MAGIC_SCT@data) %in% unlist(gene.pairs[[ii]][[1]])),
                 which(rownames(sc_data_no_int@assays$MAGIC_SCT@data) %in% unlist(gene.pairs[[ii]][[2]])))
  
  for(kk in 1:length(ind_clusters)){
    temp.sc <- data.frame(t(sc_data_no_int@assays$MAGIC_SCT@data[ind_genes, ind_clusters[[kk]]]))
    colnames(temp.sc) <- c("g1", "g2")
    
    if(kk == 1){
      xmax <- max(temp.sc$g1, na.rm = T)
      ymax <- max(temp.sc$g2, na.rm = T)
    }
    
    pdf(paste0("scatter_plot_individual_02272023/",
               unlist(gene.pairs[[ii]][[1]]),
               "_",
               unlist(gene.pairs[[ii]][[2]]),
               "_Cluster_",
               names(ind_clusters)[kk],
               ".pdf"),
        width = 3,
        height = 3)
    gg <- ggplot(temp.sc,
                 aes(x = g1,
                     y = g2)) +
      geom_point(size = 0.75) +
      theme_classic() + 
      xlab(gene.pairs[[ii]][[1]]) + ylab(gene.pairs[[ii]][[2]]) +
      xlim(0,xmax*1.1) + ylim(0,ymax*1.1) +
      ggtitle(paste("Cluster", names(ind_clusters)[kk]))
    
    print(gg)
    dev.off()
  }
}

# plot PD1 CITESeq vs PD1 gene - 03072023 ####
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

ggplot(sc_data_no_int@meta.data,
       aes(x = sc_data_no_int@assays$MAGIC_SCT@data["PDCD1",],
           y = sc_data_no_int@assays$Protein@data["ADT-CD279-TotalSeqC",])) +
  geom_point(size = 0.25) +
  theme_bw() +
  xlab("Transcript") + ylab("CITESeq")

FeaturePlot(sc_data_no_int,
            features = c("PDCD1", "ADT-CD279-TotalSeqC"),
            col = c("gray", "tomato"),
            min.cutoff = "q25")

ggplot(sc_data_no_int@meta.data,
       aes(x = log10(range01(sc_data_no_int@assays$MAGIC_SCT@data["PDCD1",])+1))) +
  stat_density(fill = scales::alpha("black", 0.2), col = "black") +
  stat_density(aes(x = log10(range01(sc_data_no_int@assays$Protein@data["ADT-CD279-TotalSeqC",])+1)), fill = scales::alpha("red", 0.5), col = "tomato") +
  theme_bw() +
  xlab("PD1 levels") + ylab("Abundance Distribution")






