library(Seurat)
library(SoupX)
library(DoubletFinder)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
library(monocle3)
library(tidyverse)
library(pheatmap)
library(velocyto.R)
library(org.Hs.eg.db)
library(dplyr)
source("~/spatial_transcriptomics/R/toolbox/seurat_functions.R")
source("~/spatial_transcriptomics/R/toolbox/entrez.R")
source("~/spatial_transcriptomics/R/toolbox/thyroid.R")
setwd("~/Desktop/organoids/")
object <- readRDS("HS_ORGANOIDSC1_seurat_mt26.Rds")

## PREPROCESSING - PARAMETERS
s <- "HS_ORGANOIDSC1"
nExp <- 3000
dimensionality <- 1:30
filterMT <- TRUE
nCores <- 8

## PREPROCESSING - STANDARD WITH MT26
## sc <- load10X("outs/")
## sc <- autoEstCont(sc)
## mat <- adjustCounts(sc)
## rho <- sc$fit$rhoEst
## doublet.object <- CreateSeuratObject(counts = mat, project = s)
## doublet.object[["percent.mt"]] <- PercentageFeatureSet(doublet.object, pattern = "^MT-")
## doublet.object <- SCTransform(doublet.object, vars.to.regress = "percent.mt")
## doublet.object <- FindVariableFeatures(doublet.object, selection.method = "vst", nfeatures = 3000)
## doublet.object <- RunPCA(doublet.object, features = VariableFeatures(object = doublet.object),seed.use = 1)
## doublet.object <- FindNeighbors(doublet.object, dims = dimensionality)
## doublet.object <- FindClusters(doublet.object, resolution = 0.8)
## doublet.object <- RunUMAP(doublet.object, dims = dimensionality, seed.use = 1)
## ## pK Identification (no ground-truth) -------------------------------------------------------
## sweep.res.list <- paramSweep_v3(doublet.object, PCs = dimensionality, sct = TRUE)
## sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
## bcmvn <- find.pK(sweep.stats)
## selected.pk <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),"pK"]))
## ## Homotypic Doublet Proportion Estimate -----------------------------------------------------
## annotations <- doublet.object@meta.data$SCT_snn_res.0.8
## homotypic.prop <- modelHomotypic(annotations)
## doublet.rates.df <- data.frame(multiplet_rate = c(0.004,0.008,0.016,0.023,0.031,0.039,
##                                                   0.046,0.054,0.061,0.069,0.076),
##                                cells_recovered = c(500,1000,2000,3000,4000,5000,6000,
##                                                    7000,8000,9000,10000))
## doublet.rate <- doublet.rates.df[doublet.rates.df$cells_recovered == nExp,"multiplet_rate"]
## n.exp.poi <- round(doublet.rate * length(doublet.object@meta.data$SCT_snn_res.0.8))
## n.exp.poi.adj <- round(n.exp.poi*(1-homotypic.prop))
## ## Run DoubletFinder -------------------------------------------------------------------------
## doublet.object <- doubletFinder_v3(doublet.object, PCs = dimensionality, pN = 0.25,
##                                    pK = selected.pk,nExp = n.exp.poi.adj,
##                                    reuse.pANN = FALSE, sct = TRUE)
## df.classif <- paste0("DF.classifications_0.25_", selected.pk, "_", n.exp.poi.adj)
## is.doublet <- doublet.object@meta.data[[df.classif]] == "Doublet"
## ## rm(doublet.object)
## gc()

## ## Save complete raw counts matrix along with QC metrics
## raw.object <- CreateSeuratObject(counts = mat, project = s)
## raw.object[["percent.mt"]] <- PercentageFeatureSet(raw.object, pattern = "^MT-")
## raw.object[["passed.qc"]] <- raw.object$percent.mt < 26 & raw.object$nFeature_RNA > 200 & !is.doublet
## cells.rm <- sum(!raw.object$passed.qc)
## cells.rm.mt <- sum(raw.object$percent.mt > 26)
## cells.rm.features <- sum(raw.object$nFeature_RNA < 200)
## cells.rm.doublet <- sum(is.doublet)
## ## rm(mat)
## gc()

## ## Generate full Seurat object with filtered counts matrix
## object <- CreateSeuratObject(counts = raw.object@assays$RNA@counts[,raw.object$passed.qc], project = s)
## object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
## object@misc$cells.removed <- c(cells.rm, cells.rm.mt, cells.rm.features, cells.rm.doublet)
## names(object@misc$cells.removed) <- c("total", "mt", "features", "doublets")
## object@misc$soupx.rho <- rho
## object <- SCTransform(object, vars.to.regress = "percent.mt")
## object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
## object <- RunPCA(object, features = VariableFeatures(object = object), seed.use = 1)
## object <- FindNeighbors(object, dims = dimensionality)
## object <- FindClusters(object, resolution = 0.1)
## object <- FindClusters(object, resolution = 0.5)
## object <- FindClusters(object, resolution = 1)
## object <- RunUMAP(object, dims = dimensionality, seed.use = 1)
## DimPlot(object)
## ggsave("figures/processing/dimplot_pre_ccycle_r.png")

## ## PREPROCESSING - CELL CYCLE REGRESSION
## s.genes <- cc.genes$s.genes
## g2m.genes <- cc.genes$g2m.genes
## RidgePlot(object, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
## ggsave("figures/processing/ftplot_ccycle.png")
## object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
## object <- SCTransform(object, assay = 'RNA', new.assay.name = 'SCT',
##                       vars.to.regress = c('percent.mt', 'S.Score', 'G2M.Score'))
## object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
## object <- RunPCA(object, features = VariableFeatures(object = object), seed.use = 1)
## object <- FindNeighbors(object, dims = 1:30)
## object <- FindClusters(object, resolution = 0.1)
## object <- FindClusters(object, resolution = 0.5)
## object <- FindClusters(object, resolution = 1)
## object <- RunUMAP(object, dims = 1:30, seed.use = 1)

## ANNOTATION
DimPlot(object,
        group.by = "SCT_snn_res.1",
        ## group.by = "SCT_snn_res.0.5",
        label = T)
## ggsave("figures/processing/dimplot_post_ccycle_r.png")

FeaturePlot(object, features = "NANOG")
VlnPlot(object, features = "nCount_RNA")
object <- SetIdent(object = object, value = object@meta.data$SCT_snn_res.1)
m <- FindAllMarkers(object, only.pos = T) 
object@active.ident <- object$SCT_snn_res.1
m[m$cluster == 10, c(2,5,6,7)][1:25,]
m[grep("^HIST", m$gene), c(2,5,6,7)]
PrintEntrez("HES1")
lapply(m[m$cluster == 11, 7][2:5], PrintEntrez)
grep("NKX2-1",rownames(object), value = T)

a <- 0:max(as.numeric(levels(m$cluster))) ; a <- as.character(a)
names(a)[a %in% c(10)] <- "Thyroid progenitors" #ES cells : HISTs, POU3F1, SOX4, SOX11, KLF4, MYCL,RASD1
names(a)[a %in% c(1)] <- "Mature thyrocytes" #PAX8, FOXE1, TSHR, TG, TPO
names(a)[a %in% c(2)] <- "Immature thyrocytes" #PAX8, HHEX, NKX2-1, FOXE1, TSHR, TPO
names(a)[a %in% c(5)] <- "Immature thyrocytes" #PAX8, HHEX, NKX2-1, TSHR, TG, TPO
names(a)[a %in% c(9)] <- "Mature thyrocytes" #PAX8, NKX2-1, FOXE1, TSHR, TG
names(a)[a %in% c(0)] <- "Thyroid progenitors" #HHEX, NKX2-1
names(a)[a %in% c(4)] <- "Thyroid progenitors" #FOXE1
names(a)[a %in% c(6)] <- "Fibroblasts" #PRRX1, TWIST1, PRRX2, PDGFRA, COLs, MMPs
names(a)[a %in% c(7)] <- "Cardiovascular cells" #MYLs x3, MYHs x2, ACTA2, ACTG2, CNN1, CALD1
names(a)[a %in% c(3)] <- "Airway cells" #SOX2-7-15, KRT19, KRT5, TP63
names(a)[a %in% c(8)] <- "Endoderm epithelial cells" #FOXA1-2-3, SOX2-6, GATA6, KRT19, AGR3, ADAM28
names(a)[a %in% c(11)] <- "Endoderm epithelial cells" #KRT18-19, AGR3, ADAM28, TWIST1, PRRX1, MMPs
names(a)[a %in% c(12)] <- "Cardiovascular cells" #MYLs x10, MYHs x3, ACTA2, ACTG2, CNN1, CALD1, TNNT2
names(a)[a %in% c(13)] <- "Airway cells" #CXCL17, MUC16

Cell.Types <- as.character(object$SCT_snn_res.1)
names(Cell.Types) <- names((object$SCT_snn_res.1))
for (i in levels(object$SCT_snn_res.1)) {Cell.Types[Cell.Types == i] <- names(a[a == i])}
object$Cell.Types <- factor(Cell.Types, levels = c("Thyroid progenitors","Immature thyrocytes",
                                                   "Mature thyrocytes","Fibroblasts",
                                                   "Airway cells","Cardiovascular cells",
                                                   "Endoderm epithelial cells"))
DimPlot(object, group.by = "Cell.Types", label = T)

## OUTPUTTING FIGURES
FeaturePlot(object, features = c("NKX2-1","PAX8","FOXE1","TSHR","TG","TPO"), ncol = 3)
ggsave("figures/humanv6/ftplot_thyroid.png", units = "px", width = 6300, height = 4200)
VlnPlot(object, features = c("NKX2-1","PAX8","FOXE1","TSHR","TG","TPO"), ncol = 2, group.by = "Cell.Types")
ggsave("figures/humanv6/vlnplot_thyroid.png", units = "px", width = 4200, height = 6300)
FeaturePlot(object, features = c("KRT5","TP63"))
ggsave("figures/humanv6/ftplot_airway.png", units = "px", width = 4200, height = 2100)
FeaturePlot(object, features = c("FOXA2","ADAM28"))
ggsave("figures/humanv6/ftplot_endo.png", units = "px", width = 4200, height = 2100)
FeaturePlot(object, features = c("POU3F1","HIST1H2BG"))
ggsave("figures/humanv6/ftplot_esc.png", units = "px", width = 4200, height = 2100)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4148027/
FeaturePlot(object, features = "DCN")
ggsave("figures/humanv6/ftplot_fibro.png", units = "px", width = 2100, height = 2100)
FeaturePlot(object, features = "ACTA2")
ggsave("figures/humanv6/ftplot_cardio.png", units = "px", width = 2100, height = 2100)
DimPlot(object, group.by = "Cell.Types") + ggtitle(label = "")
ggsave("figures/humanv6/dimplot_annot.png", units = "px", width = 2100, height = 2100)

genes <- c("POU3F1","HIST1H2BG","NKX2-1","PAX8","FOXE1","HHEX","EGFP","TG","TSHR","TPO","PRRX1","DCN","COL1A2","KRT5","TP63","ACTA2","TNNT2","FOXA1","FOXA2","ADAM28")
genes <- c("NKX2-1","PAX8","FOXE1","HHEX","EGFP","TG","TSHR","TPO","PRRX1","DCN","COL1A2","KRT5","TP63","ACTA2","TNNT2","FOXA1","FOXA2","ADAM28")
mat <- matrix(ncol = length(genes), nrow = length(unique(object$Cell.Types)))
object@active.ident <- object$Cell.Types
rownames(mat) <- levels(object@active.ident)
colnames(mat) <- genes
for (g in genes){mat[,g] <- unlist(AverageExpression(object, assays = "SCT", slot = "counts",features = g))}
pdf(file = "figures/humanv6/heatmap_markers.pdf")
pheatmap(mat, cluster_rows = F, display_numbers = F, scale = "column", cluster_cols = F, cellwidth = 16, cellheight = 16)
dev.off()

active.ligands <- list(tgf = c("TGFB1","TGFB2","TGFBR1","TGFBR2"),
          bmp = c("BMP2","BMP4","BMP7","BMPR1A","BMPR2"),
          fgf = c("FGF2","FGFR1","FGFR2"),
          wnt = c("WNT2","WNT5A","FZD1","FZD2","FZD3"),
          igf = c("IGF1","IGF2","IGF1R","IGF2R"))
for(set in names(active.ligands)){
    DotPlot(object,features = active.ligands[[set]])
    ggsave(filename = paste0("figures/humanv6/dotplot_",set,".png"), units = "px",
           width = 3000, height = 2100)
}


FeaturePlot(object, features = c("POU3F1","SOX4","SOX11","KLF4","HIST1H2BG","RASD1"), ncol = 3)
ggsave("figures/humanv6/ftplot_cycling_p.png")
RidgePlot(object, features = "percent.mt")
ggsave("figures/humanv6/ridge_mt.png")

object <- SetIdent(object = object, value = object$Cell.Types)
object.markers <- FindAllMarkers(object, only.pos = T)
object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(object, features = top10$gene) + NoLegend()
ggsave("figures/humanv6/heatmap_clusters_markers.png", units = "px", width = 11000, height = 4000)
## VELOCYTO
## ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
## bm <- as.Seurat(x = ldat)
## bm <- SCTransform(object = bm, assay = "spliced")
## bm <- RunPCA(object = bm, verbose = FALSE)
## bm <- FindNeighbors(object = bm, dims = 1:20)
## bm <- FindClusters(object = bm)
## bm <- RunUMAP(object = bm, dims = 1:20)
## bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = object)))
names(x = ident.colors) <- levels(x = object)
cell.colors <- ident.colors[Idents(object = object)]
names(x = cell.colors) <- colnames(x = object)
png(filename = "veloplot_nomt.png", width = 960, height = 1080)
show.velocity.on.embedding.cor(emb = Embeddings(object = object, reduction = "umap"), vel = Tool(object = object, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)
dev.off()

## MONOCLE
seurat <- object
subsetting <- c(0,1,2,4,5,9,10)
cells.keep <- colnames(seurat)[seurat$SCT_snn_res.1 %in% subsetting]
seurat.subset <- subset(seurat, cells = cells.keep)
cds <- as.cell_data_set(DietSeurat(seurat.subset, graphs = "umap"))
## cds <- preprocess_cds(cds, num_dim = 50)
## cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot(wrap_plots(p1, p2))
## seurat.subset <- as.Seurat(cds)
cds <- learn_graph(cds, close_loop = F)
cds <- order_cells(cds)

DimPlot(seurat.subset,
        ## group.by = "SCT_snn_res.0.1",
        group.by = "CellTypes",
        ## group.by = "orig.ident",
        label = T)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE, 
           label_branch_points = FALSE)

object <- AddMetaData(
  object = object,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Thyroid.pseudotime"
)
FeaturePlot(object, features = "Thyroid.pseudotime", cols = c("green","red")) + ggtitle(label = "")
ggsave("figures/humanv6/ftplot_ptime.png", units = "px", width = 2100, height = 2100)

## ADDING METADATA
thyroid.genes <- thyroid.genes[!thyroid.genes %in% c("CDH1","FN1")]
object <- AddModuleScore(object, features = list(thyroid = thyroid.genes,
                                                 meso = c("COL16A2","MFAP4","DCN","FHL1", "COL1A2","COL3A1"),
                                                 pluri = c("SOX2", "POU5F1", "NANOG", "MYC", "KLF4", "FUT4")))
names(object@meta.data)[names(object@meta.data) == "Cluster1"] <- "ThyroidMeta"
names(object@meta.data)[names(object@meta.data) == "Cluster2"] <- "MesoMeta"
names(object@meta.data)[names(object@meta.data) == "Cluster3"] <- "PluriMeta"
FeaturePlot(object, features = "ThyroidMeta")

####################################### CELLPHONEDB ANALYSIS #######################################
###### plotting function from Elif 
## !!! means_path and pvalues arguments
dot_plot = function(selected_rows =  rows$V1,
                    selected_columns = columns$V1,
                    filename = './out/plot_human_d45.pdf',
                    width = 10,
                    height = 13,
                    means_path = './out/means.txt',
                    pvalues_path = './out/pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_line(),
          panel.grid.major = element_line(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}

## extracting data from seurat object for cellphonedb
setwd("cpdb")
ct <- as.vector(object$Cell.Types)
names(ct) <- names(object$Cell.Types)
names(ct) <- gsub("-1","",names(ct))
ct <- gsub(" ","_",ct)
ct.df <- data.frame(Cell = names(ct), cell_type = ct)
m <- object@assays$RNA@counts
colnames(m) <- gsub(".1","",colnames(m))
g <- rownames(m)
ng <- mapIds(org.Hs.eg.db, keys = g, column = "ENSEMBL", keytype="SYMBOL")
ngc <- ng[!ng %in% NA]
mc <- m[!ng %in% NA,]
rownames(mc) <- ngc

write.table(ct.df, file="./hs_meta.txt", col.names=T, row.names=F,quote=F, sep="\t")
write.table(mc, file="./hs_counts.txt", quote=F, sep="\t")

## running cellphonedb statistical methods using python API

## Importing and analyzing results
means <- read.table("./out/means.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
sig.means <- read.table("./out/significant_means.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
pvals <- read.table("./out/pvalues.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
cnet <- read.table("./out/count_network.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
decon <- read.table("./out/deconvoluted.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
cint <- read.table("./out/interaction_count.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
    
## Generating plots
columns <- read.table("./elif/columns",stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)
rows <- read.table("./elif/rows.txt",stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)

dot_plot(filename = './out/plot_human_d45_full.pdf')

## GENES OVER PSEUDOTIME TRAJECTORY
colors <- c("TG" = "red","TSHR" = "blue","TPO" = "green",
            "NKX2.1" = "purple","PAX8" = "orange","FOXE1" = "yellow")
genes <- names(colors)
counts <- object@assays$SCT@counts

df <- data.frame(Pseudotime = object$Thyroid.pseudotime, TG = counts["TG",],
                 TSHR = counts["TSHR",], TPO = counts["TPO",], NKX2.1 = counts["NKX2-1",],
                 PAX8 = counts["PAX8",], FOXE1 = counts["FOXE1",])
for (g in genes){
    df[,g] <- (df[,g] - mean(df[,g])) / sd(df[,g])
}

gg <- ggplot(df, aes(x = Pseudotime)) + ylab("Normalized expression") +
    scale_color_manual(name = NULL,values = colors) + theme_classic()
for (gene in genes){
    ## df.loess <- loess(data = df, formula = FOXE1 ~ Pseudotime, parametric = F)
    gg <- gg + geom_smooth(mapping = aes_(y = df[,gene], color = gene), se = F,
                           method = loess)
        ## geom_hline(aes(yintercept = max(predict(df.loess))))
}
gg

ggsave(gg, filename = "figures/humanv6/ptime_thyroid_genes.png", units = "px", width = 4700, height = 2100)

## SAVING OBJECT AND CELL COUNTS
setwd("~/Desktop/organoids/")
numbers.cells <- table(object$Cell.Types)
write.table(x = numbers.cells, file = "cell_counts_human.tsv",
            sep = "\t", row.names = F, col.names =F)
saveRDS(object, file = "HS_ORGANOIDSC1_seurat_mt26.Rds")
