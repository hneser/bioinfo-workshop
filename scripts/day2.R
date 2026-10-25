#kullanılacak kitaplıklar
library(Seurat)
library(tidyverse)
library(patchwork)

#Matrix i okumak
counts <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
dim(counts)

#Örnek genlere bakalım
counts[c("CD3D", "MS4A1", "LYZ"), 1:6]

format(object.size(counts), units = "MB")
sprintf("as an ordinary (dense) matrix: %.0f MB", prod(dim(counts)) * 8 / 1e6)

#Okuma dosyasının yüzde kaçı sıfırlardan oluşur?
zero_share <- 1 - Matrix::nnzero(counts) / prod(dim(counts))
zero_share

#Seurat objesi oluşturma (en az 200 farklı RNA içersin ve her bir RNA için en az 3 hücrede içeren olsun)
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)

View(pbmc)

pbmc

#dosyayı kontrol et
head(pbmc[[]])

#QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
summary(pbmc[[]][, c("nCount_RNA", "nFeature_RNA", "percent.mt")])

#Violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Total vs Unique RNA
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Esik değerler için görsel oluşturma
meta <- pbmc[[]] |> as_tibble(rownames = "barcode")
meta |>
  pivot_longer(c(nFeature_RNA, nCount_RNA, percent.mt), names_to = "metric") |>
  ggplot(aes(value)) +
  geom_histogram(bins = 80) +
  facet_wrap(~ metric, scales = "free") +
  theme_bw()

#Filtreye karar vermek için farklı kombinasyonları denemek
thresholds <- expand_grid(max_features = c(2000, 2500, 3000, Inf), max_mt = c(5, 10, 20)) |>
  mutate(cells_kept = map2_int(max_features, max_mt,
                               \(f, m) sum(meta$nFeature_RNA > 200 & meta$nFeature_RNA < f & meta$percent.mt < m)))
thresholds |> pivot_wider(names_from = max_mt, values_from = cells_kept, names_prefix = "percent_mt_below_")

#Biyolojik kayıp incelemesi (en az 10 tane ppbp içeren hücreleri platelet olarak etiketle 
#bu platelerin içinden en az 200,300,400 ve 500 unique RNA içerenleri seç)
platelet_like <- LayerData(pbmc, layer = "counts")["PPBP", ] >= 10
tibble(min_features = c(200, 300, 400, 500)) |>
  mutate(cells_kept = map_int(min_features, \(m) sum(meta$nFeature_RNA > m & meta$nFeature_RNA < 2500 & meta$percent.mt < 5)),
         platelet_like_kept = map_int(min_features, \(m) sum(platelet_like & meta$nFeature_RNA > m)))

#Filtreleme (min 200 olsun farklı RNA sayısı da 2500 altında kalsın dersek, mitokondriyal contentimde %5 kalsın)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc

#Normalizasyon(birbirleri ile karşılaştırılabilir hale getirmek için)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#En değişken ve bilgi verici RNA kümelerini bul
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
top10

#Plot: top10 i labellanmış şekilde göreceğiz
LabelPoints(plot = VariableFeaturePlot(pbmc), points = top10, repel = TRUE)

#Scaling
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)
pbmc

#PCA (Principle Component Analysis)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
print(pbmc[["pca"]], dims = 1:3, nfeatures = 5)

#Kaç tane PCA?
DimPlot(pbmc, reduction = "pca") + ElbowPlot(pbmc, ndims = 30)

#Kaçıncı componentten sonra gürültü geliyor
DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)

#Kümeleme işlemi (bir hücrenin diğer hücre ile benzerliğini/benzemezliğini hesaplamak gerekir)
#Önce ne kadar yakınlar/benzerler?
pbmc <- FindNeighbors(pbmc, dims = 1:10)

#Ne kadar uzaklar?
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Cozunurluk neye gore kaar veriyoruz
res_table <- map_dfr(c(0.2, 0.5, 0.8, 1.2), \(r) {
  cl <- FindClusters(pbmc, resolution = r, verbose = FALSE)
  tibble(resolution = r, clusters = nlevels(Idents(cl)))
})
res_table

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()

#Marker genleri bulma
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

#Kaç marker adayı bulduk?
nrow(markers)

#Hepsi gercekten marker olabilir mi? her bir cluster için marker olabileceğini düşündüğün genleri bana göster
markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 7, with_ties = FALSE) |>
  summarise(top_markers = paste(gene, collapse = ", "))

#bilinen markerlar (canonical markers)
canonical <- c("IL7R", "CCR7", "S100A4", "CD14", "LYZ", "MS4A1", "CD8A",
               "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")
DotPlot(pbmc, features = canonical) + RotatedAxis()

#İlgili genleri farklı clusterlar üzerine oturtalım (heatmap üzerine)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

#hucreleri adlandiralim, markerlara bagli olarak
marker_sets <- list(
  "Naive CD4 T"  = c("IL7R", "CCR7", "LEF1"),
  "Memory CD4 T" = c("IL7R", "S100A4", "CD3E"),
  "CD14+ Mono"   = c("CD14", "LYZ"),
  "B"            = c("MS4A1", "CD79A"),
  "CD8 T"        = c("CD8A", "CD8B"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
  "NK"           = c("GNLY", "NKG7", "KLRD1"),
  "DC"           = c("FCER1A", "CST3"),
  "Platelet"     = c("PPBP", "PF4")
)
avg <- AverageExpression(pbmc, features = unique(unlist(marker_sets)), layer = "scale.data",
                         group.by = "seurat_clusters")$RNA

#Her cluster hangi hücre tipi(marker) yönünden güçlü? (Hangi cluster, yukarida yazilan markerlar yonunden, hangi hucre tipini daha cok ifade ediyor)
scores <- sapply(marker_sets, \(g) colMeans(avg[intersect(g, rownames(avg)), , drop = FALSE]))
round(scores, 2)

#İsimlendirmeyi tamamlayalım
assigned <- colnames(scores)[apply(scores, 1, which.max)]
names(assigned) <- rownames(scores)
assigned

#UMAP üzerinden isimlendirelim
names(assigned) <- sub("^g", "", names(assigned))
pbmc$cell_type <- unname(assigned[as.character(pbmc$seurat_clusters)])
Idents(pbmc) <- "cell_type"
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()

#Hangi hucre tipi yuzde kac bulunuyor?
composition <- pbmc[[]] |>
  as_tibble(rownames = "barcode") |>
  count(cell_type, name = "cells") |>
  mutate(percent = round(100 * cells / sum(cells), 1)) |>
  arrange(desc(cells))
composition

#bar plot
ggplot(composition, aes(x = reorder(cell_type, cells), y = cells)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Cells", title = "Cell-type composition") + theme_bw()

#Markerlarimiz, logFC duzeyinde ne kadar farklilar??
markers_named <- markers |>
  as_tibble() |>
  mutate(cluster_cell_type = unname(assigned[as.character(cluster)])) |>
  relocate(cluster_cell_type, .after = cluster)
markers_named |> group_by(cluster_cell_type) |> slice_max(avg_log2FC, n = 2, with_ties = FALSE) |>
  select(cluster_cell_type, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

#SAVE
saveRDS(pbmc, "results/pbmc3k_annotated.rds")
file.size("results/pbmc3k_annotated.rds") / 1e6


