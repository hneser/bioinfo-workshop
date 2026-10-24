#library download
library(tidyverse)
library(tximport)
library(DESeq2)

library(pheatmap)
library(ggrepel)
library(apeglm)

#Loading samples
samples <- read_csv("data/airway_samples.csv")
samples

#deney dizaynı için kondisyon faktorleme
samples <- samples |>
  select(run, cell, dex) |>
  mutate(cell = factor(cell),
         dex = factor(dex, levels = c("untrt", "trt")))
samples

levels(samples$dex)

#salmon quant verilerini DeSeq2 için kullanılabilir hale getirmek
tx2gene <- read_tsv("data/tx2gene.gencode.v50.tsv.gz")
genes   <- read_tsv("data/genes.gencode.v50.tsv.gz")

tx2gene
genes

#elimde bu datalar var mı hazır mıyım?
own <- c("SRR1039508", "SRR1039509")
files <- if_else(samples$run %in% own,
                 file.path("results/salmon", samples$run, "quant.sf"),
                 file.path("data/quant", samples$run, "quant.sf.gz"))
names(files) <- samples$run
files
all(file.exists(files))

#dosya tipi değiştirme yapıyoruz. Salmona uyumlu hale getiriyoruz.
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

names(txi)
dim(txi$counts)

#değişkenin içinde hangi bilgiler var
round(head(txi$counts), 1)


#EDA: keşifsel data analizi
lib <- tibble(run = colnames(txi$counts),
              assigned_millions = colSums(txi$counts) / 1e6,
              genes_detected = colSums(txi$counts >= 1)) |>
  left_join(samples, by = "run")
lib

#EDA: Grafik
ggplot(lib, aes(x = run, y = assigned_millions, fill = dex)) +
  geom_col() +
  labs(x = NULL, y = "Reads assigned to genes (millions)",
       title = "Library size per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Countlarına bakıcaz
log_counts <- as_tibble(log2(txi$counts + 1), rownames = "gene_id") |>
  pivot_longer(-gene_id, names_to = "run", values_to = "log2_count")
ggplot(log_counts |> filter(log2_count > 0), aes(x = log2_count, colour = run)) +
  geom_density() +
  labs(x = "log2(count + 1), genes with at least one read", y = "Density",
       title = "Count distributions") +
  theme_bw()

#DeSeq2 data oluşturmalıyız, coldata dediğimiz: metadata
coldata <- samples |> column_to_rownames("run")
coldata

#Normal DeSeq2 da olmayan ancak Salmon için zorunlu olan bir adım
DESeqDataSetFromTximport(txi, colData = coldata[c(1, 2, 3:8), ], design = ~ cell + dex)
#NOT:Bu ekrandaki 78. satırda yer alan kod Salmondan gelen datalar için istisnai bir adımdır

#adds
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ cell + dex)

#Düşük ifadeli ve bilgi vermeyecek genleri temizleme
filter_table <- expand_grid(min_count = c(1, 5, 10, 20), min_samples = c(2, 4, 6)) |>
  mutate(genes_kept = map2_int(min_count, min_samples,
                               \(m, s) sum(rowSums(counts(dds) >= m) >= s)))
filter_table |> pivot_wider(names_from = min_samples, values_from = genes_kept,
                            names_prefix = "in_at_least_")

nrow(dds)

#Filtreleme: en az 4 örnekte ifade edilen min 10 örnek seçersek
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
nrow(dds)

#sample level EDA ve PCA VST (örnekler arası varyasyona bakıyoruz)
vsd <- vst(dds, blind = TRUE)

#PCA (varyansı en yüksek 500 geni kullanarak yaptık)
pca <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)

pv <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, colour = dex, shape = cell)) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = name), size = 2.8, show.legend = FALSE) +
  labs(x = paste0("PC1: ", pv[1], "% of variance"), y = paste0("PC2: ", pv[2], "% of variance"),
       title = "PCA of the 500 most variable genes (VST)") +
  theme_bw()

#örnek uzaklığı: heatmap (uzamsal uzaklıklarına) 
sample_dist <- dist(t(assay(vsd)))
ann <- as.data.frame(colData(vsd)[, c("dex", "cell")])
pheatmap(as.matrix(sample_dist),
         clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist,
         annotation_col = ann, main = "Euclidean distances between samples (VST)")

#Differential expression
dds <- DESeq(dds)

#Dispersion plot (data seti içindeki davranışı diğer genlerle uyumlu olanları seçiyoruz)
plotDispEsts(dds, main = "Dispersion estimates")

#kontrastı neyin üzerinden yapacağını söylüyoruz, burada q-value üzerinden değerlendirme yapıcaz 
#burada LFC değerleri ile ilgileniyoruz, q-value normalizasyonu yapmalıyız
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
summary(res)


#Esik değer için inceleme
threshold_table <- expand_grid(padj_below = c(0.01, 0.05, 0.1), abs_log2FC_above = c(0, 0.585, 1)) |>
  mutate(genes = map2_int(padj_below, abs_log2FC_above,
                          \(p, l) sum(res$padj < p & abs(res$log2FoldChange) > l, na.rm = TRUE)))
threshold_table |> pivot_wider(names_from = abs_log2FC_above, values_from = genes,
                               names_prefix = "abs_log2FC_above_")

#Esik değer secimi
res_lfc1 <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05, lfcThreshold = 1)
sum(res_lfc1$padj < 0.05, na.rm = TRUE)

#gürültülerini belirleme
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", res = res, type = "apeglm")

#anoyaston ve siralama
res_tbl <- as_tibble(as.data.frame(res), rownames = "gene_id") |>
  left_join(genes, by = "gene_id") |>
  left_join(as_tibble(as.data.frame(res_shrunk), rownames = "gene_id") |>
              select(gene_id, log2FC_shrunk = log2FoldChange), by = "gene_id") |>
  relocate(gene_name, gene_type, chromosome, .after = gene_id) |>
  arrange(padj)
res_tbl |> select(gene_name, baseMean, log2FoldChange, log2FC_shrunk, padj) |> head(15)

#Up ve Down önemli gen sayisi
sig <- res_tbl |> filter(padj < 0.05, abs(log2FoldChange) > 1)
sig |> dplyr::count(direction = if_else(log2FoldChange > 0, "up in treated", "down in treated"))


#volcano plot
volcano <- res_tbl |>
  filter(!is.na(padj)) |>
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1,
         padj_plot = pmax(padj, 1e-300))
ggplot(volcano, aes(x = log2FoldChange, y = -log10(padj_plot), colour = significant)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_text_repel(data = slice_min(volcano, padj_plot, n = 12, with_ties = FALSE),
                  aes(label = gene_name), size = 3, colour = "black", max.overlaps = 30) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "log2 fold change (treated / untreated)", y = "-log10 adjusted p-value",
       title = "Dexamethasone versus untreated") +
  theme_bw()

#yukarıdan aşağı genleri sıraladığımızda kaç tane genin reliable olduğuna bakıyoruz
top6 <- res_tbl |> filter(!is.na(padj)) |> slice_head(n = 6)
norm_long <- counts(dds, normalized = TRUE)[top6$gene_id, ] |>
  as_tibble(rownames = "gene_id") |>
  pivot_longer(-gene_id, names_to = "run", values_to = "normalized_count") |>
  left_join(samples, by = "run") |>
  left_join(genes |> select(gene_id, gene_name), by = "gene_id")
norm_long

#Sample heatmap for top 30 gene
top30 <- res_tbl |> filter(padj < 0.05) |> slice_min(padj, n = 30, with_ties = FALSE)
m <- assay(vsd)[top30$gene_id, ]
m <- m - rowMeans(m)
rownames(m) <- make.unique(top30$gene_name)
pheatmap(m, annotation_col = ann, main = "Top 30 genes, VST, centred per gene")

