---
title: "SBL Cambodia Script Analysis"
output: html_document
date: "2025-11-12"
---

# Data Information

**Project:** Léa Jobert Project – SBBL Cambodia experiment  

**Experiment:**  
Changes in the microbiome associated with diseased rice in Cambodia:  
What is the legacy for the growth and resistance of the next generation?  

**Plant species:** Rice var. *Srangae Sral*  
**Compartments:** Rhizosphere, root, leaf, and soil  
**Inoculation:** *Xoo*, *XOC*, and field manipulation  
Amplicon sequencing of *ITS* and *ITS*  

---

## 1 – Common Setup

### 1.1 Install required packages
```{r message=FALSE, warning=FALSE, include=FALSE}
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(microViz); packageVersion("microViz")
library(dplyr); packageVersion("dplyr")
library(data.table); packageVersion("data.table")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(DECIPHER); packageVersion("DECIPHER")
library(picante); packageVersion("picante")
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(vctrs); packageVersion("vctrs")
library(microbiomeutilities); packageVersion("microbiomeutilities")
library(pairwiseAdonis)
library(devtools)
library(pheatmap)
library(RColorBrewer)
library(decontam); packageVersion("decontam")
library(ranacapa)
```


```{r message=FALSE, warning=FALSE, include=FALSE}
# Set working directory and file paths
setwd("your_path/ITS")
path <- "your_path/ITS"
list.files(path)

# Load the main CSV file
file_path <- file.path(path, "SBL_ITS_tout_R.csv")
df <- read.csv(file_path, sep=";", header=TRUE, stringsAsFactors=FALSE)

# Create taxa table with ASVs as row names
taxa <- df[, c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
rownames(taxa) <- taxa$ASV
taxa$ASV <- NULL
saveRDS(taxa, file.path(path, "datataxamodifie.RDS"))

# Create ASV table (samples as rows, ASVs as columns)
asv <- df[, setdiff(names(df), c("taxonomy", "ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]
rownames(asv) <- asv$ID
asv$ID <- NULL
asv <- t(asv)
colnames(asv) <- df$ASV
saveRDS(asv, file.path(path, "dataASVmodifie.RDS"))

# Load RDS and metadata
taxa <- readRDS(file.path(path, "datataxamodifie.RDS"))
metadata <- read.csv2(file.path(path, "metadataSBL.csv"))
asvITS <- readRDS(file.path(path, "dataASVmodifie.RDS"))

# Ensure sample IDs match between ASV table and metadata
all(rownames(asvITS) %in% metadata$SampleID)  # Expected TRUE

# Keep only relevant metadata columns
rownames(metadata) <- metadata$SampleID
keep.cols <- c("SampleID", "Generation", "Initialsoil", "Compartment", "Condition", "Group", "Status", "Compartment_Status")
samples.out <- rownames(asvITS)
metadata <- metadata[metadata$SampleID %in% samples.out, keep.cols]
rownames(metadata) <- metadata$SampleID

# Convert ASV table to numeric
asvITS <- apply(asvITS, 2, as.numeric)
rownames(asvITS) <- samples.out

# Create phyloseq object
otu_tab <- otu_table(as.matrix(asvITS), taxa_are_rows=FALSE)
taxa_tab <- tax_table(as.matrix(taxa))
psITS_1 <- phyloseq(otu_tab, sample_data(metadata), taxa_tab)
psITS_1  

# Taxonomic filtering
# Remove non-target kingdoms and unassigned taxa
psITS_1 <- subset_taxa(psITS_1, Kingdom != "Viridiplantae")
psITS_1 <- subset_taxa(psITS_1, Kingdom != "Alveolata")
psITS_1 <- subset_taxa(psITS_1, Kingdom != "Stramenopila")
psITS_1 <- subset_taxa(psITS_1, Kingdom != "Rhizaria")

# Remove NA in Kingdom and Phylum
psITS_2 <- subset_taxa(psITS_1, !is.na(Kingdom) & !is.na(Phylum))

# Remove mitochondria and chloroplast
psITS_3 <- subset_taxa(psITS_2, Family != "Mitochondria" | is.na(Family) & Class != "Chloroplast" | is.na(Class))
psITS_4 <- subset_taxa(psITS_3, !is.na(Kingdom) & !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
psITS <- psITS_4

# Remove unwanted conditions
psITS <- subset_samples(psITS, !(Condition %in% c("Diseased2", "Xoc")))
psITS <- prune_taxa(taxa_sums(psITS) >= 1, psITS)

# Keep taxa with at least 10 reads
psITS10 <- prune_taxa(taxa_sums(psITS) >= 10, psITS)
```


#Figure 1: Alpha and beta-diversity analyses of field sample microbiota.
## Figure 1B - Alpha Diversity (ITS) - Generation 0
```{r}
library(ggpubr)
library(cowplot)

# Subset Generation 0
psITSgen0 <- subset_samples(psITS, Generation == "Gen0")
psITSgen0 <- prune_taxa(taxa_sums(psITSgen0) >= 1, psITSgen0)

# Calculate Observed richness
alpha_div <- estimate_richness(psITSgen0, measures = c("Observed")) %>%
  mutate(
    SampleID = rownames(.),
    Compartment = sample_data(psITSgen0)$Compartment,
    Status = sample_data(psITSgen0)$Status,
    Group = paste0(Compartment, "_", Status),
    Compartment = factor(Compartment, levels = c("Leaf", "Root", "Rhizosphere")),
    Status = factor(Status, levels = c("Healthy", "Diseased")),
    Group = factor(Group, levels = c("Leaf_Healthy", "Leaf_Diseased",
                                     "Root_Healthy", "Root_Diseased",
                                     "Rhizosphere_Healthy", "Rhizosphere_Diseased"))
  )

# Colors
Couleurs_gen0tout <- c("Leaf_Healthy" = "#AFD1E9", "Leaf_Diseased" = "#369686",
                       "Root_Healthy" = "#EDDF82", "Root_Diseased" = "#E3CD3B",
                       "Rhizosphere_Healthy" = "#927763", "Rhizosphere_Diseased" = "#55463A")

Couleurs_hd <- c("Healthy" = "#67B85A", "Diseased" = "#E97428")

# Plot
p <- ggplot(alpha_div, aes(x = Group, y = Observed, fill = Group, shape = Status)) +
  geom_boxplot(aes(color = Status), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Status), position = position_jitter(0.2), size = 2.5, alpha = 0.7) +
  scale_fill_manual(values = Couleurs_gen0tout) +
  scale_color_manual(values = Couleurs_hd) +
  theme_minimal(base_size = 17) +
  ggtitle("Alpha Diversity - Observed Richness (ITS)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(labels = c("Healthy", "Diseased", "Healthy", "Diseased", "Healthy", "Diseased")) +
  xlab(NULL) 

# Wilcoxon tests
comparisons <- list(
  c("Leaf_Healthy", "Leaf_Diseased"),
  c("Root_Healthy", "Root_Diseased"),
  c("Rhizosphere_Healthy", "Rhizosphere_Diseased")
)

p_stat <- p + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",
  label = "p.signif",
  hide.ns = FALSE,
  tip.length = 0.01,
  size = 6
)

# Save
ggsave(
  filename = "your_path/ITS/Figures/AlphaDiversity_plots_gen0/AlphaDiversity_Observed_Grouped_with_stats.svg",
  plot = p_stat,
  width = 9,
  height = 6,
  units = "in",
  device = "svg"
)
```




## Figure 1F-G-H - Beta Diversity (ITS) - Generation 0
```{r}
## ITS - PCoA per compartment
## Figures 1G (Root), 1H (Rhizosphere), 1F (Leaf)

# Output directory
output_dir_ITS_pcoa <- "your_path/PCoA_par_compartiment"
dir.create(output_dir, showWarnings = FALSE)

# Colors
line_colors <- c("Healthy" = "#67B85A", "Diseased" = "#E97428")
fill_root <- c("Root_Healthy" = "#EDDF82", "Root_Diseased" = "#E3CD3B")
fill_rhizo <- c("Rhizosphere_Healthy" = "#927763", "Rhizosphere_Diseased" = "#55463A")
fill_leaf <- c("Leaf_Healthy" = "#AFD1E9", "Leaf_Diseased" = "#369686")

## Figure 1G - Root (ITS)
ps_root <- subset_samples(psITS, Compartment == "Root")
ps_root <- prune_taxa(taxa_sums(ps_root) >= 1, ps_root)
ps_root_tss <- transform_sample_counts(ps_root, function(x) x / sum(x))

bray_root <- phyloseq::distance(ps_root_tss, method = "bray")
ord_root <- ordinate(ps_root_tss, method = "PCoA", distance = "bray")
ord_root_data <- plot_ordination(ps_root_tss, ord_root, type = "samples", color = "Status", shape = "Status")$data

metadata_root <- as(sample_data(ps_root_tss), "data.frame")
permanova_root <- adonis2(bray_root ~ Status, data = metadata_root, permutations = 999)
pval_root <- formatC(permanova_root$`Pr(>F)`[1], format = "e", digits = 2)
r2_root <- formatC(permanova_root$R2[1], format = "f", digits = 2)

eig_root <- ord_root$values$Relative_eig * 100
x_lab_root <- paste0("PCoA1 (", round(eig_root[1], 1), "%)")
y_lab_root <- paste0("PCoA2 (", round(eig_root[2], 1), "%)")

ord_root_data$Status <- factor(ord_root_data$Status, levels = c("Healthy", "Diseased"))
ord_root_data$Compartment_Status <- paste0("Root_", ord_root_data$Status)

p_root <- ggplot(ord_root_data, aes(x = Axis.1, y = Axis.2, color = Status, shape = Status)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Compartment_Status), geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(color = Status), geom = "path", linewidth = 1) +
  scale_color_manual(values = line_colors) +
  scale_fill_manual(values = fill_root) +
  scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 17)) +
  labs(x = x_lab_root, y = y_lab_root,
       title = paste0("PCoA - Root\nPERMANOVA p = ", pval_root, "**, R² = ", r2_root)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.position = "right")

ggsave(file.path(output_dir, "PCoA_ITS_Root_Fig1F.svg"), plot = p_root, width = 7, height = 5, device = "svg")


## Figure 1H - Rhizosphere (ITS)
ps_rhizo <- subset_samples(psITS, Compartment == "Rhizosphere")
ps_rhizo <- prune_taxa(taxa_sums(ps_rhizo) >= 1, ps_rhizo)
ps_rhizo_tss <- transform_sample_counts(ps_rhizo, function(x) x / sum(x))

bray_rhizo <- phyloseq::distance(ps_rhizo_tss, method = "bray")
ord_rhizo <- ordinate(ps_rhizo_tss, method = "PCoA", distance = "bray")
ord_rhizo_data <- plot_ordination(ps_rhizo_tss, ord_rhizo, type = "samples", color = "Status", shape = "Status")$data

metadata_rhizo <- as(sample_data(ps_rhizo_tss), "data.frame")
permanova_rhizo <- adonis2(bray_rhizo ~ Status, data = metadata_rhizo, permutations = 999)
pval_rhizo <- formatC(permanova_rhizo$`Pr(>F)`[1], format = "e", digits = 2)
r2_rhizo <- formatC(permanova_rhizo$R2[1], format = "f", digits = 2)

eig_rhizo <- ord_rhizo$values$Relative_eig * 100
x_lab_rhizo <- paste0("PCoA1 (", round(eig_rhizo[1], 1), "%)")
y_lab_rhizo <- paste0("PCoA2 (", round(eig_rhizo[2], 1), "%)")

ord_rhizo_data$Status <- factor(ord_rhizo_data$Status, levels = c("Healthy", "Diseased"))
ord_rhizo_data$Compartment_Status <- paste0("Rhizosphere_", ord_rhizo_data$Status)

p_rhizo <- ggplot(ord_rhizo_data, aes(x = Axis.1, y = Axis.2, color = Status, shape = Status)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Compartment_Status), geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(color = Status), geom = "path", linewidth = 1) +
  scale_color_manual(values = line_colors) +
  scale_fill_manual(values = fill_rhizo) +
  scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 17)) +
  labs(x = x_lab_rhizo, y = y_lab_rhizo,
       title = paste0("PCoA - Rhizosphere\nPERMANOVA p = ", pval_rhizo, "**, R² = ", r2_rhizo)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.position = "right")

ggsave(file.path(output_dir, "PCoA_ITS_Rhizo_Fig1D.svg"), plot = p_rhizo, width = 7, height = 5, device = "svg")


## Figure 1F - Leaf (ITS)
ps_leaf <- subset_samples(psITS, Compartment == "Leaf")
ps_leaf <- prune_taxa(taxa_sums(ps_leaf) >= 1, ps_leaf)
ps_leaf_tss <- transform_sample_counts(ps_leaf, function(x) x / sum(x))

bray_leaf <- phyloseq::distance(ps_leaf_tss, method = "bray")
ord_leaf <- ordinate(ps_leaf_tss, method = "PCoA", distance = "bray")
ord_leaf_data <- plot_ordination(ps_leaf_tss, ord_leaf, type = "samples", color = "Status", shape = "Status")$data

metadata_leaf <- as(sample_data(ps_leaf_tss), "data.frame")
permanova_leaf <- adonis2(bray_leaf ~ Status, data = metadata_leaf, permutations = 999)
pval_leaf <- formatC(permanova_leaf$`Pr(>F)`[1], format = "e", digits = 2)
r2_leaf <- formatC(permanova_leaf$R2[1], format = "f", digits = 2)

eig_leaf <- ord_leaf$values$Relative_eig * 100
x_lab_leaf <- paste0("PCoA1 (", round(eig_leaf[1], 1), "%)")
y_lab_leaf <- paste0("PCoA2 (", round(eig_leaf[2], 1), "%)")

ord_leaf_data$Status <- factor(ord_leaf_data$Status, levels = c("Healthy", "Diseased"))
ord_leaf_data$Compartment_Status <- paste0("Leaf_", ord_leaf_data$Status)

p_leaf <- ggplot(ord_leaf_data, aes(x = Axis.1, y = Axis.2, color = Status, shape = Status)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Compartment_Status), geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(color = Status), geom = "path", linewidth = 1) +
  scale_color_manual(values = line_colors) +
  scale_fill_manual(values = fill_leaf) +
  scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 17)) +
  labs(x = x_lab_leaf, y = y_lab_leaf,
       title = paste0("PCoA - Leaf\nPERMANOVA p = ", pval_leaf, "**, R² = ", r2_leaf)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.position = "right")

ggsave(file.path(output_dir, "PCoA_ITS_Leaf_Fig1E.svg"), plot = p_leaf, width = 7, height = 5, device = "svg")


```



#Figure2_Top 30 microbial genera in the roots and rhizosphere of symptomatic and asymptomatic field rice.
```{r}
## Figure 2 - ITS Top 30 genera 
#script for Fig. 2 B) Rhizosphere ITS, same to be done with root (Fig 2. D)

# Load required package
library(ggplot2)
library(phyloseq)
library(microbiome)

# Output directory
output_dir <- "your_path/ITS/Fig2_TaxonomicBinning"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## 1. Subset to Generation 0 - Rhizosphere samples
psITS_gen0 <- subset_samples(psITS, Generation == "Gen0")
psITS_gen0_rhizo <- subset_samples(psITS_gen0, Compartment == "Rhizosphere")
psITS_gen0_rhizo <- prune_taxa(taxa_sums(psITS_gen0_rhizo) >= 1, psITS_gen0_rhizo)

## 2. Agglomerate taxa to the Genus level
psITS_genus <- tax_glom(psITS_gen0_rhizo, taxrank = "Genus", NArm = TRUE)

## 3. Transform to relative abundances
psITS_genus_rel <- microbiome::transform(psITS_genus, "compositional")

## 4. Select the top 30 genera across all samples
genus_sums <- taxa_sums(psITS_genus_rel)
top30_ids <- names(sort(genus_sums, decreasing = TRUE))[1:30]
psITS_top30 <- prune_taxa(top30_ids, psITS_genus_rel)

## 5. Define custom colors for the 30 genera
custom_cols <- c(
  "#C70E7B", "#A9CBB7", "#DD1C1A", "#54BCD1", "#98FB98", "#F4B95A",
  "#009F3F", "#8FDA04", "#AF6125", "#FF499E", "#B25D91", "#EFC7E6",
  "#EF7C12", "#A52A2A", "#D8BFD8", "#F0E68C", "#8A2BE2", "#DCCCA3",
  "#FF4500", "#7B68EE", "#2E8B57", "#D295BF", "#9ACD32", "#D2691E",
  "#20B2AA", "#FFD700", "#00CED1", "#DC143C", "#556B2F", "#FF6347"
)

## 6. Ensure Condition factor order and rename samples
sample_data(psITS_top30)$Condition <- factor(
  sample_data(psITS_top30)$Condition,
  levels = c("Healthy", "Diseased1") # Healthy first
)

sample_data_df <- as.data.frame(sample_data(psITS_top30))
new_names <- ifelse(
  sample_data_df$Condition == "Healthy",
  paste0("RhH", ave(rep(1, nrow(sample_data_df)), sample_data_df$Condition, FUN = seq_along)),
  paste0("RhD", ave(rep(1, nrow(sample_data_df)), sample_data_df$Condition, FUN = seq_along))
)
sample_names(psITS_top30) <- new_names

## 7. Plot stacked bar chart of top 30 genera
p <- plot_bar(psITS_top30, fill = "Genus") +
  scale_fill_manual(values = custom_cols) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right"
  ) +
  facet_wrap(~Condition, scales = "free_x") +
  labs(x = "Rhizosphere of field plants", y = "Relative abundance") +
  ggtitle("Top 30 genera - ITS, Generation 0, Rhizosphere")

## 8. Save the figure
ggsave(
  filename = file.path(output_dir, "Taxonomic_Binning_Top30Genus_ITS_Gen0_Rhizo.svg"),
  plot = p,
  dpi = 300,
  device = "svg",
  width = 12,
  height = 6
)

```



#Figure 4: Comparison of microbiota between healthy and disease rhizosphere mix. 
##Figure4. A)
```{r}
#Loading packages
library(phyloseq)
library(vegan)
library(ggplot2)

# output file
output_dir_final <- "your_path"
dir.create(output_dir_final, showWarnings = FALSE)

####### PART 1 — G0 Rhizo vs Soil #######

psITS_gen0_sol <- subset_samples(psITS10, Compartment %in% c("Rhizosphere", "Soil"))
psITS_gen0_sol <- prune_taxa(taxa_sums(psITS_gen0_sol) >= 1, psITS_gen0_sol)
psITS_gen0_solTSS <- transform_sample_counts(psITS_gen0_sol, function(x) x / sum(x))

colors_gen0_sol <- c("G0-Rhizo-D1" = "#E97428", "G0-Rhizo-H" = "#67B85A", 
                     "MixsolH" = "#3A7231", "MixsolXoo" = "#CB5E15")

# Bray-Curtis distance
bray_dist <- phyloseq::distance(psITS_gen0_solTSS, method = "bray")

# PCoA ordination
ordination_pcoa <- ordinate(psITS_gen0_solTSS, method = "PCoA", distance = "bray")
ordination_data_pcoa <- plot_ordination(psITS_gen0_solTSS, ordination_pcoa, 
                                        type = "samples", color = "Group", shape = "Compartment")$data

# Métadonnées et tests stats
metadata <- as(sample_data(psITS_gen0_solTSS), "data.frame")
permanova <- adonis2(bray_dist ~ Group, data = metadata, permutations = 999)
pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
r2val <- formatC(permanova$R2[1], format = "f", digits = 2)

dispersion <- betadisper(bray_dist, metadata$Group)
disp_test <- permutest(dispersion)
p_disp <- formatC(disp_test$tab[1, "Pr(>F)"], format = "e", digits = 2)


#  % de variance expliquée par axe
eig_vals <- ordination_pcoa$values$Relative_eig * 100
x_lab <- paste0("PCoA1 (", round(eig_vals[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(eig_vals[2], 1), "%)")



# Graph PCoA
p.PCoA_gen0_sol <- ggplot(ordination_data_pcoa, aes(x = Axis.1, y = Axis.2, color = Group, shape = Compartment)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, color = NA) +
  scale_color_manual(values = colors_gen0_sol) +
  scale_fill_manual(values = colors_gen0_sol) +
  scale_shape_manual(values = c("Soil" = 8, "Rhizosphere" = 16)) +
  labs(
    x = x_lab,
    y = y_lab,
    title = paste0("PCoA - Transmission field → rhizosphere mix\n",
                   "PERMANOVA p = ", pval, ", R² = ", r2val,
                   " | Dispersion p = ", p_disp)
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.position = "right",
        legend.box = "vertical")

print(p.PCoA_gen0_sol)
ggsave(file.path(output_dir_final, "Fig4_PCoA_ordination_champs_sol_ITS.svg"), 
       plot = p.PCoA_gen0_sol, width = 7, height = 5, device = "svg")
```




##Figure4. C)
```{r}
#test deseq2 at genus level
# loading packages
library(phyloseq)
library(DESeq2)
library(dplyr)


# subset "soil"
psITS_soil <- subset_samples(psITS10, Compartment %in% c("Soil"))
psITS_soil <- prune_taxa(taxa_sums(psITS_soil) > 0, psITS_soil)


psITS_soil_filtered <- phyloseq::filter_taxa(psITS_soil, function(x) {
  prevalence <- sum(x > 0) / length(x)
  rel_abundance <- sum(x) / sum(otu_table(psITS_soil))
  prevalence >= 0.0 & rel_abundance >= 0.001
}, prune = TRUE)
psITS_soil_filtered


# Agglom to genus level
psITS_soil_genus <- tax_glom(psITS_soil_filtered, taxrank = "Genus", NArm = TRUE)
dds <- phyloseq_to_deseq2(psITS_soil_genus, ~ Condition)
dds <- estimateSizeFactors(dds, type = "poscounts")  # méthode tolérante aux zéros
dds <- DESeq(dds)


# Comparison
res <- results(dds, contrast = c("Condition", "Diseased1", "Healthy"))  # log2FC > 0 : enrichi dans G1Xoo
res <- as.data.frame(res[order(res$padj), ])
res_sig <- res[res$padj < 0.05, ]

# Add genus names
res_sig$taxon <- rownames(res_sig)


tax_table_df <- as.data.frame(tax_table(psITS_soil_genus))
tax_table_df$taxon <- rownames(tax_table_df)
res_sig_annotated <- merge(res_sig, tax_table_df, by = "taxon")

# Barplot log2FC significatifs
library(ggplot2)

deseq_barplot <- ggplot(res_sig_annotated, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "#1B9E77"),
                    labels = c("Enriched in G1H", "Enriched in G1Xoo")) +
  labs(title = "Differentially Transmitted Genera (DESeq2)",
       y = "log2 Fold Change (G1Xoo vs G1H)",
       x = "Genus") +
  theme_minimal(base_size = 14)
deseq_barplot
#ggsave("barplot_deseq_differentially_transmitted_gen0_gen1_soil.svg", deseq_barplot, width = 10, height = 10)



# 📦 Charger les packages nécessaires
library(openxlsx)
library(dplyr)

# 📥 Extraction des tables de base
otu <- as.data.frame(otu_table(psITS_soil_genus))
if (!taxa_are_rows(psITS_soil_genus)) {
  otu <- t(otu)
}
otu <- as.data.frame(otu)
otu$taxon <- rownames(otu)

tax <- as.data.frame(tax_table(psITS_soil_genus))
tax$taxon <- rownames(tax)

res$taxon <- rownames(res)

# 🔍 Récupérer seulement les genres significatifs
taxons_enrichis <- rownames(res_sig)
taxons_all <- rownames(res)

# 🧬 Sélectionner les abundances
otu_enrichis <- otu[otu$taxon %in% taxons_enrichis, , drop = FALSE]
otu_all <- otu[otu$taxon %in% taxons_all, , drop = FALSE]

# 🧾 Fusion complète : abondances + taxo + stats
res_sig_full <- otu_enrichis %>%
  left_join(tax, by = "taxon") %>%
  left_join(res, by = "taxon")

res_all_full <- otu_all %>%
  left_join(tax, by = "taxon") %>%
  left_join(res, by = "taxon")

# ✅ Réorganiser les colonnes dans un ordre clair
colonnes_base <- c("taxon", "log2FoldChange", "pvalue", "padj",
                   "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

autres_colonnes <- setdiff(colnames(res_sig_full), colonnes_base)
ordre_colonnes <- c(colonnes_base, autres_colonnes)
res_sig_full <- res_sig_full[, ordre_colonnes]

autres_colonnes <- setdiff(colnames(res_all_full), colonnes_base)
ordre_colonnes <- c(colonnes_base, autres_colonnes)
res_all_full <- res_all_full[, ordre_colonnes]


# 📁 Créer le dossier de sortie si nécessaire
output_dir_sol <- "your_path/Deseq2_tables_genre"
dir.create(output_dir_sol, showWarnings = FALSE)

# 💾 Exporter en Excel
write.xlsx(res_sig_full, 
           file = file.path(output_dir_sol, "psITS_soil_genus_Deseq_Diseased1_vs_Healthy_significatifs.xlsx"), 
           rowNames = FALSE)

# 💾 Exporter en Excel
write.xlsx(res_all_full, 
           file = file.path(output_dir_sol, "psITS_soil_genus_Deseq_Diseased1_vs_Healthy_tout.xlsx"), 
           rowNames = FALSE)




library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)

# Lire les données DESeq depuis un fichier Excel
df <- read_excel("ypur_path/psITS_soil_genus_Deseq_Diseased1_vs_Healthy_significatifs.xlsx")
#df <- read_excel("C:/Users/2022lj002/Documents/Cambodge/SBL sequençage/Que_Xoo/ITS/Figures/Deseq2_tables_genre/psITS_soil_genus_Deseq_Diseased1_vs_Healthy_tout.xlsx")
# Afficher les premières lignes pour vérifier
head(df)


df_sig <- df %>%
  filter(padj < 0.05, !is.na(Genus))


# Moyenne des log2FoldChange par Genus
bar_data <- df_sig %>%
  group_by(Genus) %>%
  summarise(mean_lfc = mean(log2FoldChange, na.rm = TRUE),
            n = n()) %>%
  arrange(desc(mean_lfc))


# Préparation des données

# Graphique façon LEfSe
# Préparation des données
bar_data <- df_sig %>%
  group_by(Genus) %>%
  summarise(
    mean_lfc = mean(log2FoldChange, na.rm = TRUE),
    n = n()
  ) %>%
  mutate(direction = ifelse(mean_lfc > 0,
                            "Enriched in rhizosphere mix from diseased plants",
                            "Enriched in rhizosphere mix from healthy plants")) %>%
  arrange(mean_lfc)


# Barplot façon LEfSe inversée
ggplot(bar_data, aes(x = reorder(Genus, mean_lfc), y = mean_lfc, fill = direction)) +
  geom_col() +
  coord_flip() +
  
  # Ajouter les labels du côté opposé au barre
  geom_text(
    aes(
      y = 0,  # Texte centré sur l'axe zéro
      label = Genus,
      hjust = ifelse(mean_lfc > 0, 1.05, -0.05)  # Vers la gauche pour +, vers la droite pour -
    ),
    size = 6,
    fontface = "italic"
  ) +

  # Palette personnalisée
  scale_fill_manual(values = c(
    "Enriched in rhizosphere mix from diseased plants" = "#CB5E15",
    "Enriched in rhizosphere mix from healthy plants" = "#3A7231"
  )) +
  
  # Ajouter de l’espace pour laisser respirer les textes
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +

  labs(
    title = "Log2 Fold Change par genre (padj < 0.05)",
    x = NULL,
    y = "log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 13) +
theme(
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.box = "vertical",
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)


ggsave(
  filename = "your_path/Fig4_log2FoldChange_barplot_sol_genre_ITS.svg",
  plot = last_plot(),  # ou remplace par le nom de ton objet ggplot si tu l’as assigné
  width = 8,
  height = 10,
  units = "in",
  dpi = 300,
  device = "svg"
)
```





#Figure 6: Comparison of the microbiome between the rhizosphere mix and the rhizosphere of the next plant generation. A)
```{r}
####### PARTIE 2 — G1 Soil/Rhizo #######

psITS_sol_gen1 <- subset_samples(psITS10, Compartment %in% c("Soil", "Rhizo"))
psITS_sol_gen1 <- prune_taxa(taxa_sums(psITS_sol_gen1) >= 1, psITS_sol_gen1)
psITS_sol_gen1TSS <- transform_sample_counts(psITS_sol_gen1, function(x) x / sum(x))

colors_sol_gen1 <- c("G1-SoilXoo-Rhizo-Xoo" = "#D52F20", "G1-SoilH-Rhizo-Xoo" = "#3EBBA4",
                     "MixsolH" = "#3A7231", "MixsolXoo" = "#CB5E15")

# Bray distance
bray_dist <- phyloseq::distance(psITS_sol_gen1TSS, method = "bray")

# PCoA
ordination_pcoa <- ordinate(psITS_sol_gen1TSS, method = "PCoA", distance = "bray")
ordination_data_pcoa <- plot_ordination(psITS_sol_gen1TSS, ordination_pcoa, 
                                        type = "samples", color = "Group", shape = "Compartment")$data

metadata <- as(sample_data(psITS_sol_gen1TSS), "data.frame")
permanova <- adonis2(bray_dist ~ Group, data = metadata, permutations = 999)
pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
r2val <- formatC(permanova$R2[1], format = "f", digits = 2)

dispersion <- betadisper(bray_dist, metadata$Group)
disp_test <- permutest(dispersion)
p_disp <- formatC(disp_test$tab[1, "Pr(>F)"], format = "e", digits = 2)


eig_vals <- ordination_pcoa$values$Relative_eig * 100
x_lab <- paste0("PCoA1 (", round(eig_vals[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(eig_vals[2], 1), "%)")


# Graphique PCoA
p.PCoA_sol_gen1 <- ggplot(ordination_data_pcoa, aes(x = Axis.1, y = Axis.2, color = Group, shape = Compartment)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, color = NA) +
  scale_color_manual(values = colors_sol_gen1) +
  scale_fill_manual(values = colors_sol_gen1) +
  scale_shape_manual(values = c("Soil" = 8, "Rhizo" = 16)) +
  labs(
    x = x_lab,
    y = y_lab,
    title = paste0("PCoA - Rhizosphere mix → controlled conditions\n",
                   "PERMANOVA p = ", pval, ", R² = ", r2val)
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.position = "right",
        legend.box = "vertical")

print(p.PCoA_sol_gen1)
ggsave(file.path(output_dir_final, "Fig5_PCoA_ordination_sol_g1_ITS.svg"),
       plot = p.PCoA_sol_gen1, width = 8, height = 5, device = "svg")

```





#Figure 7: Transmission of rhizosphere microbial taxa and differential abundance analysis. 
##Figure 7 A) and B) ITS
```{r}
## Figure 7 – Venn Diagrams (Genus level) for ITS data

# 📦 Load required packages
library(phyloseq)
library(tidyverse)
library(ggVennDiagram)

# 🔍 Subset data to Rhizosphere & Soil samples and agglomerate at Genus level
psITSsoilrhizo_genre <- subset_samples(psITS10, Compartment %in% c("Rhizo", "Soil", "Rhizosphere"))
psITSsoilrhizo_genre <- prune_taxa(taxa_sums(psITSsoilrhizo_genre) >= 1, psITSsoilrhizo_genre)
psITSsoilrhizo_genre <- tax_glom(psITSsoilrhizo_genre, taxrank = "Genus", NArm = TRUE)

# 🧮 Function: Retrieve taxa present in a group based on prevalence & abundance
get_present_taxa_genre <- function(ps_obj, group_name, min_prevalence = 0.1, min_abundance = 0.0001) {
  
  # Keep only samples from the specified group
  samples_to_keep <- rownames(sample_data(ps_obj))[sample_data(ps_obj)$Group == group_name]
  ps_sub <- prune_samples(samples_to_keep, ps_obj)
  
  # Return empty vector if no data left
  if (nsamples(ps_sub) == 0 || ntaxa(ps_sub) == 0) return(character(0))
  
  # Extract OTU table
  otu <- as.data.frame(otu_table(ps_sub))
  if (!taxa_are_rows(ps_sub)) otu <- t(otu)
  
  # Calculate prevalence (number of samples where the taxon is present)
  prevalence <- rowSums(otu > 0)
  min_prevalence_count <- ceiling(min_prevalence * nsamples(ps_sub))
  
  # Calculate relative abundance
  total_reads <- sum(otu)
  rel_abund <- rowSums(otu) / total_reads
  
  # Keep taxa meeting both prevalence & abundance thresholds
  keep_asvs <- rownames(otu)[prevalence >= min_prevalence_count & rel_abund >= min_abundance]
  return(keep_asvs)
}

# 🔹 Extract taxa for D1 groups
g0_d1_genre   <- get_present_taxa_genre(psITSsoilrhizo_genre, "G0-Rhizo-D1", min_prevalence = 0.1, min_abundance = 0.0001)
soil_d1_genre <- get_present_taxa_genre(psITSsoilrhizo_genre, "MixsolXoo",    min_prevalence = 0.5, min_abundance = 0.001)
g1_xoo_genre  <- get_present_taxa_genre(psITSsoilrhizo_genre, "G1-SoilXoo-Rhizo-Xoo", min_prevalence = 0.5, min_abundance = 0.001)

# 📈 Venn diagram – D1 groups
venn_plot_d1 <- ggVennDiagram(
  list(
    "G0 Rhizo D1"  = g0_d1_genre,
    "Mixsol Xoo"   = soil_d1_genre,
    "G1 Rhizo Xoo" = g1_xoo_genre
  ),
  label_alpha = 0
) +
  scale_fill_gradient(low = "#FFF5F0", high = "#A81008") +
  ggtitle("Venn (Genus) – G0 Rhizo D1 / MixsolXoo / G1 Rhizo Xoo") +
  theme(text = element_text(size = 14), legend.position = "none")

# Display plot
venn_plot_d1

# Save plot
ggsave(
  "your_path/Fig7_Venn_Genus_D1_ITS.svg",
  venn_plot_d1,
  width = 5, height = 4.5
)

# 🔁 Intersection of the 3 D1 groups
genre_intersection_d1 <- Reduce(intersect, list(g0_d1_genre, soil_d1_genre, g1_xoo_genre))


# 📋 Function: Extract taxonomy + abundances for selected taxa
extract_taxo_table_genre <- function(physeq_obj, asv_list) {
  
  # Keep only taxa that exist in the phyloseq object
  asv_list <- intersect(taxa_names(physeq_obj), asv_list)
  if (length(asv_list) == 0) {
    warning("❌ No ASVs found in the phyloseq object.")
    return(data.frame())
  }
  
  # Subset phyloseq object to selected taxa
  physeq_sub <- prune_taxa(asv_list, physeq_obj)
  
  # Extract OTU table
  otu_df <- as.data.frame(otu_table(physeq_sub))
  if (!taxa_are_rows(physeq_sub)) otu_df <- t(otu_df)
  otu_df$ASV <- rownames(otu_df)
  
  # Extract taxonomy table
  tax_df <- as.data.frame(tax_table(physeq_sub))
  tax_df$ASV <- rownames(tax_df)
  
  # Merge taxonomy and abundance data
  result <- merge(tax_df, otu_df, by = "ASV")
  return(result)
}

# Extract shared genera for D1
taxo_genres_shared_d1 <- extract_taxo_table_genre(psITSsoilrhizo_genre, genre_intersection_d1)
# Optional export:
# write_csv(taxo_genres_shared_d1, "Shared_Genus_G0_Soil_G1_D1.csv")

# 🔹 Repeat process for H groups
g0_h_genre   <- get_present_taxa_genre(psITSsoilrhizo_genre, "G0-Rhizo-H", min_prevalence = 0.1, min_abundance = 0.0001)
soil_h_genre <- get_present_taxa_genre(psITSsoilrhizo_genre, "MixsolH",    min_prevalence = 0.5, min_abundance = 0.001)
g1_h_genre   <- get_present_taxa_genre(psITSsoilrhizo_genre, "G1-SoilH-Rhizo-Xoo", min_prevalence = 0.5, min_abundance = 0.001)

# 📈 Venn diagram – H groups
venn_plot_h <- ggVennDiagram(
  list(
    "G0 Rhizo H"  = g0_h_genre,
    "Mixsol H"    = soil_h_genre,
    "G1 Rhizo H"  = g1_h_genre
  ),
  label_alpha = 0
) +
  scale_fill_gradient(low = "#EDF7F5", high = "#407038") +
  ggtitle("Venn (Genus) – G0 Rhizo H / MixsolH / G1 Rhizo H") +
  theme(text = element_text(size = 14), legend.position = "none")

# Display plot
venn_plot_h

# Save plot
ggsave(
  "your_path/Fig7_Venn_Genus_H_ITS.svg",
  venn_plot_h,
  width = 5, height = 4.5
)

# 🔁 Intersection of the 3 H groups
genre_intersection_h <- Reduce(intersect, list(g0_h_genre, soil_h_genre, g1_h_genre))

# Extract shared genera for H
taxo_genres_shared_h <- extract_taxo_table_genre(psITSsoilrhizo_genre, genre_intersection_h)

# Save to CSV
write_csv(taxo_genres_shared_h, "Shared_Genus_G0_Soil_G1_H.csv")

```




##Figure 7 E) ITS
```{r}
#Deseq on transmitted genera
library(phyloseq)
library(DESeq2)

# Filter phyloseq on transmitted genera (union of the two groups, for example)
genres_transmis <- union(genre_intersection_h, genre_intersection_d1)
ps_filtered <- prune_taxa(genres_transmis, psITSsoilrhizo_genre)

# Check that the sample names and grouping variables are ok
sample_data(ps_filtered)$Group <- factor(sample_data(ps_filtered)$Group)  # Ex: "G1Xoo" vs "G1H"

# DESeq2 depuis phyloseq
dds <- phyloseq_to_deseq2(ps_filtered, ~ Group)
dds <- estimateSizeFactors(dds, type = "poscounts")  # méthode tolérante aux zéros
dds <- DESeq(dds)

# Comparaison
res <- results(dds, contrast = c("Group", "G1-SoilXoo-Rhizo-Xoo", "G1-SoilH-Rhizo-Xoo"))  
res <- as.data.frame(res[order(res$padj), ])
res_sig <- res[res$padj < 0.05, ]

# Add genera names
res_sig$taxon <- rownames(res_sig)

# Optional: merge with taxonomy table
tax_table_df <- as.data.frame(tax_table(ps_filtered))
tax_table_df$taxon <- rownames(tax_table_df)
res_sig_annotated <- merge(res_sig, tax_table_df, by = "taxon")

# Barplot des log2FC significatifs
library(ggplot2)

# Create a 'direction' column according to the log2FC sign
res_sig_annotated <- res_sig_annotated %>%
  mutate(direction = ifelse(log2FoldChange > 0,
                            "Enriched in G1Xoo",
                            "Enriched in G1H")) %>%
  arrange(log2FoldChange)

# Create a barplot with text positioned in the LEfSe style
deseq_barplot <- ggplot(res_sig_annotated, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = direction)) +
  geom_col() +
  coord_flip() +
  
 # Text positioned to the left or right of the bar, aligned with y = 0
  geom_text(
    aes(
      y = 0,
      label = Genus,
      hjust = ifelse(log2FoldChange > 0, 1.05, -0.05)
    ),
    size = 4,
    fontface = "italic"
  ) +
  
  # My colors
  scale_fill_manual(values = c(
    "Enriched in G1Xoo" = "#D52F20",
    "Enriched in G1H" = "#3EBBA4"
  )) +
  
# Add a little space on the y axis for labels
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +

  labs(
    title = "Differential Abundance of Transmitted Genera (DESeq2)",
    y = "log2 Fold Change",
    x = NULL,
    fill = "Direction"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.box = "vertical"
  )

# show plot
deseq_barplot

# Save
ggsave("your_path/Fig5_barplot_deseq_differentially_transmitted_gen0_gen1_soil_ITS.svg", deseq_barplot, width = 6, height = 15)

```







#SUPPLEMENTARY FIGURES
##S1-C)-E) rarefaction curve
```{r rarefaction curve gen0 par compartment}
library(scales)

#subset 
psITSgen0 <- subset_samples(psITS, Generation == "Gen0")
psITSgen0
psITSgen0 <- prune_taxa(taxa_sums(psITSgen0) >= 1, psITSgen0)
psITSgen0 #35098 taxa and 54 samples



# Courbe de rarefaction avec limite de l'axe des y à 10^4
RC <- ggrare(psITSgen0, step = 100, se = FALSE, color = "Compartment") + 
  labs(title = "Rarefaction curve ITS, 
       Plants from the field", 
       x = "Depth",
       y = "Number of ASVs") +
  theme_minimal(base_size = 15) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  ) +
      scale_color_manual(values = c("Leaf" = "#369686", "Root" = "#E3CD3B", "Rhizosphere" = "#55463A")) +  # Couleurs personnalisées  scale_x_continuous(limITS = c(0, 75000))

scale_x_continuous(limits = c(0, 25000))

# Afficher le graphique
print(RC)
# Enregistrer le graphique
ggsave(RC, dpi=300, device = svg, width = 7, height = 4, filename = "your_path/Supp_1_Rarefaction curve tout par compartiment gen0 ITS.svg")




#subset dans tous les sens
psITSgen1 <- subset_samples(psITS, Generation == "Gen1")
psITSgen1
psITSgen1 <- prune_taxa(taxa_sums(psITSgen1) >= 1, psITSgen1)
psITSgen1 #35098 taxa and 54 samples



# Courbe de rarefaction avec limite de l'axe des y à 10^4
RC <- ggrare(psITSgen1, step = 100, se = FALSE, color = "Compartment") + 
  labs(title = "Rarefaction curve ITS, Plants 
       grown under controlled conditions", 
       x = "Depth",
       y = "Number of ASVs") +
  theme_minimal(base_size = 15) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  ) +
      scale_color_manual(values = c( "Root" = "#E3CD3B", "Rhizo" = "#55463A")) +  # Couleurs personnalisées  scale_x_continuous(limITS = c(0, 75000))

scale_x_continuous(limits = c(0, 25000))

# Afficher le graphique
print(RC)
# Enregistrer le graphique
ggsave(RC, dpi=300, device = svg, width = 7, height = 4, filename = "your_path/Supp_1_Rarefaction curve tout par compartiment gen1 ITS.svg")


```




