



#Figure 3: Heatmap of bacterial and fungal pathogens abundance in healthy and diseased rice  leaves.
```{r}
#----------------------------------------------------
# Heatmap - Pathogen-focused dataset
#----------------------------------------------------

# Load libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)

# Load and filter data
df <- read_excel("your_path/heatmap_que_patho.xlsx") %>%
  filter(somme > 100) %>%
  mutate(
    Genus_species = as.character(Genus_species_disease),
    Kingdom = as.factor(Kingdom),
    Interaction = as.factor(Interaction)
  )

# Select abundance columns and convert to matrix
sample_cols <- c(paste0("H", 1:6), paste0("D", 1:6))
heatmap_matrix <- df %>%
  select(Genus_species, all_of(sample_cols)) %>%
  column_to_rownames("Genus_species") %>%
  as.matrix()

# Create row annotations
row_annot <- df %>%
  select(Genus_species, Kingdom, Interaction, enriched_in) %>%
  distinct(Genus_species, .keep_all = TRUE) %>%
  filter(Genus_species %in% rownames(heatmap_matrix)) %>%
  arrange(match(Genus_species, rownames(heatmap_matrix))) %>%
  mutate(Enrichment = enriched_in) %>%
  column_to_rownames("Genus_species")

# Colors for enrichment
enrichment_colors <- c(
  "Healthy" = "#67B85A",
  "Diseased" = "#E97428",
  "No_diff" = "white"
)

# Reorder rows by interaction category
interaction_order <- c("Pathogen", "Potential_pathogen", "Neutral_beneficial", "Beneficial", "Unknown")
row_annot$Interaction <- factor(row_annot$Interaction, levels = interaction_order)
ordered_taxa <- rownames(row_annot)[order(row_annot$Interaction)]
heatmap_matrix <- heatmap_matrix[ordered_taxa, ]
row_annot <- row_annot[ordered_taxa, ]

# Convert abundances to relative values within each kingdom per sample
kingdom_vec <- row_annot$Kingdom
heatmap_matrix_rel <- matrix(NA, nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix),
                             dimnames = dimnames(heatmap_matrix))

for (sample in colnames(heatmap_matrix)) {
  for (king in unique(kingdom_vec)) {
    taxa_idx <- which(kingdom_vec == king)
    col_sum <- sum(heatmap_matrix[taxa_idx, sample])
    if (col_sum > 0) {
      heatmap_matrix_rel[taxa_idx, sample] <- (heatmap_matrix[taxa_idx, sample] / col_sum) * 100
    } else {
      heatmap_matrix_rel[taxa_idx, sample] <- 0
    }
  }
}
heatmap_matrix <- heatmap_matrix_rel

# Define annotation colors
kingdom_colors <- c("Fungi" = "darkgreen", "Bacteria" = "darkblue")
interaction_colors <- c(
  "Beneficial" = "forestgreen",
  "Unknown" = "gray50",
  "Potential_pathogen" = "blue",
  "Neutral_beneficial" = "white",
  "Pathogen" = "red"
)
row_annot$Interaction <- as.character(row_annot$Interaction)

# Create heatmap
topPathoHeatmap <- Heatmap(
  heatmap_matrix,
  name = "Abundance",
  col = colorRamp2(c(0, 1, 2, 5, 10, 20, 50, 100),
                   c("gray", "#FFEDA0", "#FFD076", "#FEB24C",
                     "#F77736", "#F03B20", "#C6321B", "#9B2816")),
  right_annotation = rowAnnotation(
    Kingdom = row_annot$Kingdom,
    Interaction = row_annot$Interaction,
    Enrichment = row_annot$Enrichment,
    col = list(Kingdom = kingdom_colors,
               Interaction = interaction_colors,
               Enrichment = enrichment_colors),
    gp = gpar(col = "white", lwd = 0.5)
  ),
  column_order = sample_cols,
  row_order = rownames(heatmap_matrix),
  row_names_gp = gpar(
    fontface = ifelse(row_annot$Interaction == "Pathogen", "bold.italic", "italic"),
    fontsize = 10
  ),
  row_split = row_annot$Kingdom,
  column_names_rot = 45,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_title = NULL,
  border = TRUE,
  rect_gp = gpar(col = "gray", lwd = 0.5)
)

# Display heatmap
topPathoHeatmap

```


#Supplementary Figure S2: Heatmap showing the abundance of bacterial and fungal taxa in healthy 567 and diseased rice leaves.
```{r}
#----------------------------------------------------
# Heatmap - Full dataset
#----------------------------------------------------

# Load libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)

# Load and filter data
df <- read_excel("your_path/heatmap_tout.xlsx") %>%
  filter(somme > 100) %>%
  mutate(
    Genus_species = as.character(Genus_species),
    Kingdom = as.factor(Kingdom),
    Interaction = as.factor(Interaction)
  )

# Select abundance columns and convert to matrix
sample_cols <- c(paste0("H", 1:6), paste0("D", 1:6))
heatmap_matrix <- df %>%
  select(Genus_species, all_of(sample_cols)) %>%
  column_to_rownames("Genus_species") %>%
  as.matrix()

# Create row annotations
row_annot <- df %>%
  select(Genus_species, Kingdom, Interaction, enriched_in) %>%
  distinct(Genus_species, .keep_all = TRUE) %>%
  filter(Genus_species %in% rownames(heatmap_matrix)) %>%
  arrange(match(Genus_species, rownames(heatmap_matrix))) %>%
  mutate(Enrichment = enriched_in) %>%
  column_to_rownames("Genus_species")

# Colors for enrichment
enrichment_colors <- c(
  "Healthy" = "#67B85A",
  "Diseased" = "#E97428",
  "No_diff" = "white"
)

# Reorder rows by interaction category
interaction_order <- c("Pathogen", "Potential_pathogen", "Neutral_beneficial", "Beneficial", "Unknown")
row_annot$Interaction <- factor(row_annot$Interaction, levels = interaction_order)
ordered_taxa <- rownames(row_annot)[order(row_annot$Interaction)]
heatmap_matrix <- heatmap_matrix[ordered_taxa, ]
row_annot <- row_annot[ordered_taxa, ]

# Convert abundances to relative values within each kingdom per sample
kingdom_vec <- row_annot$Kingdom
heatmap_matrix_rel <- matrix(NA, nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix),
                             dimnames = dimnames(heatmap_matrix))

for (sample in colnames(heatmap_matrix)) {
  for (king in unique(kingdom_vec)) {
    taxa_idx <- which(kingdom_vec == king)
    col_sum <- sum(heatmap_matrix[taxa_idx, sample])
    if (col_sum > 0) {
      heatmap_matrix_rel[taxa_idx, sample] <- (heatmap_matrix[taxa_idx, sample] / col_sum) * 100
    } else {
      heatmap_matrix_rel[taxa_idx, sample] <- 0
    }
  }
}
heatmap_matrix <- heatmap_matrix_rel

# Define annotation colors
kingdom_colors <- c("Fungi" = "darkgreen", "Bacteria" = "darkblue")
interaction_colors <- c(
  "Beneficial" = "forestgreen",
  "Unknown" = "#D8D9D3",
  "Potential_pathogen" = "blue",
  "Neutral_beneficial" = "#BDCB90",
  "Pathogen" = "red"
)
row_annot$Interaction <- as.character(row_annot$Interaction)

# Create heatmap
topPathoHeatmap <- Heatmap(
  heatmap_matrix,
  name = "Abundance",
  col = colorRamp2(c(0, 1, 2, 5, 10, 20, 50, 100),
                   c("gray", "#FFEDA0", "#FFD076", "#FEB24C",
                     "#F77736", "#F03B20", "#C6321B", "#9B2816")),
  right_annotation = rowAnnotation(
    Kingdom = row_annot$Kingdom,
    Interaction = row_annot$Interaction,
    Enrichment = row_annot$Enrichment,
    col = list(Kingdom = kingdom_colors,
               Interaction = interaction_colors,
               Enrichment = enrichment_colors),
    gp = gpar(col = "white", lwd = 0.5)
  ),
  column_order = sample_cols,
  row_order = rownames(heatmap_matrix),
  row_names_gp = gpar(
    fontface = ifelse(row_annot$Interaction == "Pathogen", "bold.italic", "italic"),
    fontsize = 10
  ),
  row_split = row_annot$Kingdom,
  column_names_rot = 35,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_title = NULL,
  border = TRUE,
  rect_gp = gpar(col = "gray", lwd = 0.5)
)

# Display heatmap
topPathoHeatmap

```

