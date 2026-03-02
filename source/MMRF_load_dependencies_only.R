# Figure 3, Setup and Data Loading

library(kableExtra)
library(Seurat)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 600 * 1024^3)

# setwd("/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/")

source("source/Helper_Functions/general.R", chdir = T)
MEGASEQ <- T

output <- file.path("output", "Junia_Spatial_Analyses")

dir.create(output, recursive = T)
# objectPath <- "INTEGRATED_OBJECTS_MMRF/Integration/SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds"
# objectPath <- "objects_stripped/Full/FULL_NO_MD_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds"
newMD <- file.path("data", "metadata_032525/per_cell_metadata/032525_per_cell_md.rds") # "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data//per_cell_md.rds"
newMD <- readRDS(newMD)

# allDimReducs <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata/subumap_dimreduc.rds"

# objectPath <- file.path(dataPath, "VALIDATION_COHORT/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds")



clinical_group <- "progression_group"

reg_auc <- file.path("data", "SCENIC_FILES/auc_mtx.csv")

# output_F2 <- file.path(output, "FIG2_OLD")
# dir.create(output_F2, recursive = T)

# output_F3 <- file.path(output, "FIG3")
# dir.create(output_F3, recursive = T, showWarnings = F)

cluster_naming <- read.csv(file.path("data", "custom_MD/annotations_plotting_info_010824.csv"), header = T)
lineage_label <- data.frame(lineage_group = c("CD4", "CD8", "B", "M", "Nk", "P", "E", "Other", "LQ"), lineage_name = c("CD4+ T Cell", "CD8+ T Cell", "B Cell", "Myeloid", "Nk Cell", "Plasma", "Erythroid", "Other", "LQ/Doublet"))
cluster_naming <- dplyr::left_join(cluster_naming, lineage_label, by = "lineage_group")

# merged <- readRDS(file.path(objectPath))
# merged <- merged %>% fix_public_id()
# merged[["umap.sub"]] <- readRDS(allDimReducs) # Adds umap embeddings for all subclusters for display
# merged[["umap.sub"]] <- readRDS("/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata/validation_subumap_dimreduc.rds")
# newMD <- readRDS(newMD)
# merged <- AddMetaData(merged, newMD)

# merged@meta.data$LABEL_TRANSFER_subcluster_V03072023 <- merged@meta.data$seurat_subclusters_label_transferring_Yizhe_v1

# cluster_naming_mod <- cluster_naming
# colnames(cluster_naming_mod) <- paste0("LABEL_TRANSFER_", colnames(cluster_naming_mod))
# merged@meta.data <- merged@meta.data |>
#     dplyr::left_join(cluster_naming_mod, by = "LABEL_TRANSFER_subcluster_V03072023", suffix = c("", ".y")) |>
#     dplyr::select(-ends_with(".y")) # remove duplicates, keep existing

# merged@meta.data$INTEGRATION_subcluster_V03072023 <- merged@meta.data$seurat_subclusters_integration_Yizhe_v1
# cluster_naming_mod <- cluster_naming
# colnames(cluster_naming_mod) <- paste0("INTEGRATION_", colnames(cluster_naming_mod))
# merged@meta.data <- merged@meta.data |>
#     dplyr::left_join(cluster_naming_mod, by = "INTEGRATION_subcluster_V03072023", suffix = c("", ".y")) |>
#     dplyr::select(-ends_with(".y")) # remove duplicates, keep existing



# merged@meta.data$subcluster_V03072023 <- merged@meta.data$INTEGRATION_subcluster_V03072023
# merged@meta.data$lineage_group <- merged@meta.data$INTEGRATION_lineage_group
# merged@meta.data$lineage_order <- merged@meta.data$INTEGRATION_lineage_order
# merged@meta.data$color <- merged@meta.data$INTEGRATION_color
# merged@meta.data$color_sub <- merged@meta.data$INTEGRATION_color_sub
# merged@meta.data$primary_id_color <- merged@meta.data$INTEGRATION_primary_id_color

# merged@meta.data$subcluster_V03072023_compartment <- merged@meta.data$compartment

# merged@meta.data$cellID_short <- merged@meta.data$INTEGRATION_cellID_short
# merged@meta.data$cellID_long <- merged@meta.data$INTEGRATION_cellID_long

# rownames(merged@meta.data) <- merged$cellname
# merged$display_name <- merged$cellID_short

# cluster_group <- "subcluster_V03072023"
# merged$seurat_cluster_compartment <- merged$subcluster_V03072023_compartment

# merged@meta.data$seurat_cluster_subcluster <- merged@meta.data[, cluster_group]
# merged@meta.data <- merged@meta.data %>% mutate(treatment_simple = dplyr::case_when(
#     d_tx_induction_cat %in% c("imid_pi_steroid", "chemo_imid_pi_steroid") ~ "Triplet",
#     .default = "Doublet"
# ))

# merged@meta.data <- merged@meta.data %>% mutate(isDoublet = dplyr::case_when(
#     .data[["INTEGRATION_doublet_pred"]] %in% c("dblet_cluster", "poss_dblet_cluster") ~ "Doublet",
#     .default = "Singlet"
# ))

# merged@meta.data$lineage_group <- factor(merged@meta.data$lineage_group, levels = c("CD4", "CD8", "B", "M", "Nk", "P", "E", "Other", "LQ"))


source("source/Race/colors_and_panel_sizes.R", local = T)
