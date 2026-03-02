# Outputs to data/h5ad/Combined_Full_Object/Combined_Full_Object.h5ad

source("source/MMRF_load_data.R", local = T)

data_dir_out <- file.path(dataPath, "h5ad", "Combined_Full_Object")
data_file <- file.path(data_dir_out, "Combined_Full_Object.h5ad")


dir.create(data_dir_out, recursive = TRUE)
source("./source/Race/Trajectory/prepare_h5ad.r", local = T)
convert_MMRF_to_h5ad(merged, output_path = data_file)
