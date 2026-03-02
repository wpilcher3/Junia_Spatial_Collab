# Junia_Spatial_Collab
Misc scripts for analysis

This project uses the following docker images with all required dependencies (except pyreadr) preinstalled
Scanpy/Python: docker.io/wcp7un/mmrf_race_docker_python:v3.13_d120225b
Seurat/R: docker.io/wcp7un/mmrf_race_docker_r:v4.4.1_d110525

The .devcontainer directory includes files for launching an interactive environment inside the above containers using the vscode devcontainer plugin. Update mount paths in .devcontainer/Python/devcontainer.json and .devcontainer/R/devcontainer.json as needed. 

This repository uses some auto generated code to quickly screen and test analysis strategies. Use results at your own peril.

source/Helper_Functions contains various functions for repetitive analyses.
source/Convert_R_Object_To_Scanpy.R takes a Seurat object and converts it to an .h5ad file for anlaysis with Scanpy.
source/Scanpy_DEG_ULM contains functions for computing outcome information from the .h5ad file.
source/R_Visualization contains R scripts to visualize the outputs
source/MMRF_load_data.R contains an R script to load an MMRF Seurat object and integrate the metadata. This will also load all libraries and functions in source/Helper_Functions/. You should update this to point towards an .rds file with all metadata already loaded.
source/MMRF_load_dependencies_only.R is the same as source/MMRF_load_data.R, but only loads the packages and functions, and does not load data.
source/color_and_panel_sizes.R is called from source/MMRF_load_data.R to set some global parameters 

General workflow:

1) Convert Seurat object to a .h5ad file
2) Obtain counts data from bulk RNA along with sample metadata, using the same naming schemes to the MMRF IA metadata objects.
3) Setup parameters in 'papermill_execute...'
4) Run papermill_execute...
