import papermill as pm
import os
sub_directory = f'Junia_Spatial_Analyses/'

output_dir = "output/" + sub_directory

CUSTOM_PATHWAYS = {
    # "PC_1": ["NUSAP1", "H3C8", "PRIM2", "TOP2A", "CDK1", "AURKB", "PLK1", "H3C10", "HMGB2", "H4C3", "TPX2", "H2BC5", "KPNA2", "LIG1", "TK1", "STMN1", "UBE2S", "TYMS", "LMNB1", "EZH2", "UBE2T", "HMGN2", "TUBA1B", "CKS1B", "TUBB4B", "TUBB", "H2AZ1", "TUBA1C", "H2AZ1"],
    "PC_1_UP_2020A": ["HIST1H1B3", "NUSAP1", "HIST1H3G", "HMGN2", "TOP2A", "CDK1", "AURKB", "PLK1", "HIST1H3H", "HMGB2", "HIST1H4C6", "TPX2", "HIST1H2BD", "KPNA2", "LIG1", "TK1", "STMN1", "UBE2S", "TYMS", "LMNB1", "EZH2", "UBE2T", "HMGN2", "TUBA1B", "CKS1B", "TUBB4B", "TUBB", "H2AFZ", "TUBA1C", "H2AZ1"],
    # "PC_2_UP": ["NEAT1", "NLRP1", "OGT", "PTEN", "DISC1", "SNRPD1", "CDH16", "RAS4A", "C9JAP3", "NDUFA3", "PTK2", "IRF3", "POP6", "VPS29", "SIGMAR1", "JAK3", "EI-BP1L1", "OXR1", "SDHA"],
    "PC_2_UP_2020A": ["NEAT1", "NLRP1", "OGT", "PTEN", "DISC1", "SNRPD1", "CDH16", "RAS4A", "C9JAP3", "NDUFA3", "PTK2", "IRF3", "POP6", "VPS29", "SIGMAR1", "JAK3", "EIBP1L1", "OXR1", "SDHA"],
    # "PC_2_Down": ["H3C2", "H3C8", "H3C10", "H4C3", "STMN1", "HMGB2", "TUBA1B", "TUBB", "H2AZ1", "PTG82", "DLL1", "PTK2", "UBE2S", "CHMP3", "CHMP3", "NDUFA3", "RAS4A", "SNRPD1", "SIGRP2", "IER5", "PTK2", "UBE2S", "SHSF10", "VPS29", "OGT", "DDX39A", "NEAT1", "PRPF3", "TUBB4B", "MDM4", "MAU2", "CDKN1B", "SIGMAR1", "SDHA", "UNK1", "OXR1"],
    "PC_2_DOWN_2020A": ["HIST1H1B3", "HIST1H3G", "HIST1H3H", "HIST1H4C6", "STMN1", "HMGB2", "TUBA1B", "TUBB", "H2AFZ", "PTG82", "DLL1", "PTK2", "UBE2S", "CHMP3", "CHMP3", "NDUFA3", "RAS4A", "SNRPD1", "SIGRP2", "IER5", "PTK2", "UBE2S", "SHSF10", "VPS29", "OGT", "DDX39A", "NEAT1", "PRPF3", "TUBB4B", "MDM4", "MAU2", "CDKN1B", "SIGMAR1", "SDHA", "UNK1", "OXR1"],
    # "PC_3_DOWN": ["H3C2", "H3C8", "H2BC5", "H3C10", "H4C3", "HMGB2", "DISC1", "PTEN", "SNRPD1", "SIGRP2", "IER5"],
    "PC_3_DOWN_2020A": ["HIST1H1B3", "HIST1H3G", "HIST1H2BD", "HIST1H3H", "HIST1H4C6", "HMGB2", "DISC1", "PTEN", "SNRPD1", "SIGRP2", "IER5"]
}


# High vs Low AFR, across Plasma cells
# Settings
# cluster_list = ["NK_adp", "NK_CD56dim", "NK_resident", "NK_CD56bright"]
file_name = "Plasma_Cells_scRNA"
subdir_name = "Plasma_Cells"
analysis_description = "Analysis of scRNA-seq plasma cells"
compartment_description = "Plasma Cells"
subset_mode = "DICT"
subset_dict = dict(
    lineage_group=dict(include=["P"], exclude=None)
)
# Execute
os.makedirs(output_dir + subdir_name + '/', exist_ok=True)
NOTEBOOK_OUTPUT = f'{output_dir}/{subdir_name}/{file_name}_REPORT.ipynb'
pm.execute_notebook(
    'source/Scanpy_ULM/_template_DEG_scRNA.ipynb',
    NOTEBOOK_OUTPUT,
    parameters=dict(
        CLUSTERS_ANALYSIS=None,
        FILE_NAME=file_name,
        ADDITIONAL_SUBDIR=subdir_name,
        ANALYSIS_DESC=analysis_description,
        COMPARTMENT_DESCRIPTION=compartment_description,
        SUBSET_MODE=subset_mode,
        SUBSET_MODE_CUSTOM_DICT=subset_dict,
        CUSTOM_PATHWAYS=CUSTOM_PATHWAYS
    )
)
os.system(
    f'jupyter nbconvert --output-dir={output_dir + subdir_name + "/"} --to html {NOTEBOOK_OUTPUT}')
