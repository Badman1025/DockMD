import sys

IS_WINDOWS = sys.platform.startswith("win")

if IS_WINDOWS:
    TOOLS = {
        "pymol":  "pymol",
        "obabel": "obabel",
        "vina":   "vina",
        "gmx":    "gmx",
        "acpype": "acpype",
    }
else:
    TOOLS = {
        "pymol":  "pymol",
        "obabel": "obabel",
        "vina":   "vina",
        "gmx":    "gmx",
        "acpype": "acpype",
    }

PROTEIN_PDB     = "6GU6-CDK1-A1B1.pdb"
PROTEIN_NAME    = "CDK1"
PROTEIN_CHAINS  = "AB"
BOX_SIZE        = 20
EXHAUSTIVENESS  = 8
VINA_SEED       = 42
VINA_CPU        = 0
FIXED_CENTER    = None

LIGAND_NAME = "JAK-Inhibitor-I"

RAW_DATA_DIR = "1.Raw_data"
RESULT_DIR   = "2.Result"

DOCKING_DIR   = "Docking_Results"
DOCKING_TEMP  = "Docking_Temp"
MD_SETUP_DIR  = "MD_Setup"
ANALYSIS_DIR  = "Analysis"
