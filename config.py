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

PROTEIN_PDB     = "1bcu-protein.pdb"
PROTEIN_NAME    = "1bcu-protein"
PROTEIN_CHAINS  = "HL"
BOX_SIZE        = 20
EXHAUSTIVENESS  = 8
VINA_SEED       = 42
VINA_CPU        = 0
FIXED_CENTER    = (9.558, 20.381, 50.273)

LIGAND_NAME = "1bcu-ligand"

RAW_DATA_DIR = "1.Raw_data"
RESULT_DIR   = "2.Result"

DOCKING_DIR   = "Docking_Results"
DOCKING_TEMP  = "Docking_Temp"
MD_SETUP_DIR  = "MD_Setup"
ANALYSIS_DIR  = "Analysis"

LLM_BASE_URL = "https://api.siliconflow.cn/v1"
LLM_MODEL    = "Qwen/Qwen2.5-7B-Instruct"
LLM_API_KEY_ENV = "DASHSCOPE_API_KEY"