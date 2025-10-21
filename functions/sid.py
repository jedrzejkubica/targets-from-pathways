import pandas as pd
from collections import defaultdict

GENE_MAP = None
DISEASE_MAP = None
DISEASE_ASSOC_DATA = None

def load_data() -> None:
    """Load necessary data files."""
    gene_data_file = "data/gene_data.txt"
    disease_data_file = "data/disease_data.txt"
    disease_assoc_file = "data/input/auto-input/filtered_associations.tsv"

    global GENE_MAP, DISEASE_MAP, DISEASE_ASSOC_DATA
    DISEASE_ASSOC_DATA = pd.read_csv(disease_assoc_file, sep="\t")
    GENE_MAP = _load_data_into_map(gene_data_file)
    DISEASE_MAP = _load_data_into_map(disease_data_file)


def _load_data_into_map(fpath:str) -> dict:
    """Loads a data file into a mapping dictionary.
    """
    idmap = defaultdict(list)
    df = pd.read_csv(fpath, sep="\t", dtype=str)
    for row in df.itertuples():
        key, value = row[1], row[2]
        if value not in idmap[key]:
            idmap[key] = value
    return idmap


def validate_input(gene_name:str, disease_name:str) -> tuple:
    """
    Validates the input gene name and disease name against predefined mappings.
    Args:
        gene_name (str): The name of the gene to validate.
        disease_name (str): The name of the disease to validate.
    Returns:
        tuple: A tuple containing the gene ID and disease ID if both are valid.
    Raises:
        ValueError: If the gene name is not found in the GENE_MAP.
        ValueError: If the disease name is not found in the DISEASE_MAP.
    """
    if gene_name not in GENE_MAP:
        raise ValueError(f"Gene name '{gene_name}' not found in GENE_MAP.")
    if disease_name not in DISEASE_MAP:
        raise ValueError(f"Disease name '{disease_name}' not found in DISEASE_MAP.")
    
    gene_id = GENE_MAP[gene_name]
    disease_id = DISEASE_MAP[disease_name]
    return gene_id, disease_id


def get_genes_associated_with_disease(disease_id:str) -> pd.DataFrame:
    """
    Retrieve genes associated with a specific disease based on its ID.
    Args:
        disease_id (str): The unique identifier of the disease.
    Returns:
        pd.DataFrame: A DataFrame containing the `targetId` and `score` columns
        for genes associated with the specified disease, sorted by `score` in
        descending order.
    """
    df = DISEASE_ASSOC_DATA[DISEASE_ASSOC_DATA['diseaseId'] == disease_id]
    df_selected = df[['targetId', 'score']].sort_values(by='score', ascending=False)
    return df_selected
