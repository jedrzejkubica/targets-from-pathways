import argparse
from functions import sid

# command line options
parser = argparse.ArgumentParser(description="Identifies secondary targets associated with a given disease on the same pathway as the primary target gene.")
parser.add_argument('--gene_name', type=str, required=True, help='Primary target name')
parser.add_argument('--disease_name', type=str, required=True, help='Didease indication space for primary target')

args = parser.parse_args()

# load all required data
sid.load_data()

# validate input gene and disease names
gene_id, disease_id = sid.validate_input(args.gene_name, args.disease_name)

# get genes associated with the disease
disease_genes_df = sid.get_genes_associated_with_disease(disease_id) # the output of this would be the input to gsea

# genes go into gsea

