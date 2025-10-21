# targets-from-pathways

## Contributors
Siddharth Sethi: sidharth.sethi@astx.com

## Codebase

### Data required to load
> **_Note:_** All data files below are just test files with minimal data for our example target genes and disease

- Gene metadata (gene_symbol, ensembl_id)  
  This [file](data/gene_data.txt) serves as a mapping between gene symbols and their corresponding Ensembl IDs. It also helps to validate user input for gene symbols. The genes selected in this file are the top 5 genes associated with male infertility sorted by Chembl score. These genes would be used as input genes for our example.

- Disease metadata (disease_name, disease_efo_id,)  
  This [file](data/disease_data.txt) serves as a mapping between disease names and their corresponding EFO IDs. It also helps to validate user input for disease names.

- Disease association by data source (disease_efo_id, ensembl_id, data_source, association_score)  
  This [file](data/disease_association_data.txt) contains associations between diseases and genes, along with the data source and association score. It is used to fetch genes associated with a specific disease based on user-defined criteria.

- Reactome pathway data from OT (pathway_id, pathway_name, ensembl_id)
- Reactome pathway gene set from REACTOME
- Reactome pathway interactions file

### Required functions
1. validate_input(gene_symbol:str, disease_name:str)
   - Validate and fetch gene ensembl_id and disease EFO_id
2. get_genes_associated_with_disease(disease_efo_id, data_source=None, score_threshold=None)
3. execute_gsea()
4. get_pathways_from_genes(ensembl_ids:List[str])
5. get_genes_from_pathway(pathway_id:List[str])
6. create_gene_matrix(ensembl_ids:List[str], pathway_ids:List[str])
7. parse reactome (?)

### main_script

I have added a main_script.py file that demonstrates how our final script will work. Please populate the remaining pipeline in the main script as required.

```python
python main_script.py --gene_name="SLC22A11" --disease_name="male infertility"
```