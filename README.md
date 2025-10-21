> Open Targets Hackathon, October 21-22, 2025

# targets-from-pathways


## Aim


## Contributors
- Jędrzej Kubica (jedrzej.kubica@univ-grenoble-alpes.fr)
- Siddharth Sethi (sidharth.sethi@astx.com)
- Polina
- Elvis


## Introdution


### TODO Initial plan:
- validate_input(gene_symbol:str, disease_name:str)
   - Validate and fetch gene ensembl_id and disease EFO_id
- get_genes_associated_with_disease(disease_efo_id, data_source=None, score_threshold=None)
- execute_gsea()
- get_pathways_for_gene(ensembl_ids:List[str])
- get_genes_for_pathway(pathway_id:List[str])
- create_gene_matrix(ensembl_ids:List[str], pathway_ids:List[str])
- GSEA
- parse reactome interactions
- find all genes that are both disease-specific and on the same pathways as target
- visualize
- prioritze targets (network propagation?)


## How to use this repo


## Methods

### Flowchart

### Data

Please, put all data files in a directory called data/

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

Reactome gene-to-pathway mapping:
```
wget https://download.reactome.org/94/Ensembl2Reactome_PE_All_Levels.txt
```

Reactome functional interactions as described in Wu et al., 2010:
```
wget http://cpws.reactome.org/caBigR3WebApp2025/FIsInGene_04142025_with_annotations.txt.zip
```


### Scripts

Part 1:

I have added a main_script.py file that demonstrates how our final script will work. Please populate the remaining pipeline in the main script as required.

```
python main_script.py --gene_name="SLC22A11" --disease_name="male infertility"
```

Part 2:

Parse Reactome gene-to-pathway mapping, outputs genes that are both disease-specific and on the same pathways as target ("BTG4", random for now)

```
python run.py --pathway_mapping_file data/Ensembl2Reactome_PE_All_Levels.txt --interactions_file data/FIsInGene_04142025_with_annotations.txt 1>interactions_reactome.tsv
```


## Results


## Future directions


## Python environment

required packages: networkx, pandas


## References
1. Wu, G., Feng, X. & Stein, L. A human functional protein interaction network and its application to cancer data analysis. Genome Biol 11, R53 (2010). https://doi.org/10.1186/gb-2010-11-5-r53