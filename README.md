> Open Targets Hackathon, October 21-22, 2025


<h1 align="left">Project #18: Target prioritisation by pathways </h1>

<p align="center">
  <img width="600" height="400" src="figures/team_logo.png">
</p>


The aim of the project was to develop a disease-specific pathway assessment tool. We built a bioinformatics workflow that uses an associated target and disease-specific biological pathways to find novel targets implicated in disease. This tool can help researchers understand the biological pathways underlying disease and identify novel intervention points considering both efficacy and safety profiles along with target tractability. Here we present a prototype developed during the hackathon.


## Contributors
- Jędrzej Kubica (jedrzej.kubica@univ-grenoble-alpes.fr)
- Polina Rusina (polina@ebi.ac.uk)
- Siddharth Sethi (sidharth.sethi@astx.com)
- Elvis Poku-Adusei (elvispokkad@yahoo.com)


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
- find all genes that are both disease-specific and on the same pathways as target
- prioritze targets (network propagation? anothe scoring formula?)
- parse reactome interactions and map new targets to pathways, visualize


## How to use this repo

```
git clone git@github.com:jedrzejkubica/targets-from-pathways.git
cd targets-from-pathways
```

Prepare input data and run the pipeline as described below. See results in scores.tsv.


## Methods

We introduced a new methodology for target prioritization. The user inputs a disease and a target of interest (perhaps the medication for this target doesn't work). First, we get a gene list for the disease from Open Targets Platform (Associations - indirect (by data source) ; http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/association_by_datasource_indirect) and perform GSEA using the blitzgsea Python package (https://github.com/MaayanLab/blitzgsea). We obtain pathways associated with the disease. Second, we find also all pathways with the target on it also using Reactome database (gene-to-pathway mapping ; https://download.reactome.org/94/Ensembl2Reactome_PE_All_Levels.txt). That's how we get pathways associated with the disease and the target and we find all genes that are on these pathways. The, we prioritize these genes as new potential targets for the disease by:
- scoring using pathway selectivity (formula below)
- scoring using network propagation (Random Walk with Restart, MutliXrank (Baptista et al, 2020) ; https://github.com/anthbapt/multixrank)
Finally, we obtain a list of scored genes, the higher the score, the more likely it is to be associated with the disease and the target (based on pathways only).

### Flowchart

### Data

Please, put all data files in a directory called `data/`

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

Reactome gene-to-pathway mapping file:
```
wget https://download.reactome.org/94/Ensembl2Reactome_PE_All_Levels.txt
```

Reactome functional interaction file (see Wu et al., 2010) to build a functional interaction network:
```
wget http://cpws.reactome.org/caBigR3WebApp2025/FIsInGene_04142025_with_annotations.txt.zip
```


### Run pipeline

Part 1:

Use Pathway Targets Extractor (https://github.com/MaayanLab/blitzgsea), a Python module for Gene Set Enrichment Analysis, to extract pathway-target gene pairs.


Part 2:

main_script.py demonstrates how our final script should work. Please populate the remaining pipeline in the main script as required.

```
python main_script.py --gene_name="SLC22A11" --disease_name="male infertility"
```

Part 3:

Parse Reactome gene-to-pathway mapping, outputs genes that are both disease-specific and on the same pathways as target (using "BTG4" for now as example)

```
python run.py --pathway_mapping_file data/Ensembl2Reactome_PE_All_Levels.txt --interactions_file data/FIsInGene_04142025_with_annotations.txt --target BTG4 1>scores.tsv
```

Part 4:

For every gene found in the overlap between disease-specific and target-specific list:

`score = (number of disease pathways containing the gene + number of target pathways containing the gene) / (total disease pathways + total target pathways)`

Network propagation on the functional interaction network:

use run.py (uncomment lines in main() if necessary) to generate interactions_reactome.tsv

Run RWR:
```
prepare config.yml
mkdir -p reactome/multiplex/1
cp ~/open-targets-hackathon/targets-from-pathways/results/interactions_reactome.tsv
prepare seeds.txt (target + disease genes same as for GSEA)
python run_multirank.py (provided in this repo)
```

## Results

Here are top 10 genes based on pathway selectivity (using BTG4 as target as example):

| Gene   | Score |
|--------|--------|
| TUBB4B | 0.76   |
| RPS2   | 0.72   |
| POLR2D | 0.68   |
| PABPC1 | 0.68   |
| EIF4G1 | 0.68   |
| EIF4E  | 0.68   |
| EIF4B  | 0.68   |
| EIF4A2 | 0.68   |
| EIF4A1 | 0.68   |
| BIRC5  | 0.68   |
| ZP2    | 0.64   |


Reactome functional interaction network:
![interactions_reactome](figures/interactions_reactome.png)

Here are top 10 genes based on network propagation (using target + disease genes as seeds):


## Future directions

Short-term:
- refine and test the pathway selectivity scoring formula
- finalize and assess network propagation with Reactome functional interaction network (interactions_reactome.tsv)
- get insights about scoring of new targets + interpretability
- assess the method's performance, compare with similar methodologies

Long-term:
- adapt for other target types (genes, proteins, miRNA)
- identify drug repurposing opportunities
- effect of target knockout within the network
- given a desired effect (inhibition, activation) recommend top drugs
- use efficy and safety (Open Targets Pharmacovigilance)
- integrate pipeline within Open Targets Platform (https://platform.opentargets.org/)


## Python environment

required packages:
- networkx
- pandas
- blitzgsea (https://github.com/MaayanLab/blitzgsea)
- multixrank (https://github.com/anthbapt/multixrank)


## Special thank you to the Organizers of the Open Targets Hackathon!


## References
1. Wu, G., Feng, X. & Stein, L. A human functional protein interaction network and its application to cancer data analysis. Genome Biol 11, R53 (2010). https://doi.org/10.1186/gb-2010-11-5-r53
2. Baptista, A., Gonzalez, A. & Baudot, A. Universal multilayer network exploration by random walk with restart. Commun Phys 5, 170 (2022). https://doi.org/10.1038/s42005-022-00937-9