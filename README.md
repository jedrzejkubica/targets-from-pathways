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

perhaps the medication for this target doesn't work

Initial plan:
- get genes associated with disease
- get disease pathways using GSEA
- get genes for disease pathways
- get genes for target pathways
- find all genes that are both disease-specific and on the same pathways as target
- prioritze targets (network propagation on reactome functional interaction network? anothe scoring formula?)
- visualize pathways 


## How to use this repo

```
git clone git@github.com:jedrzejkubica/targets-from-pathways.git
cd targets-from-pathways
```

Prepare input data and run the pipeline as described below. See results in `scores.tsv`.


## Methods

We introduced a new methodology for target prioritization. The user inputs a disease and a target of interest. First, we get a gene list for the disease from Open Targets Platform (https://platform.opentargets.org/) and perform fast Gene Set Enrichment Analysis (GSEA) using the blitzgsea Python package (https://github.com/MaayanLab/blitzgsea). We obtain pathways associated with the disease. Second, we find all pathways with the target on it also using Reactome database. That's how we get pathways associated with the disease and the target, and we find all genes that are on these pathways. Then, we prioritize these genes as new potential targets for the disease using pathway selectivity and network propagation (Random Walk with Restart, MutliXrank (Baptista et al, 2020) ; https://github.com/anthbapt/multixrank). Finally, we obtain a list of scored genes, the higher the score, the more likely it is to be associated with the disease and the target (based on pathways and network propagation).


### Flowchart
<img width="901" height="261" alt="Flow chart drawio" src="https://github.com/user-attachments/assets/df33bbc8-1842-4937-9062-806452cb6636" />


### Data

Create data/ folder
```
mkdir data
cd data/
```

Download:
- Open Targets associations parquets, source: Open Targets Platform (Associations - indirect (by data source))
  ```
  wget --recursive --no-parent --no-host-directories --cut-dirs 6 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/association_by_datasource_indirect .
  ```

- Open Targets targets parquets, source: Open Targets Platform (Target)
  ```
  wget --recursive --no-parent --no-host-directories --cut-dirs 6 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/target .
  ```

- ReactomePathways.gmt file, source: https://reactome.org/download/current/
  ```
  wget https://reactome.org/download/current/ReactomePathways.gmt.zip
  unzip ReactomePathways.gmt.zip
  ```

- Reactome gene-to-pathway mapping file, source: https://reactome.org/download-data
  ```
  wget https://download.reactome.org/94/Ensembl2Reactome_PE_All_Levels.txt
  ```
  
- Reactome functional interactions (described in Wu et al., 2010), source: https://reactome.org/download-data
  ```
  wget http://cpws.reactome.org/caBigR3WebApp2025/FIsInGene_04142025_with_annotations.txt.zip
  unzip FIsInGene_04142025_with_annotations.txt.zip
  ```


### Run pipeline

During the hackathon we focused on male infertiltiy (EFO_0004248) and ESR1 as target of interest. 

1) GSEA:
```
python run_gsea.py --target_parquets_dir ../data/target/ --associations_parquets_dir ../data/association_by_datasource_indirect/ --disease EFO_0004248 --datatype genetic_association --gmt_file ../data/ReactomePathways_merged.gmt --pval_threshold 0.05 --fdr_threshold 1 1>enriched_pathways.txt
```

2) Pathway-based scoring:
```
python run_reactome.py --pathway_mapping_file ../data/Ensembl2Reactome_PE_All_Levels.txt --disease_pathways_file ../gsea/disease_pathways.txt --interactions_file ../data/FIsInGene_04142025_with_annotations.txt --target ESR1 1>scores.tsv
```

3) Network propagation scoring:

TODO run multixrank


## Results

Here are top 10 genes based on pathway selectivity (using ESR1 as target as example):

| GENE          | SCORE      |
|---------------|------------|
| ESR1          | 0.84       |
| PGR           | 0.5333333333333333 |
| UBC(77-152)   | 0.4666666666666667 |
| UBC(609-684)  | 0.4666666666666667 |
| UBC(533-608)  | 0.4666666666666667 |
| UBC(457-532)  | 0.4666666666666667 |
| UBC(381-456)  | 0.4666666666666667 |
| UBC(305-380)  | 0.4666666666666667 |
| UBC(229-304)  | 0.4666666666666667 |


Reactome functional interaction network:
![interactions_reactome](figures/interactions_reactome.png)

Here are top 10 genes based on network propagation (using ESR1 + disease genes as seeds):



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


### Python environment

required packages: networkx, pandas, blitzgsea (https://github.com/MaayanLab/blitzgsea), multixrank (https://github.com/anthbapt/multixrank)


## Special thank you to the Organizers of the Open Targets Hackathon!


## References
1. Wu, G., Feng, X. & Stein, L. A human functional protein interaction network and its application to cancer data analysis. Genome Biol 11, R53 (2010). https://doi.org/10.1186/gb-2010-11-5-r53
2. Baptista, A., Gonzalez, A. & Baudot, A. Universal multilayer network exploration by random walk with restart. Commun Phys 5, 170 (2022). https://doi.org/10.1038/s42005-022-00937-9
