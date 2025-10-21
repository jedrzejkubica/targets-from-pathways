# targets-from-pathways

## Codebase

### Data required to load
- Gene metadata (gene_symbol, ensembl_id)
- Disease metadata (disease_name, disease_efo_id,)
- Disease association by data source (disease_efo_id, ensembl_id, data_source, association_score)
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