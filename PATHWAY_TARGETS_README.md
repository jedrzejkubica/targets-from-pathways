# Pathway Targets Extractor

A Python module for running GSEA (Gene Set Enrichment Analysis) and extracting pathway-target gene pairs in a simplified format.

## Features

- **Run GSEA Analysis**: Execute GSEA using blitzgsea with custom GMT pathway files
- **Extract Pathway-Targets**: Generate tab-separated output with pathways and their associated target genes
- **Flexible Input**: Accept pandas DataFrames or file paths (TSV format)
- **Simple Output**: Two-column format (pathways, targets) for easy downstream analysis

## Installation

Ensure you have the required dependencies:

```bash
pip install pandas blitzgsea
```

## Input Format

The input file/DataFrame must contain:
- `symbol`: Gene symbols (matching those in your GMT file)
- `globalScore`: Numerical scores for ranking genes

Example:
```
symbol	globalScore
CDKN2B	8.694865
ABO	8.224087
FTO	7.753310
SH2B3	6.811755
APOE	6.340978
```

## Output Format

The output is a tab-separated file with two columns:
- `pathways`: Pathway name
- `targets`: Comma-separated list of target genes from the propagated pathway set

Example:
```
pathways	targets
Apoptosis	TP53,BCL2,CASP3,BAX,CASP9,...
Cell Cycle	CDK1,CCNB1,CDC20,PLK1,BUB1,...
```

## Usage

### 1. Run GSEA and Extract Targets (Main Use Case)

```python
from pathway_targets_extractor import run_gsea_and_extract_targets

# From a file
result = run_gsea_and_extract_targets(
    input_df='scores.tsv',
    gmt_file='ReactomePathways.gmt',
    processes=4
)

# Save to file
result.to_csv('pathway_targets.tsv', sep='\t', index=False)
```

### 2. Using a DataFrame

```python
import pandas as pd
from pathway_targets_extractor import run_gsea_and_extract_targets

# Create or load your DataFrame
scores_df = pd.DataFrame({
    'symbol': ['CDKN2B', 'ABO', 'FTO'],
    'globalScore': [8.69, 8.22, 7.75]
})

result = run_gsea_and_extract_targets(
    input_df=scores_df,
    gmt_file='ReactomePathways.gmt'
)
```

### 3. Extract from Existing GSEA Results

If you already have GSEA results with `Term` and `propagated_edge` columns:

```python
from pathway_targets_extractor import extract_pathway_targets

result = extract_pathway_targets('existing_gsea_results.tsv')
result.to_csv('pathway_targets.tsv', sep='\t', index=False)
```

### 4. Direct File-to-File Conversion

```python
from pathway_targets_extractor import extract_pathway_targets_to_file

extract_pathway_targets_to_file(
    input_path='gsea_results.tsv',
    output_path='pathway_targets.tsv'
)
```

## Functions

### `run_gsea_and_extract_targets(input_df, gmt_file, processes=4)`

Main function that runs GSEA and extracts pathway-target pairs.

**Parameters:**
- `input_df`: DataFrame or path to TSV file with 'symbol' and 'globalScore' columns
- `gmt_file`: Path to GMT file with pathway definitions
- `processes`: Number of parallel processes for GSEA (default: 4)

**Returns:** DataFrame with 'pathways' and 'targets' columns

### `extract_pathway_targets(input_df, pathway_col="Term", targets_col="propagated_edge")`

Extracts pathway-target pairs from existing GSEA results.

**Parameters:**
- `input_df`: DataFrame or path to TSV file with GSEA results
- `pathway_col`: Column name for pathways (default: "Term")
- `targets_col`: Column name for targets (default: "propagated_edge")

**Returns:** DataFrame with 'pathways' and 'targets' columns

### `load_custom_gmt(path)`

Utility function to parse GMT files.

**Parameters:**
- `path`: Path to GMT file

**Returns:** Dictionary mapping pathway names to lists of genes

## GMT File Format

GMT (Gene Matrix Transposed) files should have the format:
```
PathwayName	Description	Gene1	Gene2	Gene3	...
```

Example Reactome GMT file:
```
Apoptosis [R-HSA-109581]	Apoptosis	TP53	BCL2	CASP3	BAX	...
Cell Cycle [R-HSA-1640170]	Cell Cycle	CDK1	CCNB1	CDC20	...
```

## Example Workflow

```python
from pathway_targets_extractor import run_gsea_and_extract_targets
import pandas as pd

# 1. Run GSEA
result_df = run_gsea_and_extract_targets(
    input_df='geneset_disease_zscore.tsv',
    gmt_file='ReactomePathways_merged.gmt',
    processes=4
)

# 2. Filter or analyze results
print(f"Found {len(result_df)} pathways")
print(result_df.head())

# 3. Save output
result_df.to_csv('disease_pathway_targets.tsv', sep='\t', index=False)

# 4. Get statistics
result_df['target_count'] = result_df['targets'].str.split(',').str.len()
print(f"Average targets per pathway: {result_df['target_count'].mean():.1f}")
```

## Notes

- The function automatically removes pathways with no targets
- Gene scores are coerced to numeric; NaN values are dropped
- The `propagated_edge` column contains all genes from the pathway definition (not just leading edge genes)

## License

See LICENSE file in the repository root.

