# LUSC Bioinformatics Analysis

This project performs bioinformatics analysis of Lung Squamous Cell Carcinoma (LUSC).

## Data Source
- TCGA-LUSC

## Analysis Workflow
1. Data preprocessing
2. Differential expression analysis
3. GO enrichment analysis
4. KEGG pathway analysis

## conclusion
GO enrichment analysis indicated that differentially expressed genes (DEGs) were mainly involved in processes such as epidermis development, epithelial cell differentiation, leukocyte migration, chemotaxis, extracellular matrix organization, and DNA replication.

KEGG pathway analysis revealed significant enrichment in immune-related pathways, including systemic lupus erythematosus, neutrophil extracellular trap formation, complement and coagulation cascades, as well as extracellular matrix (ECM)–receptor interaction and cell cycle pathways.

These findings suggest that the DEGs are closely associated with squamous epithelial differentiation, immune and inflammatory responses, ECM remodeling, and enhanced cell proliferation in lung squamous cell carcinoma.

## Software
- R 4.5.2

packages:
- tidyverse, rjson, ggVolcano, ggplot2, stringr, pheatmap, limma, R.utils, GOplot, stringi, clusterProfiler, org.Hs.eg.db, enrichplot
