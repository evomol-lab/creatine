# Creatine
Files for a bioinformatic analyses on creatine effects on kidney tissues.

# Creatine review
## Running Title
*Creatine effects on kidney tissues and renal function: new insights to old questions.*
## First Steps
### Creatine-related proteins on Stitch Database
Our first step was to use the term "*creatine*" in a search using the [Stitch Database](http://stitch.embl.de/) ([PUBMED](https://pubmed.ncbi.nlm.nih.gov/26590256/)), which list the known and predicted interactions between chemicals and proteins. The search was limited to the organism *Homo sapiens*. We found the proteins listed in the file [stitch_interactions-1st.tsv](stitch_interactions-1st.tsv). Then, using the advanced settings, we restricted the search for integration of RNA-Seq expression data from kidney tissues, using the [Human Protein Atlas](https://www.proteinatlas.org/) database. The file [stitch_interactions-2nd-ProteinAtlas.tsv](stitch_interactions-2nd-ProteinAtlas.tsv) contain the updated list, where the proteins AKT1, AKT2 and AKT3 appeared as a tissue-specific interactions. We also used the restriction from [Tissues Database](https://tissues.jensenlab.org/Search), using kidney as a term [stitch_interactions-3rd.tsv](stitch_interactions-3rd.tsv).
#### First observations
- The genes GATM and SLC6A8 were present in all three above searches.
- The genes from AKT family were present using expression data from the [Human Protein Atlas](https://www.proteinatlas.org/).
- The genes PSMB5, PSMD3, SGK1 and PRPS1 appeared when using expression data from the [Tissues Database](https://tissues.jensenlab.org/Search).

The first suggested protein list is:
```
SLC2A4
IGF1
GATM
SLC6A8
AKT1
AKT2
AKT3
PSMB5
PSMD3
SGK1
PRPS1
```
Since our aim is to review gene expression and protein alterations induced by creatine that are specific to renal tissues, we excluded the genes that are expected to interact with creatine, such as CK and CKMT. Interestingly, these genes were not found in the 3rd list (from Tissues database).

### Using the Comparative Toxicogenomics Database
Our second step was to perform searches using the [Comparative Toxicogenomics Database](https://ctdbase.org/).



### Resultado enrichr-KG
The subnetwork shows the following associations: 
- **From Gene Ontology:** PSMB5, and PSMD3 belong to the biological process regulation of cellular amino acid metabolic process (GO:0006521). PSMB5, and PSMD3 belong to the biological process regulation of cellular ketone metabolic process (GO:0010565). SLC6A8, and GATM belong to the biological process creatine metabolic process (GO:0006600). PSMB5, and PSMD3 belong to the biological process regulation of cellular amine metabolic process (GO:0033238). PSMB5, and PSMD3 belong to the biological process negative regulation of cell cycle G2/M phase transition (GO:1902750).
- **From KEGG**: The gene products PSMB5, and PSMD3 are members of the KEGG pathway Proteasome. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Huntington disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Prion disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Parkinson disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Spinocerebellar ataxia.
- **From Jensen lab:** The disease Arts syndrome is associated with the gene PRPS1. The disease Intellectual disability is associated with the genes SLC6A8, and PRPS1. The disease Gyrate atrophy is associated with the gene GATM. The disease Purine nucleoside phosphorylase deficiency is associated with the gene PRPS1. The disease AGAT deficiency is associated with the gene SLC6A8. 
- **From DisGeNET:** The disease Mammary Carcinoma, Animal is associated with the following genes: SGK1, and GATM. The disease Drug Resistant Epilepsy is associated with the following genes: SLC6A8, and SGK1. The disease Congenital pes cavus is associated with the following genes: SLC6A8, and PRPS1. The disease Neonatal Hypotonia is associated with the following genes: SLC6A8, and PRPS1. The disease Creatine deficiency, X-linked is associated with the following genes: SLC6A8, and GATM.