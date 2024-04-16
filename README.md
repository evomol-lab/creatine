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

The first suggested protein list is ([Stitch-1st](Stitch-1st)):
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
Our second step was to perform searches using the [Comparative Toxicogenomics Database](https://ctdbase.org/). As a result, we obtained a list with 44 entries [CTD_gene_results_20240409072043.tsv](CTD_gene_results_20240409072043.tsv). After filtering the repetitive entries, we obtained the following gene list:

```
CKB
CKBE    Locus
CKBP1   Pseudogene
CKM
CKMT1A
CKMT1B
CKMT2
SLC16A12
SLC6A10P    Pseudogene
SLC6A10PB   Pseudogene
SLC6A8
```
Interestingly, three pseudogenes and one *locus* were retrieved. The informations from each gene symbol were obtained from [GeneCards](https://www.genecards.org/) database. The validated symbols were the following ones [CTD-genes-list-Edited-1st.tsv](CTD-genes-list-Edited-1st.tsv):

```
CKB
CKM
CKMT1A
CKMT1B
CKMT2
SLC16A12
SLC6A8
```
After that, we used the batch query tool from CTD, using NCBI gene symbols from the two above lists as a input. Then we got the following results:

- For [Stitch-1st](Stitch-1st) list:
    - Curated Gene/diseases associations: [CTD_gene_diseases_curated_1713265713318.tsv](CTD_gene_diseases_curated_1713265713318.tsv).
    - Inferred Gene/diseases associations: [CTD_gene_diseases_inferred_1713266126277.tsv](CTD_gene_diseases_inferred_1713266126277.tsv).
    - Curated Gene/pathways associations: [CTD_gene_pathways_curated_1713266155612.tsv](CTD_gene_pathways_curated_1713266155612.tsv).
    - Chemical/Genes interactions at the expression level: [CTD_gene_cgixns_1713266196128.tsv](CTD_gene_cgixns_1713266196128.tsv).
    - Curated Gene/Chemicals interactions: [CTD_gene_chems_curated_1713267587656.tsv](CTD_gene_chems_curated_1713267587656.tsv).
- For [CTD-genes-list-Edited-1st.tsv](CTD-genes-list-Edited-1st.tsv):
    - Curated Gene/diseases associations: [CTD_gene_diseases_curated_CTD-list.tsv](CTD_gene_diseases_curated_CTD-list.tsv).
    - Inferred Gene/diseases associations: [CTD_gene_diseases_inferred_CTD-list.tsv](CTD_gene_diseases_inferred_CTD-list.tsv).
    - Curated Gene/pathways associations: [CTD_gene_pathways_curated_CTD-list.tsv](CTD_gene_pathways_curated_CTD-list.tsv).
    - Chemical/Genes interactions at the expression level: [CTD_gene_cgixns_CTD-list.tsv](CTD_gene_cgixns_CTD-list.tsv).
    - Curated Gene/Chemicals interactions: [CTD_gene_chems_curated_CTD-list.tsv](CTD_gene_chems_curated_CTD-list.tsv).

#### Filtering CTD's results
We started filtering the inferred gene/diseases associations, since the curated ones are enriched to genetic diseases (inborn). For the Stitch-1st and CTD lists, we obtained a file with 153952 and 50868 rows, respectively. Then, we filtered the results using the presence of the terms `kidney` and `renal`, obtaining 2596 and 818 rows, for Stitch and CTD gene lists. The higher number of entries in the Stitch list is probably associated with the initial filtering process for expression in kidney tissues. The following table describes the number of entries per gene, in each list:

| Stitch |                   | CTD   |                   |
|:------:|:-----------------:|:-----:|:-----------------:|
| **Gene**   | **Number of entries** | **Gene**  | **Number of entries** |
| SLC2A4     |  283 | CKB      | 214 |
| IGF1       |  431 | CKMT1A   | 39  |
| GATM       |  172 | CKMT1B   | 37  |
| SLC6A8     |  157 | CKMT2    | 90  |
| AKT1       |  691 | SLC16A12 | 114 |
| AKT2       |  192 | SLC6A8   | 157 |
| AKT3       |  111 |       |                   |
| PSMB5      |  119 |       |                   |
| PSMD3      |  84  |       |                   |
| SGK1       |  268 |       |                   |
| PRPS1      |  109 |       |                   |

From the filtered files, we obtained the unique DiseaseIDs, and unique DiseaseIDs/gene.

![Stitch Genes with associated DiseaseIDs related to Kidney/Renal diseases](Gene-DId-Stich.png)
Stitch Genes with associated DiseaseIDs related to Kidney/Renal diseases
![CTD Genes with associated DiseaseIDs related to Kidney/Renal diseases](Gene-DId-CTD.png)
CTD Genes with associated DiseaseIDs related to Kidney/Renal diseases
![Presence/Absence of DiseaseID per gene from the Stitch-list](Facetted-DId.png)
Presence/Absence of DiseaseID per gene from the Stitch-list
![Presence/Absence of DiseaseID per gene from the CTD-list](Facetted-CTD.png)
Presence/Absence of DiseaseID per gene from the CTD-list

Here are the number of DiseaseIDs and counts, common to the two genes lists:
```
DiseaseID      
MESH:D007676       77
MESH:D007674       77
MESH:D007680       77
MESH:D007683       77
MESH:D058186       77
MESH:D007681       66
MESH:D007673       66
MESH:D007669       60
MESH:D052177       40
MESH:C537152       32
MESH:D007690       32
MESH:D021782       32
MESH:D012080       10
MESH:C538445        9
MESH:D000092702     8
```
The list from the Stitch database has these unique DiseaseID counts:
```
DiseaseID
MESH:D016891    3
MESH:C531755    1
```

### Resultado enrichr-KG
The subnetwork shows the following associations: 
- **From Gene Ontology:** PSMB5, and PSMD3 belong to the biological process regulation of cellular amino acid metabolic process (GO:0006521). PSMB5, and PSMD3 belong to the biological process regulation of cellular ketone metabolic process (GO:0010565). SLC6A8, and GATM belong to the biological process creatine metabolic process (GO:0006600). PSMB5, and PSMD3 belong to the biological process regulation of cellular amine metabolic process (GO:0033238). PSMB5, and PSMD3 belong to the biological process negative regulation of cell cycle G2/M phase transition (GO:1902750).
- **From KEGG**: The gene products PSMB5, and PSMD3 are members of the KEGG pathway Proteasome. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Huntington disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Prion disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Parkinson disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Spinocerebellar ataxia.
- **From Jensen lab:** The disease Arts syndrome is associated with the gene PRPS1. The disease Intellectual disability is associated with the genes SLC6A8, and PRPS1. The disease Gyrate atrophy is associated with the gene GATM. The disease Purine nucleoside phosphorylase deficiency is associated with the gene PRPS1. The disease AGAT deficiency is associated with the gene SLC6A8. 
- **From DisGeNET:** The disease Mammary Carcinoma, Animal is associated with the following genes: SGK1, and GATM. The disease Drug Resistant Epilepsy is associated with the following genes: SLC6A8, and SGK1. The disease Congenital pes cavus is associated with the following genes: SLC6A8, and PRPS1. The disease Neonatal Hypotonia is associated with the following genes: SLC6A8, and PRPS1. The disease Creatine deficiency, X-linked is associated with the following genes: SLC6A8, and GATM.