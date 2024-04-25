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
MESH:D007676       77   Kidney Failure, Chronic
MESH:D007674       77   Kidney Diseases
MESH:D007680       77   Kidney Neoplasms
MESH:D007683       77   Kidney Tubular Necrosis, Acute
MESH:D058186       77   Acute Kidney Injury
MESH:D007681       66   Kidney Papillary Necrosis
MESH:D007673       66   Kidney Cortex Necrosis
MESH:D007669       60   Kidney Calculi
MESH:D052177       40   Kidney Diseases, Cystic
MESH:C537152       32   Hypomagnesemia 2, renal [Supplementary Concept]
MESH:D007690       32   Polycystic Kidney Diseases
MESH:D021782       32   Multicystic Dysplastic Kidney
MESH:D012080       10   Chronic Kidney Disease-Mineral and Bone Disorder
MESH:C538445        9   Clear-cell metastatic renal cell carcinoma [Supplementary Concept]
MESH:D000092702     8   Chronic Kidney Diseases of Uncertain Etiology
```
The list from the Stitch database has these unique DiseaseID counts:
```
DiseaseID
MESH:D016891    3   Polycystic Kidney, Autosomal Dominant
MESH:C531755    1   Kidney disorder involving deposition of calcium and oxalate or phosphate in the renal tubules [Supplementary Concept]
```

## Gene enrichment analyses

We merged the two above lists and use them as input for a gene enrichment analysis (GEA), using enrichR and enrichR-KG.

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
CKB
CKM
CKMT1A
CKMT1B
CKMT2
SLC16A12
```
The first filter used was CellTypes enrichment analysis. As a selection criterion, we chose the databases where the words "kidney" or "renal" were present with associated terms (P-value<0.05), from the all tested databases. Using the `.tsv` tables as input, we executed the following command:

```
awk -F'\t' '$1 ~ /Kidney/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' listName
```
We also searched against Kidney cell lines:
- **ARCHS4 cell-lines database**: 293F, FLPIN TREX 293, and HEK293.


Only three databases retrieved results within the selection criteria. The table below describes the terms and their associated values:

| Database | CellTypes | P-value | adjusted P-value | OddsRatio | Combined Score | Genes |
|:--------:|:---------:|:-------:|:----------------:|:---------:|:--------------:|:-----:|
| **CellMarker 2024** | Collecting Duct Intercalated Cell Kidney Mouse | 0.003593 | 0.07843 | 25.49 | 143.46 | CKMT1B, SGK1 |
| **HuBMAP ASCTplusB augmented 2022** | Kidney Interstitial fibroblast - Kidney | 8.808317003142788E-5 | 0.0013050178709848516 | 42.18246110325318 | 393.86730230343085 | GATM, SLC16A12, IGF1 |
| | Connecting Tubule Intercalated Cell Type A - Kidney | 0.0033964046824717384 | 0.022655611815195546 | 26.246864686468648 | 149.2144192699383 | CKMT2, SLC16A12 |
| | Kidney Outer Medulla Collecting Duct Intercalated cell - Kidney | 0.0034612740273215416 | 0.022655611815195546 | 25.988235294117647 | 147.2524218660231 | CKMT2, SLC16A12 |
| **MoTrPAC 2023** | T59-Kidney Consensus | 0.0027283289091587022 | 0.024323568734636137 | 12.568017057569296 | 74.20240192548566 | SLC2A4, IGF1, CKB|
| | T59-Kidney Female 4W Up | 0.004243143988996933 | 0.026702092577683224 | 312.171875 | 1705.2235013634931 | SLC2A4 |
| **ARCHS4 cell-lines database** | HEK293 | 0.04417020235184285 | 0.8857370659521143 | 3.0671199442119943 | 9.568509033514891 | PRPS1, CKMT1A, AKT3, CKMT1B, CKB |
| **ProteomicsDB_2020** | Renal SN-12C BTO:0004221 P0001751 | 0.04134985878495573 | 0.44706319361150326 | 6.841535776614311 | 21.794986607296053 | PRPS1, AKT1 |
| **ARCHS4 Tissues** | KIDNEY (BULK TISSUE) | 2.952357524640197E-4 | 0.025095038959441675 | 6.807240516079338 | 55.327456253044275 | PRPS1, SLC6A8, PSMB5, CKMT1A,PSMD3, CKMT1B, SLC2A4, CKB |

>From Copilot:
>
>**GATM (Glycine Amidinotransferase)**:
>*GATM plays a crucial role in the synthesis of guanidinoacetic acid (GAA), a metabolite associated with pancreatic ductal adenocarcinoma (PDAC) liver metastasis¹>[6]. High expression of GATM has been positively correlated with advanced N stage in PDAC¹[6]. Knockdown of GATM significantly reduced the intracellular level of >GAA, suppressed epithelial-mesenchymal transition (EMT), and inhibited PDAC liver metastasis¹[6]. These findings suggest that GATM-mediated de novo GAA synthesis >promotes PDAC metastasis¹[6].*
>
>**SLC16A12 (Solute Carrier Family 16 Member 12)**:
>*The SLC16A family members, including SLC16A12, play crucial roles in tumorigenesis and tumor progression²[1]. However, the exact role of SLC16A12 in human >pancreatic cancer remains unclear²[1]. The SLC16A family of monocarboxylate transporters (MCTs) is known to play an important role in cell metabolism, such as >aerobic glycolysis and pH homeostasis²[1]. Current research progress on the MCT family, including SLC16A12, is mainly based on the findings of MCT1-4²[1]. These >proteins exhibit broad substrate specificity and are also involved in the transport of other monocarboxylic metabolites²[1].*
>
>*Please note that while these genes and their proteins have been associated with pancreatic cancer, the exact mechanisms and their potential as therapeutic >targets are still under investigation. For more detailed information, please refer to the respective scientific research articles²[1]¹[6].*
>
>**References:**
[1]: https://www.nature.com/articles/s41598-020-64356-y.pdf.
[2]: https://www.mdpi.com/2072-6694/10/4/103. 
[3]: https://www.medrxiv.org/content/10.1101/2024.03.03.24303664v1.full.pdf. 
[4]: https://columbiasurgery.org/pancreas/genetics-pancreatic-cancer. 
[5]: https://link.springer.com/article/10.1007/s11605-022-05553-0. 
[6]: https://jeccr.biomedcentral.com/articles/10.1186/s13046-023-02698-x. 
[7]: https://www.nature.com/articles/s41392-023-01662-7.pdf. 
[8]: https://www.nebraskamed.com/cancer/pancreatic/the-role-that-genes-play-in-pancreatic-cancer. 
[9]: https://doi.org/10.3390/cancers10040103.





### Resultado enrichr-KG
The subnetwork shows the following associations: 
- **From Gene Ontology:** PSMB5, and PSMD3 belong to the biological process regulation of cellular amino acid metabolic process (GO:0006521). PSMB5, and PSMD3 belong to the biological process regulation of cellular ketone metabolic process (GO:0010565). SLC6A8, and GATM belong to the biological process creatine metabolic process (GO:0006600). PSMB5, and PSMD3 belong to the biological process regulation of cellular amine metabolic process (GO:0033238). PSMB5, and PSMD3 belong to the biological process negative regulation of cell cycle G2/M phase transition (GO:1902750).
- **From KEGG**: The gene products PSMB5, and PSMD3 are members of the KEGG pathway Proteasome. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Huntington disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Prion disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Parkinson disease. The gene products PSMB5, and PSMD3 are members of the KEGG pathway Spinocerebellar ataxia.
- **From Jensen lab:** The disease Arts syndrome is associated with the gene PRPS1. The disease Intellectual disability is associated with the genes SLC6A8, and PRPS1. The disease Gyrate atrophy is associated with the gene GATM. The disease Purine nucleoside phosphorylase deficiency is associated with the gene PRPS1. The disease AGAT deficiency is associated with the gene SLC6A8. 
- **From DisGeNET:** The disease Mammary Carcinoma, Animal is associated with the following genes: SGK1, and GATM. The disease Drug Resistant Epilepsy is associated with the following genes: SLC6A8, and SGK1. The disease Congenital pes cavus is associated with the following genes: SLC6A8, and PRPS1. The disease Neonatal Hypotonia is associated with the following genes: SLC6A8, and PRPS1. The disease Creatine deficiency, X-linked is associated with the following genes: SLC6A8, and GATM.