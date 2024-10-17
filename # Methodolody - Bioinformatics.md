# Methodolody - Bioinformatics

We also evaluated the expression of creatine-related genes using datasets and gene expression series from the Gene Expression Omnibus (GEO) database (Barret et al. 2012). We used each gene symbol as the main search term, applying the following criteria: (i) Organism, *Homo sapiens*; (ii) Differential expression, Up/down genes; (iii) DataSet keyword, Kidney and Renal tissues. The exclusion criteria were: undetected expression of the genes, experimental designs that did not compare disease/metabolic/chemical treatment conditions, or with mutant cell lines. The gene expression series were analyzed using the R packages GEO2R (https://www.ncbi.nlm.nih.gov/geo/geo2r/), GEOquery (Davis and Meltzer, 2007), limma (Smyth, 2005), and DESeq2 (Love et al. 2014) employing a case-control approach. We generated specific plots for evaluated genes using ggplot2 (Wickham et al. 2016) and Python programming language modules numpy (Harris et al. 2020) and matplotlib (Hunter, 2007). The scripts used are available on GitHub.


Barrett T, Wilhite SE, Ledoux P, Evangelista C, Kim IF, Tomashevsky M, Marshall KA, Phillippy KH, Sherman PM, Holko M, Yefanov A, Lee H, Zhang N, Robertson CL, Serova N, Davis S, Soboleva A.
NCBI GEO: archive for functional genomics data sets--update.
Nucleic Acids Res. 2013 Jan;41(Database issue):D991-5.

Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2. (Publisher link).

J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.

Smyth, G. K. (2005). Limma: linear models for microarray data. In: Bioinformatics and Computational Biology Solutions using R and Bioconductor, R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer, New York, pages 397-420.

Sean Davis and Paul S. Meltzer (2007). GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics 23(14): 1846-1847.

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.



Jared L. Johnson, Tomer M. Yaron, Emily M. Huntsman, Alexander Kerelsky, Junho Song, Amit Regev, Ting-Yu Lin, Katarina Liberatore, Daniel M. Cizin, Benjamin M. Cohen, Neil Vasan, Yilun Ma, Konstantin Krismer, Jaylissa Torres Robles, Bert van de Kooij, Anne E. van Vlimmeren, Nicole Andrée-Busch, Norbert F. Käufer, Maxim V. Dorovkov, Alexey G. Ryazanov, Yuichiro Takagi, Edward R. Kastenhuber, Marcus D. Goncalves, Benjamin D. Hopkins, Olivier Elemento, Dylan J. Taatjes, Alexandre Maucuer, Akio Yamashita, Alexei Degterev, Mohamed Uduman, Jingyi Lu, Sean D. Landry, Bin Zhang, Ian Cossentino, Rune Linding, John Blenis, Benjamin E. Turk, Michael B. Yaffe, and Lewis C. Cantley. An atlas of substrate specificities for the human serine/threonine kinome. Nature 613, 759–766 (2023). DOI: 10.1038/s41586-022-05575-3

Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, Ma’ayan A. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications 9. Article number: 1366 (2018), doi:10.1038/s41467-018-03751-6


