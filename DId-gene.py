import pandas as pd

# Load the TSV file
df = pd.read_csv('CTD_gene_diseases_inferred_Kidney-CTD.tsv', sep='\t')

# Group by 'GeneSymbol' and get unique 'DiseaseID' for each 'GeneSymbol'
grouped = df.groupby('GeneSymbol')['DiseaseID'].unique()

# Convert the Series to a DataFrame and reset the index
result = pd.DataFrame(grouped).reset_index()

# Save the results to a new TSV file
result.to_csv('results3.tsv', sep='\t', index=False)

