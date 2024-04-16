import pandas as pd

# Load the TSV file
df = pd.read_csv('CTD_gene_diseases_inferred_Kidney-CTD.tsv', sep='\t')

# Get the unique values of the 'DiseaseID' column
unique_values = df['DiseaseID'].unique()

# Save the unique values to a new text file
with open('unique_values.txt', 'w') as f:
    for value in unique_values:
        f.write(str(value) + '\n')
