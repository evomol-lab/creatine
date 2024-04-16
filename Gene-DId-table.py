import pandas as pd
import re

# Load the TSV file
df = pd.read_csv('results3.tsv', sep='\t')

# Function to convert the 'DiseaseID' column from string to list
def convert_to_list(s):
    # Remove the square brackets and split the string into a list
    return re.findall(r"'(.*?)'", s)

# Convert the 'DiseaseID' column from string to list
df['DiseaseID'] = df['DiseaseID'].apply(convert_to_list)

# Explode the 'DiseaseID' column to have a separate row for each 'DiseaseID'
df_exploded = df.explode('DiseaseID')

# Save the results to a new TSV file
df_exploded.to_csv('results4.tsv', sep='\t', index=False)
