import pandas as pd

# Load the TSV files
df1 = pd.read_csv('results2.tsv', sep='\t')
df2 = pd.read_csv('results4.tsv', sep='\t')

# Get the DiseaseIDs that are present in both files
common_disease_ids = pd.merge(df1['DiseaseID'], df2['DiseaseID'], how='inner')

# Count the number of each DiseaseID
disease_id_counts = common_disease_ids.value_counts()

# Sort the DiseaseIDs by their counts
sorted_disease_ids = disease_id_counts.sort_values(ascending=False)

# Print the sorted DiseaseIDs and their counts
print(sorted_disease_ids)

# Get the DiseaseIDs that are present only in one of the files
distinct_disease_ids_df1 = df1[~df1['DiseaseID'].isin(df2['DiseaseID'])]['DiseaseID']
distinct_disease_ids_df2 = df2[~df2['DiseaseID'].isin(df1['DiseaseID'])]['DiseaseID']

# Count the number of each DiseaseID
disease_id_counts_df1 = distinct_disease_ids_df1.value_counts()
disease_id_counts_df2 = distinct_disease_ids_df2.value_counts()

# Sort the DiseaseIDs by their counts
sorted_disease_ids_df1 = disease_id_counts_df1.sort_values(ascending=False)
sorted_disease_ids_df2 = disease_id_counts_df2.sort_values(ascending=False)

# Print the sorted DiseaseIDs and their counts
print("Distinct DiseaseIDs in file1:")
print(sorted_disease_ids_df1)
print("\nDistinct DiseaseIDs in file2:")
print(sorted_disease_ids_df2)

