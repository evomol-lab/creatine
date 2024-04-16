# Load the necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Load the TSV file
df <- read_tsv('results.tsv')

# Count the distinct 'DiseaseID' for each 'GeneSymbol'
df_count <- df %>%
  group_by(GeneSymbol) %>%
  summarise(DiseaseCount = n_distinct(DiseaseID))

# Plot the counts
ggplot(df_count, aes(x = GeneSymbol, y = DiseaseCount)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Gene Symbol", y = "Distinct Disease ID Count", title = "Count of Distinct DiseaseID per GeneSymbol")

# Load the necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Load the TSV file
df2 <- read_tsv('results2.tsv')

# Count the number of each 'DiseaseID' for each 'GeneSymbol'
df_count <- df2 %>%
  group_by(GeneSymbol, DiseaseID) %>%
  summarise(Count = n())

# Create a facetted bar plot
ggplot(df_count, aes(x = DiseaseID, y = Count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ GeneSymbol, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Disease ID", y = "Count", title = "Count of Each DiseaseID for Each GeneSymbol")
