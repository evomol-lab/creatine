import matplotlib.pyplot as plt

# Plot a stacked bar plot
plt.figure(figsize=(10, 6))

# Define positions and mutation classes
positions = data['Position']
mutation_classes = ['ambiguous', 'benign', 'pathogenic']

# Create the stacked bar plot
plt.bar(positions, data['ambiguous'], label='Ambiguous', color='skyblue')
plt.bar(positions, data['benign'], bottom=data['ambiguous'], label='Benign', color='lightgreen')
plt.bar(positions, data['pathogenic'], bottom=data['ambiguous']+data['benign'], label='Pathogenic', color='salmon')

# Add labels and title
plt.xlabel('Position')
plt.ylabel('Number of Mutations')
plt.title('Number of Each Mutation Class per Position')
plt.legend()

# Show plot
plt.tight_layout()
plt.show()
