# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")

install.packages("ggplot2")

# Import necessary libraries
library(GEOquery)
library(ggplot2)

# Function to read GEO expression soft file
read_geo_file <- function(file_path) {
    # Read the file
    gset <- getGEO(filename=file_path, destdir=".")
    
    # Check the class of the returned object
    print(class(gset))
    
    # Extract expression data based on the class of the object
    if (class(gset) == "list") {
        data <- exprs(gset[[1]])
    } else if (class(gset) == "ExpressionSet") {
        data <- exprs(gset)
    } else if (class(gset) == "GDS") {
        data <- Table(gset)
    } else {
        stop("Unexpected class of object returned by getGEO: ", class(gset))
    }
    
    return(data)
}

# Function to plot gene expression
plot_gene_expression <- function(data, gene_name) {
    # Get the expression values for the given gene
    gene_data <- data[gene_name, ]
    
    # Create a data frame for plotting
    plot_data <- data.frame(Sample = colnames(data), Expression = gene_data)
    
    # Create the plot
    plot <- ggplot(plot_data, aes(x = Sample, y = Expression)) +
        geom_line() +
        labs(title = paste('Expression values for', gene_name), x = 'Sample', y = 'Expression value')
    
    # Print the plot
    print(plot)
}

# Use the functions
file_path <- 'GDS3712.soft'  # Replace with your file path
data <- read_geo_file(file_path)

gene_name <- 'GATM'  # Replace with your gene name
plot_gene_expression(data, gene_name)
