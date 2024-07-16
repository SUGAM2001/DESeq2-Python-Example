#import important library

import pandas as pd
import numpy as np
from rpy2.robjects import pandas2ri,r,Formula
from rpy2.robjects.packages import import
from rpy2.robjects import conversion, default_converter
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# Import DESeq2
deseq2 = importr('DESeq2')
pandas2ri.activate()

# load the dataset and convert it into R dataframe
file_data = "C:/Users/sugam patel/Downloads/filtered_gene_expression_counts.csv"
exp_count = pd.read_csv(file_data)
print(exp_count)  

#Convert it to R DataFrame
r_df = pandas2ri.py2rpy(exp_count)
print(r_df)

# Ensure the countData has correct column names
exp_count.columns = ['count.C1', 'count.C2', 'count.C3', 'count.T1', 'count.T2', 'count.T3']

# Convert to R DataFrame
with conversion.localconverter(default_converter + pandas2ri.converter):
    r_df = pandas2ri.py2rpy(exp_count)

# # Import DESeq2
# deseq2 = importr('DESeq2')

# Create the design formula
design = Formula('~ condition')

# Create the colData DataFrame in Python
col_data = pd.DataFrame({
    'condition': ['control', 'control', 'control', 'treatment', 'treatment', 'treatment']
}, index=['count.C1', 'count.C2', 'count.C3', 'count.T1', 'count.T2', 'count.T3'])

# Convert colData to R DataFrame and ensure 'condition' is a factor
with conversion.localconverter(default_converter + pandas2ri.converter):
    r_col_data = pandas2ri.py2rpy(col_data)
r_col_data = r['within'](r_col_data, 'condition <- as.factor(condition)')

# Ensure row names of r_col_data match column names of r_df
r['rownames'](r_col_data, r['colnames'](r_df))

# Create the DESeqDataSet
dds = deseq2.DESeqDataSetFromMatrix(countData=r_df, colData=r_col_data, design=design)

# Run the DESeq function
dds = deseq2.DESeq(dds)

# Get the results and convert to a data frame in R
res = deseq2.results(dds)
r_res_df = r['as.data.frame'](res)

# Convert the R data frame to a Pandas data frame
with conversion.localconverter(default_converter + pandas2ri.converter):
    res_df = pandas2ri.rpy2py(r_res_df)

# # Save the results to a CSV file
 res_df.to_csv('C:/Users/sugam patel/Downloads/deseq2_results.csv')

# Set the criteria for the significant gene
significant_genes = res_df[res_df['padj'] < 0.05]

# Volcano plot
plt.figure(figsize=(10, 6))
sns.scatterplot(x='log2FoldChange', y=-np.log10(res_df['padj']), data=res_df, hue=(res_df['padj'] < 0.05))
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-Value')
plt.title('Volcano Plot of Differential Expression')
plt.legend(title='Significant')
plt.show()

# MA Plot
plt.figure(figsize=(10, 6))
sns.scatterplot(x=res_df['baseMean'], y=res_df['log2FoldChange'], hue=(res_df['padj'] < 0.05))
plt.axhline(0, color='grey', lw=1)
plt.xlabel('Mean of Normalized Counts')
plt.ylabel('Log2 Fold Change')
plt.title('MA Plot of Differential Expression')
plt.legend(title='Significant')
plt.xscale('log')
plt.show()

# Histogram of p-values
plt.figure(figsize=(10, 6))
sns.histplot(res_df['pvalue'].dropna(), bins=50, kde=True)
plt.xlabel('P-value')
plt.ylabel('Frequency')
plt.title('Histogram of P-values')
plt.show()

# Ensure the index of res_df is the same as the index of exp_count
# This step assumes that your gene IDs or names are in the index
exp_count.index = exp_count.index.astype(str)
res_df.index = res_df.index.astype(str)

# Extract counts of significant genes
significant_genes = res_df[res_df['padj'] < 0.05]
significant_gene_names = significant_genes.index

# Filter the original count data for significant genes
# Ensure the significant genes are present in the exp_count DataFrame
significant_counts = exp_count.loc[exp_count.index.intersection(significant_gene_names)]

# Standardize the counts
scaler = StandardScaler()
scaled_counts = scaler.fit_transform(significant_counts)

# Create the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(scaled_counts, cmap='viridis', xticklabels=exp_count.columns, yticklabels=significant_counts.index)
plt.title('Heatmap of Significant Genes')
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.show()






