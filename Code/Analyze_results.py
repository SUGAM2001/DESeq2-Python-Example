import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# Load the DESeq2 results
res_df = pd.read_csv('C:/Users/sugam patel/Downloads/deseq2_results.csv', index_col=0)


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
