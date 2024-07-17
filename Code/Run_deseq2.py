#import important library

import pandas as pd
import numpy as np
from rpy2.robjects import pandas2ri,r,Formula
from rpy2.robjects.packages import import
from rpy2.robjects import conversion, default_converter

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

