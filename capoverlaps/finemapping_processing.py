
# from liftover import ChainFile
# # Liftover chain
# chain37to38 = 'genetic_department/liftover/hg19ToHg38.over.chain.gz'
# # Liftover input data
# converter = ChainFile(chain37to38, one_based=True)
# finemapping['chromosome'] = finemapping['chromosome'].astype(int)
# finemapping[['chr_38', 'pos_38']] = finemapping.apply(lambda x: pd.Series(converter.convert_coordinate(x['chromosome'], x['position'])[0][:2]) if converter.convert_coordinate(x['chromosome'], x['position']) else pd.Series([None, None]), axis=1)
# finemapping = finemapping.dropna(subset=['chr_38', 'pos_38'])
# finemapping


import pandas as pd

# Paths to your BED files
path_to_file1 = 'comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_hg38_clustered_sorted_final.bed'
path_to_file2 = 'comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_Suzuki_finemapping.bed'

# Read the BED files without naming columns
df1 = pd.read_csv(path_to_file1, sep='\t', header=None)
df2 = pd.read_csv(path_to_file2, sep='\t', header=None)

# Merge the DataFrames on the column with index 3 (assuming 'rsid' is the fourth column)
merged_df = pd.merge(df1, df2, left_on=3, right_on=4, suffixes=('_file1', '_file2'))

merged_df.head()

# Save the merged data to a new BED file
#merged_df.to_csv('merged_data.bed', sep='\t', index=False, header=False)
# Now, merged_df contains rows from both dataframes matched by 'rsid'

#
#gene_count_trait = results.groupby(['trait'])['gene_name'].nunique().reset_index(name='gene_count')
