import pandas as pd

def normalize_positions(peptide_df, chr_df):
    # Merge the peptide info with chromosome info based on Chr column
    merged_df = pd.merge(peptide_df, chr_df, left_on='Chr', right_on='Chr')

    # Define a function to compute normalized position
    def compute_normalized_position(row):
        if row['Start'] < row['cen']:
            # Peptide is between tel1 and cen
            return (row['Start'] - row['tel1']) / (row['cen'] - row['tel1'])
        else:
            # Peptide is between cen and tel2
            return (row['Start'] - row['cen']) / (row['tel2'] - row['cen'])

    # Apply the function to each row
    merged_df['Normalized_Position'] = merged_df.apply(compute_normalized_position, axis=1)
    return merged_df[['Peptide', 'Chr', 'Start', 'Strand', 'Normalized_Position']]

# Load the chromosome and peptide data from previously specified paths
chr_info_df = pd.read_excel('chrinfo.xlsx')
peptide_info_df = pd.read_excel('cp-cp-sp.xlsx')

# Apply the normalization function
normalized_peptides_df = normalize_positions(peptide_info_df, chr_info_df)

# Bin the normalized positions into intervals of 0.01
bins = pd.interval_range(start=0, end=1, freq=0.01, closed='left')
normalized_peptides_df['Bin'] = pd.cut(normalized_peptides_df['Normalized_Position'], bins=bins)

# Calculate the count of peptides in each bin
peptide_distribution = normalized_peptides_df['Bin'].value_counts().sort_index()

# Save the normalized positions to a new Excel file
output_path = 'normalized_peptides_sp.xlsx'
normalized_peptides_df.to_excel(output_path, index=False)

# Save the peptide distribution to another Excel file
output_distribution_path = 'peptide_distribution_sp.xlsx'
peptide_distribution.to_excel(output_distribution_path)

output_path, output_distribution_path
