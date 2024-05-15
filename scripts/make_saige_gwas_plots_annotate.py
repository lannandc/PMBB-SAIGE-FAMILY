from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional list argument for phenotype
    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')

    # Add a non-optional list argument for cohort
    parser.add_argument('-c', '--cohort', required=True, help='cohort')
    parser.add_argument('-s', '--sumstats', required=True, help='Path to summary statistics file')
    parser.add_argument('-t', '--phenoTable', required=True)
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None,
                        help='Path to output directory. Default: current working directory')
    parser.add_argument('-a', '--annot', required=False, default=None)

    return parser


args = make_arg_parser().parse_args()
cohort = args.cohort
pheno = args.phenotype
annot_file = args.annot

sumstats = args.sumstats
pheno_table = args.phenoTable
output_dir = args.outDir
output_dir += '/' if output_dir[-1] != '/' else ''

pheno_df = pd.read_csv(pheno_table)
pheno_df = pheno_df[pheno_df['PHENO'] == pheno]
pheno_df = pheno_df[pheno_df['COHORT'] == cohort]


trait_type = 'bin' if pheno_df['Cases'].count() != 0 else 'quant'

# specify outdir if given
if output_dir:
    output_manhattan = f'{output_dir}/{cohort}.{pheno}.manhattan_vertical.png'
    output_qq = f'{output_dir}/{cohort}.{pheno}.qq.png'
else:
    output_manhattan = f'{cohort}.{pheno}.manhattan_vertical.png'
    output_qq = f'{cohort}.{pheno}.qq.png'

plot_title = f'SAIGE GWAS Annotated {cohort}: {pheno.replace("_", " ")}'

if trait_type == 'bin':
    plot_title += f'\nCases = {pheno_df.iloc[0].Cases:,.0f}, Controls = {pheno_df.iloc[0].Controls:,.0f}'
else:
    plot_title += f'\nN = {pheno_df.iloc[0].N:,.0f}'

mp = ManhattanPlot(sumstats, test_rows=None, title=plot_title)
mp.load_data()

mp.df['chromosome_noCHR'] = mp.df['chromosome'].str.replace('chr', '').astype(int)
mp.clean_data(col_map={'chromosome_noCHR': '#CHROM', 'base_pair_location': 'POS', 'p_value': 'P', 'variant_id': 'ID'})

mp.get_thinned_data()
print(mp.thinned)
print(len(mp.thinned))

# edge case protection
# if ~np.any(mp.thinned['P'] < 5E-8):
#     p_thresh = np.quantile(mp.thinned['P'], 10 / len(mp.thinned))
# else:
# 		p_thresh = 5E-8

annot_thresh = 1E-5 if np.any(mp.thinned['P'].min() < 1E-5) else np.nanquantile(mp.thinned['P'], 10 / len(mp.thinned))


mp.update_plotting_parameters(vertical=True, sig=annot_thresh if not np.any(mp.thinned['P'] < 5E-8) else 5E-8, 
                              sug=annot_thresh, annot_thresh=annot_thresh, merge_genes=True)

mp.full_plot(save=output_manhattan, rep_boost=False, extra_cols={
             'beta': 'beta', 'RSID': 'RSID'}, number_cols=['beta', 'RSID'], keep_chr_pos=False)
plt.clf()

print(f"Saved Manhattan plot to: {output_manhattan}")

mp.qq_plot(save=output_qq)

print(f"Saved qq plot to: {output_qq}")





# mp.df.rename(columns=columns_map_inv, inplace=True)
# get rows where POS==UR, select POS col; 
# select rows from mp where POS=UR, select MarkerID col; split by : and grab first;
# use split string to find the index in genefile, then get "START" and CHR col. 
# assign numpy array to POS column in selected rows of mp
# mp.df.loc[mp.df['POS'] == 'UR', 'POS'] = gene_file.loc[mp.df.loc[mp.df['POS'] == 'UR', 'MarkerID'].str.split(':', expand=True)[0], 'START'].values
# mp.df.loc[mp.df['CHR'] == 'UR', 'CHR'] = gene_file.loc[mp.df.loc[mp.df['CHR'] == 'UR', 'MarkerID'].str.split(':', expand=True)[0], 'CHR'].values


"""
# remap columns from column map file
infile = 'colnames.txt'
# Initialize an empty dictionary
columns_map = {}
with open(infile, 'r') as file:
    for line in file:
        # Split the line on '=' and strip whitespace and quotes
        key, value = line.split('=')
        key = key.strip()
        value = value.strip().strip("'")

        # Add to the dictionary
        columns_map[key] = value
# reverse code the columns map
columns_map_inv = {v: k for k, v in columns_map.items()}
"""