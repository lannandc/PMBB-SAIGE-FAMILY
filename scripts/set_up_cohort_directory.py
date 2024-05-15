import pandas as pd
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    parser.add_argument('-d', '--data', required=True, help='.csv Phenotype and covariate file')
    parser.add_argument('-c', '--cohort', required=True, help='Cohort to set up')
    parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments')
    parser.add_argument('-i', '--id', required=True, help='Column with sample IDs')
    parser.add_argument('--step1Fam', required=True)
    parser.add_argument('--exomeFam', required=True)

    return parser

args = make_arg_parser().parse_args()

id_col = args.id
cohort = args.cohort
step1_fam = args.step1Fam
exome_fam = args.exomeFam

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
samples = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

print(data)
print(samples)

step1_fam = pd.read_table(step1_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype=str)
exome_fam = pd.read_table(exome_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype=str)

#getting shared index labels across all 4 datasets
keep_samples = data.index.intersection(samples.index).intersection(step1_fam.index).intersection(exome_fam.index)

data, samples = data.loc[keep_samples], samples.loc[keep_samples]
cohort_samples = samples.index[samples[cohort] == 1]

open('sample_list.txt', 'w+').write('\n'.join(cohort_samples))

data = data.loc[cohort_samples]
data.index.name = id_col
data.to_csv('saige_pheno_covars.txt', sep='\t')
