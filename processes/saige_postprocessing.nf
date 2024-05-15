
process merge_and_filter_saige_gene_regions_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr_list), path(chr_inputs)

        path merge_regions_script
    output:
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_regions.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_regions.filtered.saige.csv")
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_cauchy_tests.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_cauchy_tests.filtered.saige.csv")
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          -p ${pheno} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Sumstats/

        touch ${pheno}.gene_regions.saige.gz
        touch ${pheno}.gene_regions.filtered.saige.csv
        """
}

process merge_and_filter_saige_gene_singles_phewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"

    input:
        // variables
        tuple val(cohort_dir), val(chr), path(chr_inputs)
        path merge_singles_script
    output:
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_singles.saige.gz")
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_singles.filtered.saige.csv")
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_ultrarare_tests.saige.gz")
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_ultrarare_tests.filtered.saige.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_singles_script} \
          --chr ${chr} \
          --phewas \
          --rareVars \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        touch ${pheno}.gene_phewas_singles.saige.gz
        touch ${pheno}.gene_phewas_singles.filtered.saige.csv
        """
}

process merge_and_filter_saige_gene_regions_phewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"

    input:
        // variables
        tuple val(cohort_dir), val(chr), path(pheno_inputs)

        path merge_regions_script
    output:
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_regions.saige.gz")
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_regions.filtered.saige.csv")
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_cauchy_tests.saige.gz")
        tuple val(cohort_dir), val(chr), path("chr${chr}.gene_phewas_cauchy_tests.filtered.saige.csv")
    shell:
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          --phewas \
          --regions \
          --chr ${chr} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${pheno_inputs.join(' ')}
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Sumstats/

        touch chr${chr}.gene_phewas_regions.saige.gz
        touch chr${chr}.gene_phewas_regions.filtered.saige.csv
        touch chr${chr}.gene_phewas_cauchy_tests.saige.gz
        touch chr${chr}.gene_phewas_cauchy_tests.filtered.saige.csv
        """
}

process merge_and_filter_saige_gene_singles_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr_list), path(chr_inputs)
        path merge_singles_script
    output:
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_singles.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_singles.filtered.saige.csv")
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_ultrarare_tests.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${pheno}.exwas_ultrarare_tests.filtered.saige.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_singles_script} \
          -p ${pheno} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        touch ${pheno}.gene_phewas_singles.saige.gz
        touch ${pheno}.gene_phewas_singles.filtered.saige.csv
        """
}

process merge_and_filter_saige_gwas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    memory '35GB'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr_list), path(chr_inputs)

        path merge_regions_script
    output:
        tuple val(cohort_dir), val(pheno), path("${pheno}.gwas.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${pheno}.gwas.filtered.saige.csv")
    shell:
        """
        echo "${params.gwas_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          -p ${pheno} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Sumstats/

        touch ${pheno}.gwas.saige.gz
        touch ${pheno}.gwas.filtered.saige.csv
        """
}


process merge_and_filter_saige_variantphewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr), path(chr_inputs)

        path merge_regions_script
    output:
        tuple val(cohort_dir), val(pheno), path("chr${chr}.variant_phewas.saige.gz")
        tuple val(cohort_dir), val(pheno), path("chr${chr}.variant_phewas.filtered.saige.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          -p ${pheno} \
          -c colnames.txt \
          --chr ${chr} \
          --phewas \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Sumstats/

        touch ${pheno}.variant_phewas.saige.gz
        touch ${pheno}.variant_phewas.filtered.saige.csv
        """
}

process make_summary_regions_output {
    publishDir "${launchDir}/Summary/"

    input:
        path(filtered_regions, stageAs: '?/*')
    output:
        path('saige_exwas_suggestive_regions.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${filtered_regions}".replace('[', '').replace(']', '').split()]
        output = 'saige_exwas_suggestive_regions.csv'

        for f in input_list:
            dfs.append(pd.read_csv(f))

        region_col = '${params.regions_col_names.keySet().toList().contains('Region') ? params.regions_col_names['Region'] : 'Region' }'
        pd.concat(dfs).sort_values(by=[region_col]).to_csv(output, index=False)
        """
    stub:
        """
        mkdir -p ${launchDir}/Summary/

        touch saige_exwas_suggestive_regions.csv
        """
}

process make_summary_singles_output {
    publishDir "${launchDir}/Summary/"

    input:
        path(filtered_singles, stageAs: '?/*')
    output:
        path('saige_exwas_suggestive_singles.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${filtered_singles}".replace("[", "").replace("]", "").split()]
        output = "saige_exwas_suggestive_singles.csv"

        for f in input_list:
            dfs.append(pd.read_csv(f))

        chr_col = '${params.singles_col_names.keySet().toList().contains('CHR') ? params.singles_col_names['CHR'] : 'CHR' }'
        pos_col = '${params.singles_col_names.keySet().toList().contains('POS') ? params.singles_col_names['POS'] : 'POS' }'
        pd.concat(dfs).sort_values(by=[chr_col, pos_col]).to_csv(output, index=False)
        """
    stub:
        """
        mkdir -p ${launchDir}/Summary/

        touch saige_exwas_suggestive_singles.csv
        """
}

process gwas_make_biofilter_positions_input {
    publishDir "${launchDir}/Annotations/"

    input:
        path(filtered_sumstats)
    output:
        path('gwas_biofilter_input_positions.txt')
    script:
        """
        #! ${params.my_python}
        import pandas as pd

        dfs = []
        input_list = '${filtered_sumstats.join(' ')}'.split()
        for f in input_list:
            temp = pd.read_table(f)
            temp['chromosome_noCHR'] = temp['chromosome'].str.replace('chr', '').astype(int)
            temp = temp.rename(columns={'chromosome_noCHR': 'CHR', 'base_pair_location': 'POS', 'p_value': 'P', 'variant_id': 'ID'})
            dfs.append(pd.read_table(f))
        all = pd.concat(dfs)
        keep_cols = ['CHR', 'variant_id', 'POS']
        all[keep_cols].to_csv('gwas_biofilter_input_positions.txt', header=False, index=False, sep=' ')
        """
    stub:
        '''
        touch gwas_biofilter_input_positions.txt
        '''
}

process make_summary_table_with_annot {
    publishDir "${launchDir}/Summary/"

    input:
        path all_filtered_sumstats
        tuple val(data_nickname), path(biofilter_annots)
    output:
        path('meta_all_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_filtered_sumstats.join(' ')}'.split()
        output = "meta_all_suggestive.csv"

        for f in input_list:
            dfs.append(pd.read_table(f))

        all_meta = pd.concat(dfs).sort_values(by='variant_id')

        annot_df = pd.read_csv('${biofilter_annots}', index_col='Var_ID')
        all_meta[['Gene', 'RSID']] = annot_df.loc[all_meta['variant_id'], ['Gene', 'RSID']].values


        all_meta.to_csv(output, index=False)
        """
    stub:
        '''
        touch meta_all_suggestive.csv
        '''
}