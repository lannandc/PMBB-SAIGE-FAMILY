process set_up_cohort {
    publishDir "${launchDir}/${cohort}/"

    input:
        val cohort
        path cohort_script
        path pheno_covar_table
        path cohort_table
        path(step1_fam, stageAs: 'Plink_QC/*')
        path(exome_fam, stageAs: 'Exome/*')
    output:
        tuple val(cohort), path('sample_list.txt')
        tuple val(cohort), path('saige_pheno_covars.txt')
    shell:
        """
        ${params.my_python} ${cohort_script} \
          --data ${pheno_covar_table} \
          --cohort ${cohort} \
          --samples ${cohort_table} \
          --step1Fam ${step1_fam} \
          --exomeFam ${exome_fam} \
          --id ${params.id_col}
        """
    stub:
        """
        touch sample_list.txt
        touch saige_pheno_covars.txt
        """
}

process check_needed_scripts {
    shell:
        """
        ls ${params.step1_script}
        ls ${params.step2_script}
        """
}
//TODO change exome to genetic data label for readability
process make_pheno_summaries {
    publishDir "${launchDir}/Summary/"

    input:
        val cohort_list
        val bin_pheno_list
        val quant_pheno_list
        path(step1_fam, stageAs: 'Plink_QC/*')
        path(exome_fam, stageAs: 'Exome/*')
        path pheno_covar_table
        path cohort_table
        

        path(pheno_table_script)
    output:
        path('pheno_summaries.csv')

    shell:
        """
        ${params.my_python} ${pheno_table_script} \
          -c ${cohort_list.join(' ')} \
          -b ${bin_pheno_list.join(' ')} \
          -q ${quant_pheno_list.join(' ')} \
          --step1Fam ${step1_fam} \
          --exomeFam ${exome_fam} \
          --data ${pheno_covar_table} \
          --samples ${cohort_table} \
          --id ${params.id_col}
        """
}

process check_bin_cohort_pheno_combo {
    input:
        tuple val(cohort), val(pheno)
        val(pheno_table)
    output:
        tuple val(cohort), val(pheno), val(num_cases)
    script:
        // println ("${cohort}_${pheno}")
        pheno_count_file = new File(pheno_table.toString())
        rows = pheno_count_file.readLines().tail()*.split(',')
        rows = rows.findAll{it[0] == cohort}
        rows = rows.findAll{it[1] == pheno}
        match_row = rows.get(0)
        str_num_cases = match_row[4]
        // str_num_cases = rows.get(0)[4]
        // num_cases = str_num_cases == '' ? 0 : str_num_cases.toDouble()
        num_cases = str_num_cases == '' ? 0 : str_num_cases.toDouble()
        """
        echo "${cohort} ${pheno}"
        """
}

process check_quant_cohort_pheno_combo {
    input:
        tuple val(cohort), val(pheno)
        val(pheno_table)
    output:
        tuple val(cohort), val(pheno), val(num_cases)
    script:
        pheno_count_file = new File(pheno_table.toString())
        rows = pheno_count_file.readLines().tail()*.split(',')
        rows = rows.findAll{it[0] == cohort}
        rows = rows.findAll{it[1] == pheno}
        str_num_cases = rows.get(0)[2]
        num_cases = str_num_cases == '' ? 0 : str_num_cases.toDouble()
        """
        """
}
