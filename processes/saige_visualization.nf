import groovy.json.JsonBuilder

process make_pheno_covar_summary_plots {
    publishDir "${launchDir}/Plots/"

    input:
        val cohort_list
        val bin_pheno_list
        val quant_pheno_list
        path(step1_fam, stageAs: 'Plink_QC/*')
        path(exome_fam, stageAs: 'Exome/*')
        path pheno_covar_table
        path cohort_table

        path(pheno_covar_plots_script)
    output:
        path('*.png')
    shell:
        """
        ${params.my_python} ${pheno_covar_plots_script} \
          -c ${cohort_list.join(' ')} \
          -b ${bin_pheno_list.join(' ')} \
          -q ${quant_pheno_list.join(' ')} \
          --step1Fam ${step1_fam} \
          --exomeFam ${exome_fam} \
          --data ${pheno_covar_table} \
          --samples ${cohort_table} \
          --id ${params.id_col}
        """
    stub:
        '''
        touch stub_plot.png
        '''
}

process make_saige_gwas_plots {
    publishDir "${launchDir}/Plots/"
    memory '35GB'
    input:
        // from singles_merge_output
        tuple val(cohort), val(pheno), path(sumstats)
        //
        path(gwas_plot_script)

        path(pheno_table)
    output:
        path("${cohort}.${pheno}*.png")
    shell:
        """
        echo "${params.gwas_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${gwas_plot_script} \
          -c ${cohort} \
          -p ${pheno} \
          -s ${sumstats} \
          -t ${pheno_table}
        """
    stub:
        """
        touch ${cohort}.${pheno}.singles.manhattan_vertical.png
        touch ${cohort}.${pheno}.singles.qq.png
        """
}

process make_saige_exwas_singles_plots {
    publishDir "${launchDir}/Plots/"

    input:
        tuple val(cohort), val(pheno), path(singles_sumstats)

        path(gene_file)
        path(exwas_singles_plot_script)

        path(pheno_table)
    output:
        path("${cohort}.${pheno}*.png")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${exwas_singles_plot_script} \
          -c ${cohort} \
          -p ${pheno} \
          -g ${gene_file} \
          -s ${singles_sumstats} \
          -t ${pheno_table} \
          -m ${params.grouptest_maf} \
          -a "${params.grouptest_annotation}"
        """
    stub:
        """
        touch ${cohort}.${pheno}.singles.manhattan_vertical.png
        touch ${cohort}.${pheno}.singles.qq.png
        """
}

process make_saige_exwas_regions_plots {
    publishDir "${launchDir}/Plots/"

    input:
        tuple val(cohort), val(pheno), path(regions_sumstats)

        path(gene_file)
        path(exwas_regions_plot_script)

        path(pheno_table)
    output:
        path("${cohort}.${pheno}*.png")
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${exwas_regions_plot_script} \
          -c ${cohort} \
          -p ${pheno} \
          -g ${gene_file} \
          -s ${regions_sumstats} \
          -t ${pheno_table} \
          -m ${params.grouptest_maf} \
          -a "${params.grouptest_annotation}"
        """
    stub:
        """
        touch ${cohort}.${pheno}.regions.manhattan.png
        touch ${cohort}.${pheno}.regions.qq.png
        """
}

process make_saige_gene_phewas_regions_plots {
    publishDir "${launchDir}/Plots/"

    input:
        tuple val(cohort), val(chr), path(sumstats_file)
        path phenotype_descriptions
        path phewas_regions_plot_script
    output:
        path('*.png')
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${phewas_regions_plot_script} \
          -d ${phenotype_descriptions} \
          -s ${sumstats_file} \
          -c ${cohort} \
          --cols colnames.txt
        """
    stub:
        '''
        touch test.png
        '''
}

process make_exwas_report_methods_blurb {
    publishDir "${launchDir}/Summary"

    input:
        path exwas_methods_script
    output:
        path('saige_exwas_methods.html')
    shell:
        """
        echo '${new JsonBuilder(params).toPrettyString().replace(';', '|')}' > pipeline_params.txt
        ${params.my_python} ${exwas_methods_script}
        """
    stub:
        '''
        touch saige_exwas_methods.html
        '''
}

process make_exwas_report_src {
    publishDir "${launchDir}/Report/", mode: 'copy', overwrite: true

    input:
        path singles_plots
        path regions_plots
        path pheno_summary_plots
        path pheno_summary_table
        path singles_summary
        path regions_summary
        path methods_blurb
    output:
        path('src/', type: 'dir')
    shell:
        """
        mkdir src/

        for f in \$(echo "${singles_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
        echo \$f
        cp \$f src/
        done

        for f in \$(echo "${regions_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
        echo \$f
        cp \$f src/
        done

        for f in \$(echo "${pheno_summary_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
        echo \$f
        cp \$f src/
        done

        cp ${pheno_summary_table} src/
        cp ${singles_summary} src/
        cp ${regions_summary} src/

        cp ${methods_blurb} src/

        """
    stub:
        '''
        mkdir src
        touch src/test.png
        touch src/test.csv
        '''
}

process make_exwas_report {
    publishDir "${launchDir}/Report/", mode: 'copy', overwrite: true

    input:
        path report_source_dir
        path generate_html_script
    output:
        path('index.html')
        path('*.html')
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > regions_colnames.txt
        echo "${params.singles_col_names.collect().join('\n')}" > singles_colnames.txt
        ${params.my_python} ${generate_html_script} \
          -p ${params.bin_pheno_list.join(' ') + ' ' + params.quant_pheno_list.join(' ')} \
          -c ${params.cohort_list.join(' ')} \
          --pval ${params.p_cutoff_summarize} \
          --regionsColnames regions_colnames.txt \
          --singlesColnames singles_colnames.txt
        """
    stub:
        '''
        touch index.html
        '''
}
process make_gwas_plots_with_annot {
    publishDir "${launchDir}/Plots/"

    input:
        tuple val(pheno), path(sumstats), val(data_nickname), path(biofilter_annots)
        path plotting_script
    output:
        path "${pheno}.${analysis}.{manhattan.png,qq.png,qq.csv}"
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --pheno ${pheno} \
          --analysis ${analysis} \
          --sumstats ${sumstats} \
          --annot ${biofilter_annots}
        """
    stub:
        """
        touch ${pheno}.${analysis}.manhattan.png
        touch ${pheno}.${analysis}.qq.png
        touch ${pheno}.${analysis}.qq.csv
        """
}