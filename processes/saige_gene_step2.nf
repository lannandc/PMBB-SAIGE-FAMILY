process filter_chr_group_file {
    publishDir "${launchDir}/PheWAS_GroupFiles/"

    input:
        tuple val(chr), path(chr_group_file)
        path(gene_list_file)
    output:
        tuple val(chr), path("filtered_group_file.${chr}.txt")
    script:
        """
        #! ${params.my_python}

        keep_genes = open('${gene_list_file}').read().splitlines()

        keep_lines = []

        with open('${chr_group_file}') as f:
            for line in f:
                if line.split('\\t', 2)[0] in keep_genes:
                    keep_lines.append(line.strip('\\n'))
        
        open('filtered_group_file.${chr}.txt', 'w+').write('\\n'.join(keep_lines) + '\\n')
        """
    stub:
        """
        touch filtered_group_file.${chr}.txt
        """
}
process filter_gene_file_plink {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"

    input:
        tuple val(chr), path(chr_group_file)
        path(gene_list_file)
    output:
        tuple val(chr), path("_file.${chr}.txt")
    script:
        """
        #! ${params.my_python}

        keep_genes = open('${gene_list_file}').read().splitlines()

        keep_lines = []

        with open('${chr_group_file}') as f:
            for line in f:
                if line.split('\\t', 2)[0] in keep_genes:
                    keep_lines.append(line.strip('\\n'))
        
        open('filtered_group_file.${chr}.txt', 'w+').write('\\n'.join(keep_lines) + '\\n')
        """
    stub:
        """
        touch filtered_group_file.${chr}.txt
        """

}

process call_saige_gene_step2_bin {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // parameters
        val exome_prefix
        val group_file_prefix

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}"), path("${pheno}.${chr}.singleAssoc.txt"), path("${pheno}.${chr}.markerList.txt")
    shell:
        """

        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --is_Firth_beta=${params.use_firth} \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
         --SAIGEOutputFile=${pheno}.${chr} \
           > ${launchDir}/${cohort_dir}/Saige_Gene_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        touch ${pheno}.${chr}.markerList.txt
        """
}

process call_saige_gene_step2_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // parameters
        val exome_prefix
        val group_file_prefix

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}"), path("${pheno}.${chr}.singleAssoc.txt"), path("${pheno}.${chr}.markerList.txt") 
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --SAIGEOutputFile=${pheno}.${chr} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
           > ${launchDir}/${cohort_dir}/Saige_Gene_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        touch ${pheno}.${chr}.markerList.txt
        """
}

process call_saige_gene_step2_bin_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // parameters
        val exome_prefix
        val group_file_prefix

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}"), path("${pheno}.${chr}.singleAssoc.txt"), path("${pheno}.${chr}.markerList.txt")
    shell:
        """

        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --is_Firth_beta=${params.use_firth} \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
         --SAIGEOutputFile=${pheno}.${chr} \
           > ${launchDir}/${cohort_dir}/Saige_Gene_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        touch ${pheno}.${chr}.markerList.txt
        """
}

process call_saige_gene_step2_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // parameters
        val exome_prefix
        val group_file_prefix

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}"), path("${pheno}.${chr}.singleAssoc.txt"), path("${pheno}.${chr}.markerList.txt") 
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --SAIGEOutputFile=${pheno}.${chr} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
           > ${launchDir}/${cohort_dir}/Saige_Gene_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Gene_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        touch ${pheno}.${chr}.markerList.txt
        """
}