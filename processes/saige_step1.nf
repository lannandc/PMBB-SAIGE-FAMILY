
process plink_qc_for_step1 {
    publishDir "${cohort_dir}/Saige_Step1/"

    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)
    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'
    shell:
        """
        stdbuf -e0 -o0 plink --make-bed \
          --bfile ${input_bed.toString().replace('.bed', '')} \
          --keep-fam ${sample_list} \
          --maf ${params.maf} --geno ${params.geno} --not-chr X,Y,MT --hwe ${params.hwe} --snps-only \
          --out step1_plink > plink_qc.log
        """
    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

process plink_qc_for_step1_saige_gene {
    publishDir "${cohort_dir}/Saige_Step1/"
    /* commmon variants : ld pruning 50 5 0.4
    / 10-20 : plink --counts --snplist -> 1000
    / 20-430 : plink --counts --snplist -> 1000
    / plink --extract
    */
    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)
    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'
    shell:
        """

        # get counts from whole plink file
        plink --freq counts \
          --keep-fam ${sample_list} \
          --bfile ${input_bed.toString().replace('.bed', '')} \
          --out step1_counts

        awk '{ if ((\$5 >= 10 && \$5 <= 20) || (2*\$6-\$5 >= 10 && 2*\$6-\$5 <= 20)) { print \$2} }' ./step1_counts.acount > temp_mac10to20_ids.txt
        awk '{ if ((\$5 >= 20) || (2*\$6-\$5 >= 20)) { print \$2} }' step1_counts.acount > temp_mac20ormore_ids.txt

        num_mac10=\$(< temp_mac10to20_ids.txt wc -l)
        num_mac20=\$(< temp_mac20ormore_ids.txt wc -l)

        if [ "\$num_mac10" -lt 1000 ]; then
            echo "Error: Less than 1000 SNPs with MAC between 10 and 20"
        elif [ "\$num_mac20" -lt 1000 ]; then
            echo "Error: Less than 1000 SNPs with MAC greater than 20"
        else
            # Copy the filtered IDs to the output file
            shuf -n 1000 temp_mac20ormore_ids.txt > ./step1.markerid.list
            shuf -n 1000 temp_mac10to20_ids.txt >> ./step1.markerid.list
            echo "Filtered IDs written to step1.markerid.list"
        fi

        rm temp_mac10to20_ids.txt

        # ld pruning,maf filter, geno filter, hwe filter and write variant ids to step1.prune.in
        plink --indep-pairwise 50 5 0.4 \
          --bfile ${input_bed.toString().replace('.bed', '')} \
          --keep-fam ${sample_list} \
          --maf ${params.maf} --geno ${params.geno} --not-chr X,Y,MT --hwe ${params.hwe} --snps-only \
          --out step1

        # create snp list with variants from all categories

        shuf -n 200000 step1.prune.in >> step1.markerid.list

        #  make input file for step 1 using variant list
        stdbuf -e0 -o0 plink --make-bed \
          --bfile ${input_bed.toString().replace('.bed', '')} \
          --keep-fam ${sample_list} \
          --extract step1.markerid.list \
          --not-chr X,Y,MT \
          --snps-only \
          --out step1_plink > plink_qc.log

        rm step1_counts.acount
        rm step1.prune.in
        """
    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

String get_covar_list_args(String cohort, cohort_cat_covars, cohort_cont_covars) {
    String output = ''

    if (cohort_cont_covars.size() > 0 || cohort_cat_covars.size() > 0) {
        output += '--covarColList='
        if (cohort_cont_covars.size() > 0) {
            output += cohort_cont_covars.join(',')
        }
        if (cohort_cat_covars.size() > 0) {
            output += ',' + cohort_cat_covars.join(',') + ' '
        }
        else {
            output += ' '
        }
    }

    if (cohort_cat_covars.size() > 0) {
        output += '--qCovarColList=' + cohort_cat_covars.join(',') + ' '
    }

    return output
}

process call_saige_step1_bin {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"

    cpus 15

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --plinkFile=step1_plink \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --SampleIDIncludeFile=${samples} \
         --traitType=binary \
         --outputPrefix=${pheno} \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.log
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        """
}

process call_saige_step1_quant {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"

    cpus 15

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --plinkFile=step1_plink \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --SampleIDIncludeFile=${samples} \
         --outputPrefix=${pheno} \
         --traitType=quantitative \
         --invNormalize=TRUE \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.log
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        """
}

process call_saige_step1_bin_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"

    cpus 1

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --plinkFile=${params.step1_plink_prefix.split('/')[-1]} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --SampleIDIncludeFile=${samples} \
         --traitType=binary \
         --outputPrefix=${pheno} \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"

    cpus 1

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --plinkFile=${params.step1_plink_prefix.split('/')[-1]} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --SampleIDIncludeFile=${samples} \
         --outputPrefix=${pheno} \
         --traitType=quantitative \
         --invNormalize=TRUE \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}
