/*
This is a collection of SAIGE Step 2 processes
*/

process filter_snps_plink {
    publishDir "${launchDir}"

    input:
        path(snplist)
    output:
      path("${params.step2_input_prefix}.{bed,bim,fam}")
    shell:
        """
        plink --bfile ${params.step1_plink_prefix} \
          --extract range ${params.snplist} \
          --make-bed  \
          --out ${params.step2_input_prefix} 
        """
    stub:
        """
        mkdir -p ${launchDir}/Filtered_Files
        touch ${launchDir}/Filtered_Files/${params.step2_input_prefix}.bed
        touch ${launchDir}/Filtered_Files/${params.step2_input_prefix}.bim
        touch ${launchDir}/Filtered_Files/${params.step2_input_prefix}.fam
        """
}


process filter_snps_bgen {
    input:
        path(bgen_prefix)
        path(snplist)
    output:
        path('${launchDir}/Filtered_Files/${params.step2_input_prefix}.bgen')
    shell:
        """
        mkdir -p ${launchDir}/Filtered_Files

        awk -F ',' '{print ${1} ":" ${2} "-" ${3}}' ${params.snplist} > bgen_formatted_snplist.txt
        bgenix -g ${params.bgen_file} -incl-range bgen_formatted_snplist.txt > ${launchDir}/Filtered_Files/${params.step2_input_prefix}.bgen
        bgenix -g ${launchDir}/Filtered_Files/${params.step2_input_prefix}.bgen -index
        """
    stub:
        """
        mkdir -p ${launchDir}/Filtered_Files
        touch ${launchDir}/Filtered_Files/${params.step2_input_prefix}.bgen
        touch ${launchDir}/Filtered_Files/${params.step2_input_prefix}.bgi
        """

}


process call_saige_step2_PLINK_binary {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}")
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --is_output_moreDetails=TRUE \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=0.05 \
         --LOCO=TRUE \
         --SAIGEOutputFile=${pheno}.${chr} \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_PLINK_quant {
  publishDir "${launchDir}/${cohort_ßdir}/Saige_Step2_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}")
    shell:
        """
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${pheno}.${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_PLINK_binary_with_sparse_GRM {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}")
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

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
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=0.05 \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${pheno}.${chr} \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_PLINK_quant_with_sparse_GRM {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"

    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}")
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

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
         --SAIGEOutputFile=${pheno}.${chr} \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
          > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/
        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_VCF_binary {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        // inputs
        tuple path(vcfFile), path(vcfFileIndex)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --vcfFile=${vcfFile}      \
         --vcfFileIndex=${vcfFileIndex}       \
         --vcfField=${params.vcfField}       \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_VCF_quant{
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        // inputs
        tuple path(vcfFile), path(vcfFileIndex)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --vcfFile=${vcfFile}      \
         --vcfFileIndex=${vcfFileIndex}       \
         --vcfField=${params.vcfField}       \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
         --LOCO=FALSE \
         --is_output_moreDetails=FALSE \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_VCF_binary_with_sparse_GRM {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        // inputs
        tuple path(vcfFile), path(vcfFileIndex)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --vcfFile=${vcfFile}      \
         --vcfFileIndex=${vcfFileIndex}       \
         --vcfField=${params.vcfField}       \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_VCF_quant_with_sparse_GRM{
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        // inputs
        tuple path(vcfFile), path(vcfFileIndex)
        tuple path(sparse_grm), path(sparse_grm_samples)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --vcfFile=${vcfFile}      \
         --vcfFileIndex=${vcfFileIndex}       \
         --vcfField=${params.vcfField}       \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
         --LOCO=FALSE \
         --is_output_moreDetails=FALSE \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_BGEN_binary {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        //inputs
        tuple path(bgenFile), path(bgenFileIndex)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        echo "GOT HERE" > s2test.txt
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bgenFile=${bgenFile} \
         --bgenFileIndex=${bgenFileIndex} \
         --AlleleOrder=ref-first \
         --sampleFile=${params.samplefile} \
         --chrom=chr${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}.txt
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_BGEN_quant {
 publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
      
        tuple path(bgenFile), path(bgenFileIndex)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        echo ${step1_rda}

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bgenFile=${bgenFile}    \
         --bgenFileIndex=${bgenFileIndex} \
         --sampleFile=${params.samplefile} \
         --chrom=chr${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
         --LOCO=FALSE \
         --is_output_moreDetails=TRUE \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}.txt
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_BGEN_binary_with_sparse_GRM {
  publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        //inputs
        tuple path(bgenFile), path(bgenFileIndex)
        tuple path(sparse_grm), path(sparse_grm_samples)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bgenFile=${bgenFile} \
         --bgenFileIndex=${bgenFileIndex} \
         --AlleleOrder=ref-first \
         --sampleFile=${params.samplefile}
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}

process call_saige_step2_BGEN_quant_with_sparse_GRM {
 publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
      
        tuple path(bgenFile), path(bgenFileIndex)
        tuple path(sparse_grm), path(sparse_grm_samples)
    
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${pheno}.${chr}.txt")
    
    shell:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bgenFile=${bgenFile}    \
         --bgenFileIndex=${bgenFileIndex} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${pheno}.${chr}.txt \
         --LOCO=TRUE \
         --is_output_moreDetails=TRUE \
           > ${launchDir}/${cohort_dir}/Saige_Step2_Results/${pheno}.${chr}.log
        """
    
    stub:
        """
        mkdir -p ${launchDir}/${cohort_dir}/Saige_Step2_Results/

        touch ${pheno}.${chr}
        touch ${pheno}.${chr}.singleAssoc.txt
        """
}