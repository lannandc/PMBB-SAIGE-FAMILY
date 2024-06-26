log.info """\
    NEXTFLOW - DSL2 - SAIGE ExWAS - P I P E L I N E
    ==================================================
    run as                  : ${workflow.commandLine}
    run location            : ${launchDir}
    started at              : ${workflow.start}
    python exe              : ${params.my_python}

    Cohorts, Phenotypes, and Chromosomes
    ==================================================
    cohort_list             : ${params.cohort_list}
    sex-stratified_cohorts  : ${params.sex_strat_cohort_list}
    bin_pheno_file          : ${params.bin_pheno_list_file}
    quant_pheno_file        : ${params.quant_pheno_list_file}
    gene_list_file          : ${params.gene_list_file}
    chromosome_list         : ${params.chromosome_list}
    cat_covars              : ${params.cat_covars}
    cont_covars             : ${params.cont_covars}
    data_csv                : ${params.data_csv}
    GPU                     : ${params.GPU}
    cohort_sets             : ${params.cohort_sets}

    Input File Prefixes
    ==================================================
    step1_plink_prefix      : ${params.step1_plink_prefix}
    step1_sparse_grm        : ${params.step1_sparse_grm}
    step1_sparse_grm_samples: ${params.step1_sparse_grm_samples}
    sparse_grm              : ${params.use_sparse_GRM}
    step2_plink_prefix      : ${params.plink_prefix}
    group_file_prefix       : ${params.group_file_prefix}
    gene_location_file      : ${params.gene_location_file}

    SAIGE Step 1 Plink QC Parameters
    ==================================================
    min maf (maf)           : ${params.maf}
    max missingness (geno)  : ${params.geno}
    hardy-weinberg (hwe)    : ${params.hwe}

    SAIGE-GENE Parameters
    ==================================================
    minMAF                  : ${params.min_maf}
    minMAC                  : ${params.min_mac}
    maxMAF_in_groupTest     : ${params.grouptest_maf}
    annotation_in_groupTest : ${params.grouptest_annotation}
    is_Firth_beta           : ${params.use_firth}
    pCutoffforFirth         : ${params.firth_cutoff}
    LOCO                    : ${params.LOCO}

    Other Parameters
    ==================================================
    step1_script            : ${params.step1_script}
    step2_script            : ${params.step2_script}
    pheno_file_id_col       : ${params.id_col}
    p_cutoff_summarize      : ${params.p_cutoff_summarize}

    """.stripIndent()

include {
    set_up_cohort
    make_pheno_summaries
    check_bin_cohort_pheno_combo
    check_quant_cohort_pheno_combo
    } from '../processes/saige_preprocessing.nf'

include {
    plink_qc_for_step1
    call_saige_step1_bin
    call_saige_step1_bin_with_sparse_GRM
    call_saige_step1_quant
    call_saige_step1_quant_with_sparse_GRM
    } from '../processes/saige_step1.nf'

include {
    call_saige_step2_BGEN_binary
    call_saige_step2_BGEN_quant
    call_saige_step2_PLINK_binary
    call_saige_step2_PLINK_quant
    filter_snps_plink
    filter_snps_bgen
    } from '../processes/saige_variant_step2.nf'

include {
    merge_and_filter_saige_gene_regions_output
    merge_and_filter_saige_gene_singles_output
    merge_and_filter_saige_gene_singles_phewas_output
    merge_and_filter_saige_variantphewas_output
    make_summary_regions_output
    make_summary_singles_output
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots
    make_saige_exwas_regions_plots
    make_saige_exwas_singles_plots
    make_exwas_report_src
    make_exwas_report
    make_exwas_report_methods_blurb
    } from '../processes/saige_visualization.nf'

workflow {
    // user specified file type either PLINK or BGEN
    ftype = params.ftype
    // cohort channel has all cohorts
    // user specified file type either PLINK or BGEN
    ftype = params.ftype
    // cohort channel has all cohorts
    cohort = Channel.fromList(params.cohort_list)
    // chromosome channel has the relevant chromosomes (only the needed ones depending on genes)
    chromosome = Channel.fromList(params.chromosome_list)

    // read list of binary phenos then make a binary pheno channel
    bin_pheno_list = (new File(params.bin_pheno_list_file)).readLines()
    bin_pheno = Channel.fromList(bin_pheno_list)
     MIN_BIN_CASES = 50

    // read list of quantitative phenos then make a quantitative pheno channel
    quant_pheno_list = (new File(params.quant_pheno_list_file)).readLines()
    quant_pheno = Channel.fromList(quant_pheno_list)

    // since we use hard-calls with ExWAS, we assume plink 1.9 format
    plink_suffixes_list = ['.bed', '.bim', '.fam']
    step1_fam = params.step1_plink_prefix + '.fam'
    step2_fam = params.plink_prefix + '.fam'

    // Input files for cohort set-up process
    cohort_setup_script = "${launchDir}/scripts/set_up_cohort_directory.py"
    pheno_covar_table = "${params.data_csv}"
    cohort_table = "${params.cohort_sets}"

    // Cohort set-up process which returns sample lists and phenotype input tables
    (cohort_sample_lists, cohort_pheno_tables) = set_up_cohort(
        cohort, cohort_setup_script, pheno_covar_table, cohort_table, step1_fam, step2_fam
        )

    // Input files for summary table process
    pheno_table_script = "${launchDir}/scripts/make_pheno_summary_table.py"
    pheno_covar_plots_script = "${launchDir}/scripts/make_pheno_covar_summary_plots.py"
    
    // make pheno summary table, conditionally handle empty phenotype lists
    pheno_table = make_pheno_summaries(
            cohort.collect(),
            (bin_pheno_list.size() == 0)  ? '[]' : bin_pheno.toSortedList(),
            (quant_pheno_list.size() == 0) ? '[]' : quant_pheno.toSortedList(),
            step1_fam, step2_fam,
            pheno_covar_table, cohort_table,
            pheno_table_script
            )

    pheno_plots = make_pheno_covar_summary_plots(
            cohort.collect(),
            (bin_pheno_list.size() == 0)  ? '[]' : bin_pheno.toSortedList(),
            (quant_pheno_list.size() == 0) ? '[]' : quant_pheno.toSortedList(),
            step1_fam, step2_fam,
            pheno_covar_table, cohort_table,
            pheno_covar_plots_script
            )

    // Filter out binary cohort-phenotype combinations with too few cases
    all_cohort_bin_pheno_combos = cohort.combine(bin_pheno)
    cohort_bin_pheno_case_ct = check_bin_cohort_pheno_combo(all_cohort_bin_pheno_combos, pheno_table)
    keep_cohort_bin_pheno_combos = cohort_bin_pheno_case_ct.filter { cohort, pheno, cases -> cases >= 50 }.map { cohort, pheno, cases -> new Tuple(cohort, pheno) }

    // Filter out quantitative cohort-phenotype combinations with too few samples
    all_cohort_quant_pheno_combos = cohort.combine(quant_pheno)
    cohort_quant_pheno_ct = check_quant_cohort_pheno_combo(all_cohort_quant_pheno_combos, pheno_table)
    keep_cohort_quant_pheno_combos = cohort_quant_pheno_ct.filter { cohort, pheno, count -> count >= 500 }.map { cohort, pheno, count -> new Tuple(cohort, pheno) }

    if (!params.use_sparse_GRM) {
        // set up inputs for plink_qc
        plink_inputs_tuple = new Tuple(*plink_suffixes_list.collect { ext -> params.step1_plink_prefix + ext })
        // call plink_qc for each cohort
        (cohort_plinkset, logs) = plink_qc_for_step1(cohort_sample_lists, params.step1_plink_prefix, plink_inputs_tuple)
        // qc_output.view()

        /*
        QC -> Step 1 Channel Emission Tuples
        QC Out:      cohort, (plink bed, bim, fam)
        Combine:     phenotype
        Step 1 In:   cohort, (plink bed, bim, fam), phenotype
        */
        cohort_plinkset_sample_pheno = cohort_plinkset.join(cohort_sample_lists).join(cohort_pheno_tables)

        step1_bin_input = cohort_plinkset_sample_pheno.combine(bin_pheno).map { cohort, plink, sample_file, pheno_file, pheno -> new Tuple(cohort, pheno, plink, sample_file, pheno_file,params.step) }.join(keep_cohort_bin_pheno_combos, by: [0, 1])
        step1_quant_input = cohort_plinkset_sample_pheno.combine(quant_pheno).map { cohort, plink, sample_file, pheno_file, pheno -> new Tuple(cohort, pheno, plink, sample_file, pheno_file) }.join(keep_cohort_quant_pheno_combos, by: [0, 1])

        (step1_bin_output, logs) = call_saige_step1_bin(step1_bin_input)
        (step1_quant_output, logs) = call_saige_step1_quant(step1_quant_input)
    }
    else {
        cohort_plinkset = cohort.map { cohort -> \
            new Tuple(cohort, new Tuple(*plink_suffixes_list.collect { ext -> params.step1_plink_prefix + ext }))
        }

        /*
        Step 1 Input Emission Tuples
        Need: (cohort, phenotype, (plink set), sample_file, pheno_file)
        Have:
        (cohort, (plink_set))
        (cohort, sample_file)
        (cohort, pheno_file)
        (phenotype)
        Steps:
        1. (cohort, (plink_set)).join(cohort, sample_file) -> (cohort, (plink set), sample_file)
        2. .join(cohort, pheno_file) -> (cohort, (plink set), sample_file, pheno_file)
        3. .combine(bin_pheno or quant_pheno) -> (cohort, (plink set), sample_file, pheno_file, phenotype)
        4. .map() -> rearrange to get cohort, pheno first
        5. .join() -> keep only cohort/pheno combos with enough sample size
        */

        cohort_plinkset_sample_pheno = cohort_plinkset.join(cohort_sample_lists).join(cohort_pheno_tables)

        step1_bin_input = cohort_plinkset_sample_pheno.combine(bin_pheno).map { cohort, plink, sample_file, pheno_file, pheno -> new Tuple(cohort, pheno, plink, sample_file, pheno_file) }.join(keep_cohort_bin_pheno_combos, by: [0, 1])
        step1_quant_input = cohort_plinkset_sample_pheno.combine(quant_pheno).map { cohort, plink, sample_file, pheno_file, pheno -> new Tuple(cohort, pheno, plink, sample_file, pheno_file) }.join(keep_cohort_quant_pheno_combos, by: [0, 1])

        sparse_grm_input = new Tuple(params.step1_sparse_grm, params.step1_sparse_grm_samples)

        (step1_bin_output, logs) = call_saige_step1_bin(step1_bin_input)
        (step1_quant_output, logs) = call_saige_step1_quant(step1_quant_input)
    }

    /*
    Step 1 -> Step 2 Channel Emission Tuples
    Step 1 Out:  cohort, phenotype, rda, var
    Combine:     chromosome
    Step 2 In:   cohort, phenotype, rda, var, chromosome
    */

     /*
    Step 1 -> Step 2 Channel Emission Tuples
    Step 1 Out:  cohort, phenotype, rda, var
    Combine:     chromosome
    Step 2 In:   cohort, phenotype, rda, var, chromosome
    */
    step2_bin_input = step1_bin_output.combine(chromosome)
    
    step2_quant_input = step1_quant_output.combine(chromosome)


    if(ftype.equals("PLINK")){
      //pulling out cohort, etc.
      filter_snps_plink(params.snplist)

      step2_prefix = "${launchDir}/${params.step2_input_prefix}"

      genome_plink_file_bin = step2_bin_input.map {
        cohort, pheno, rda, var, chr -> \
        new Tuple(*plink_suffixes_list.collect { ext -> step2_prefix + ext })
      } 
      step2_bin_output = call_saige_step2_PLINK_binary(step2_bin_input, genome_plink_file_bin)
      
      genome_plink_file_quant = step2_quant_input.map {
        cohort, pheno, rda, var, chr -> \
        new Tuple(*plink_suffixes_list.collect { ext -> step2_prefix + ext })
      }
      step2_quant_output = call_saige_step2_PLINK_quant(step2_quant_input, genome_plink_file_quant)
    }
 
    else if(ftype.equals("BGEN")){

      filter_snps_bgen(params.bgen_prefix, params.snplist)

      step2_prefix = "${launchDir}/${params.step2_input_prefix}"

      genome_bgen_file_bin = step2_bin_input.map {
        cohort, pheno, rda, var, chr -> \
        new Tuple(*bgen_suffixes_list.collect { ext -> step2_prefix + ext })
        }

      step2_bin_output = call_saige_step2_BGEN_binary(step2_bin_input,genome_bgen_file_bin)
       
       
      genome_bgen_file_quant = step2_quant_input.map {
        cohort, pheno, rda, var, chr -> \
          new Tuple(*bgen_suffixes_list.collect { ext -> params.step2_prefix + ext })
      }
      
      step2_quant_output = call_saige_step2_BGEN_quant(step2_quant_input,genome_bgen_file_quant)
    }
    //********************Improper File Type****************************
    else{
      throw new Exception("Improper file type for step 2, please refer to your .config file \
                           and ensure the ftype is PLINK, or BGEN.")
    }
    //****************************************************

    /*
    Step 2 -> Merged Sumstats Channel Emission Tuples
    Step 2 Out:  cohort, phenotype, chromosome, regions, singles
    Group By:    cohort, phenotype
    Merge In:    cohort, phenotype, [chr_list], [region_list], [singles_list]
      - then map to split singles vs regions
    */
    // Collect saige output into channels for merge
    step2_all_output = step2_bin_output.concat(step2_quant_output)
    step2_grouped_output = step2_all_output.groupTuple(by: [0, 2]) 
    merge_singles_script = "${launchDir}/scripts/merge_and_filter_saige_variant.py"
    (singles_merge_output, filtered_singles_output) = merge_and_filter_saige_variantphewas_output(step2_grouped_output, merge_singles_script)

    // collect a list of just the filtered output files, don't need a wildcards anymore
    summary_singles_input = filtered_singles_output.map { cohort, pheno, filtered -> filtered }.collect()
    singles_summary = make_summary_singles_output(summary_singles_input)

/*
    // plots and report post-processing
    gene_file = "${params.gene_location_file}"
    singles_plots_script = "${launchDir}/scripts/make_saige_exwas_singles_plots.py"
    singles_plots = make_saige_exwas_singles_plots(singles_merge_output, gene_file, singles_plots_script, pheno_table)
*/
}
