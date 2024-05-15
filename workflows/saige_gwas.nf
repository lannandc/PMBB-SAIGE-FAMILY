#! /appl/nextflow-23.04.1.5866/bin/nextflow
nextflow.enable.dsl=2

MIN_BIN_CASES = 50
MIN_QUANT_N = 500

/*
Chris Carson modifying and adding to the work of 
Florian Wuennemann and Lindsay Guare for LPC at Penn
*/

log.info """\
    NEXTFLOW - DSL2 - SAIGE GWAS - P I P E L I N E
    ==================================================
    run as                  : ${workflow.commandLine}
    run location            : ${launchDir}
    started at              : ${workflow.start}
    python exe              : ${params.my_python}
    filetype(bgen,plink)    : ${params.ftype}

    Cohorts, Phenotypes, and Chromosomes
    ==================================================
    cohort_list             : ${params.cohort_list}
    bin_pheno_list          : ${params.bin_pheno_list}
    quant_pheno_list        : ${params.quant_pheno_list}
    chromosome_list         : ${params.chromosome_list}
    cat_covars              : ${params.cat_covars}
    cont_covars             : ${params.cont_covars}
    
    Input File Prefixes
    ==================================================
    use_sparse_GRM          : ${params.use_sparse_GRM}
    step1_plink_prefix      : ${params.step1_plink_prefix}
    genome_plink_prefix     : ${params.genome_plink_prefix}
    genome_bgen_prefix      : ${params.genome_bgen_prefix}

    SAIGE Step 1 Plink QC Parameters
    ==================================================
    min maf (maf)           : ${params.maf}
    max missingness (geno)  : ${params.geno}
    hardy-weinberg (hwe)    : ${params.hwe}

    SAIGE-GENE Parameters
    ==================================================
    minMAF                  : ${params.min_maf}
    minMAC                  : ${params.min_mac}
    pCutoffforFirth         : ${params.firth_cutoff}   

    Other Parameters
    ==================================================
    step1_script            : ${params.step1_script}
    step2_script            : ${params.step2_script}
    bgen_samplefile         : ${params.samplefile}
    pheno_file_id_col       : ${params.id_col}
    p_cutoff_summarize      : ${params.p_cutoff_summarize}
    gwas_col_names          : ${params.gwas_col_names}
    annotate                : ${params.annotate}
    biofilter_loki          : ${params.biofilter_loki}
    biofilter_script        : ${params.biofilter_script}

    """.stripIndent()

include {
    set_up_cohort;
    make_pheno_summaries
    check_bin_cohort_pheno_combo
    check_quant_cohort_pheno_combo
    } from '../processes/saige_preprocessing.nf'

include {
    plink_qc_for_step1;
    call_saige_step1_bin;
    call_saige_step1_bin_with_sparse_GRM;
    call_saige_step1_quant;
    call_saige_step1_quant_with_sparse_GRM
    } from '../processes/saige_step1.nf'

include {
    call_saige_gene_step2_bin;
    call_saige_gene_step2_bin_with_sparse_GRM;
    call_saige_gene_step2_quant;
    call_saige_gene_step2_quant_with_sparse_GRM
    } from '../processes/saige_gene_step2.nf'

include{
    call_saige_step2_PLINK_binary;
    call_saige_step2_PLINK_quant;
    call_saige_step2_PLINK_binary_with_sparse_GRM;
    call_saige_step2_PLINK_quant_with_sparse_GRM;
    call_saige_step2_BGEN_binary;
    call_saige_step2_BGEN_quant;
    call_saige_step2_BGEN_binary_with_sparse_GRM;
    call_saige_step2_BGEN_quant_with_sparse_GRM;
} from '../processes/saige_variant_step2.nf'

include {
    merge_and_filter_saige_gene_regions_output;
    merge_and_filter_saige_gene_singles_output;
    merge_and_filter_saige_gwas_output;
    make_summary_singles_output;
    gwas_make_biofilter_positions_input;
    make_summary_table_with_annot;
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots;
    make_saige_gwas_plots;
    make_gwas_plots_with_annot;
    } from '../processes/saige_visualization.nf'

include {BIOFILTER_POSITIONS} from '/project/pmbb_codeworks/projects/geno_pheno_workbench_dev/Biofilter_Wrapper/biofilter_wrapper.nf'

workflow {
    cohort = Channel.fromList(params.cohort_list)
    bin_pheno = Channel.fromList(params.bin_pheno_list)
    quant_pheno = Channel.fromList(params.quant_pheno_list)
    chromosome = Channel.fromList(params.chromosome_list)
    ftype = params.ftype

    plink_suffixes_list = ['.bed', '.bim', '.fam']
    bgen_suffixes_list = ['.bgen','.bgen.bgi']
    cohort_setup_script = "${launchDir}/scripts/set_up_cohort_directory.py"
    pheno_covar_table = "${params.data_csv}"
    cohort_table = "${params.cohort_sets}"

    step1_fam = params.step1_plink_prefix + '.fam'
    genome_fam = params.genome_plink_prefix + '.fam'
    
    //Preprocessing cohort pheno tables
    
    (cohort_sample_lists, cohort_pheno_tables) = set_up_cohort(cohort, cohort_setup_script, pheno_covar_table, cohort_table, step1_fam, genome_fam)

    pheno_table_script = "${launchDir}/scripts/make_pheno_summary_table.py"
    pheno_covar_plots_script = "${launchDir}/scripts/make_pheno_covar_summary_plots.py"

   // make pheno summary table, conditionally handle empty phenotype lists

    pheno_table = make_pheno_summaries(
            cohort.collect(),
            (params.bin_pheno_list.size() == 0)  ? '[]' : bin_pheno.toSortedList(),
            (params.quant_pheno_list.size() == 0)? '[]' : quant_pheno.toSortedList(),
            step1_fam, genome_fam,
            pheno_covar_table, cohort_table,
            pheno_table_script
            )
    pheno_plots = make_pheno_covar_summary_plots(
            cohort.collect(),
            (params.bin_pheno_list.size() == 0)  ? '[]' : bin_pheno.toSortedList(),
            (params.quant_pheno_list.size() == 0)? '[]' : quant_pheno.toSortedList(),
            step1_fam, genome_fam,
            pheno_covar_table, cohort_table,
            pheno_covar_plots_script
            )
    // Filter out binary cohort-phenotype combinations with too few cases
    all_cohort_bin_pheno_combos = cohort.combine(bin_pheno)
    cohort_bin_pheno_case_ct = check_bin_cohort_pheno_combo(all_cohort_bin_pheno_combos, pheno_table)
    keep_cohort_bin_pheno_combos = cohort_bin_pheno_case_ct.filter { cohort, pheno, cases -> cases >= MIN_BIN_CASES }.map { cohort, pheno, cases -> new Tuple(cohort, pheno) }

    // Filter out quantitative cohort-phenotype combinations with too few samples
    all_cohort_quant_pheno_combos = cohort.combine(quant_pheno)
    cohort_quant_pheno_ct = check_quant_cohort_pheno_combo(all_cohort_quant_pheno_combos, pheno_table)
    keep_cohort_quant_pheno_combos = cohort_quant_pheno_ct.filter { cohort, pheno, count -> count >= MIN_QUANT_N }.map { cohort, pheno, count -> new Tuple(cohort, pheno) }

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

        step1_bin_input = cohort_plinkset_sample_pheno.combine(bin_pheno).map { cohort, plink, sample_file, pheno_file, pheno -> new Tuple(cohort, pheno, plink, sample_file, pheno_file) }.join(keep_cohort_bin_pheno_combos, by: [0, 1])
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

        (step1_bin_output, logs) = call_saige_step1_bin_with_sparse_GRM(step1_bin_input, sparse_grm_input)
        (step1_quant_output, logs) = call_saige_step1_quant_with_sparse_GRM(step1_quant_input, sparse_grm_input)
    }

    /*
    Step 1 -> Step 2 Channel Emission Tuples
    Step 1 Out:  cohort, phenotype, rda, var
    Combine:     chromosome
    Step 2 In:   cohort, phenotype, rda, var, chromosome
    */
    step2_bin_input = step1_bin_output.combine(chromosome)
    
    step2_quant_input = step1_quant_output.combine(chromosome)

    //****************************************Conditional clauses for STEP 2 ***********************************************
    
    //***********************PLINK********************
    if(ftype.equals("PLINK")){
      //pulling out cohort, etc.
      genome_plink_file_bin = step2_bin_input.map {
        cohort, pheno, rda, var, chr -> \
        new Tuple(*plink_suffixes_list.collect { ext -> params.genome_plink_prefix + ext })
      } 
      //step2_bin_input.view{"Step 2 Bin Input: ${it}"}
      //genome_plink_file_bin.view{"Genome Plink files: ${it}"}
      step2_bin_output = call_saige_step2_PLINK_binary(step2_bin_input, genome_plink_file_bin)
    
      genome_plink_file_quant = step2_quant_input.map {
        cohort, pheno, rda, var, chr -> \
        new Tuple(*plink_suffixes_list.collect { ext -> params.genome_plink_prefix + ext })
      }
      step2_quant_output = call_saige_step2_PLINK_quant(step2_quant_input, genome_plink_file_quant)
    }

    //***********************BGEN*********************
    else if(ftype.equals("BGEN")){
        genome_bgen_file_bin = step2_bin_input.map {
            cohort, pheno, rda, var, chr -> \
            new Tuple(*bgen_suffixes_list.collect { ext -> params.genome_bgen_prefix + "chr" + chr +  ext })
        }
       //genome_bgen_file_bin.view{"Genome bin: ${it}"}
       //step2_bin_input.view{"Step 2.2 input: ${it}"}
       step2_bin_output = call_saige_step2_BGEN_binary(step2_bin_input,genome_bgen_file_bin)

       genome_bgen_file_quant = step2_quant_input.map {
        cohort, pheno, rda, var, chr -> \
          new Tuple(*bgen_suffixes_list.collect { ext -> params.genome_bgen_prefix + "chr" + chr + ext })
      }
      step2_quant_output = call_saige_step2_BGEN_quant(step2_quant_input,genome_bgen_file_quant)
    }

    //********************OOPS****************************
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
    step2_grouped_output = step2_all_output.groupTuple(by: [0, 1], size: params.chromosome_list.size())   
    merge_singles_script = "${launchDir}/scripts/merge_and_filter_saige_results.py"
    step2_grouped_output.view{"GO: ${it}"}
    (singles_merge_output, filtered_singles_output) = merge_and_filter_saige_gwas_output(step2_grouped_output, merge_singles_script)
    
    filtered_singles_output_list = filtered_singles_output.collect() //make list of filter paths
    // ANNOTATIONS 
    if (params['annotate']) {
            //combining all the filtered scripts and renaming/grabbing columns that are relevant into one file
            biofilter_input = gwas_make_biofilter_positions_input(filtered_singles_output_list)
            //adding nickname
            bf_input_channel = Channel.of('gwas_annotate').combine(biofilter_input)
            //calling biofilter
            biofilter_annots = BIOFILTER_POSITIONS(bf_input_channel)
            plotting_script = "${launchDir}/scripts/make_saige_gwas_plots_annotate.py"
            annot_singles = singles_merge_output.combine(biofilter_annots)
            annot_singles.view{"Annot Singles: ${it}"}
            plots = make_gwas_plots_with_annot(annot_singles, plotting_script)
            make_summary_table_with_annot(filtered_singles_output_list, biofilter_annots)
        }
        else {
           gwas_plots_script = "${launchDir}/scripts/make_saige_gwas_plots.py"
           gwas_plots = make_saige_gwas_plots(singles_merge_output, gwas_plots_script, pheno_table)
    
           make_summary_table(filtered_sumstats_list)
        }    

    
}