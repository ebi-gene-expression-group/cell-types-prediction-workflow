#!/usr/bin/env nextflow 

// read query data from 10x directory into SCE object 
QUERY_10X_DIR = Channel.fromPath(params.query_10x_dir)
process read_query_sce {

    conda "${baseDir}/envs/dropletutils.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 5
    memory { 16.GB * task.attempt }

    input:
        file(query_10x_dir) from QUERY_10X_DIR

    output:
        file("query_sce.rds") into QUERY_SCE

    """
    dropletutils-read-10x-counts.R\
        --samples ${query_10x_dir}\
        --col-names ${params.col_names}\
        --output-object-file query_sce.rds
    """
}

process scpred_preprocess {
    // get scpred-specific matrices from query SCE object
    conda "${baseDir}/envs/scpred.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 5
    memory { 16.GB * task.attempt }

    input:
        file(query_sce) from QUERY_SCE

    output:
        file("query_expr_mat.rds") into SCPRED_QUERY_MAT

    """
    scpred_preprocess_data.R\
        --input-sce-object ${query_sce}\
        --normalised-counts-slot ${params.norm_counts_slot}\
        --output-matrix-object query_expr_mat.rds
    """
}

// load pre-trained scpred classifiers 
SCPRED_MODELS = Channel.fromPath(params.scpred_models).map{ f -> tuple("${f.simpleName}", f) }
process scpred_run {
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/scpred.yaml" 

    input:
        file(query_mat) from SCPRED_QUERY_MAT.first()
        set val(acc), file(model) from SCPRED_MODELS

    output:
        // expect classifier name to correspond to training data set accession code
        set val(acc), file("${acc}_predicted.csv") into SCPRED_OUTPUT

    """
    scpred_predict.R\
        --input-object ${model}\
        --pred-data ${query_mat}\
        --threshold-level ${params.pred_threshold}\
        --output-path ${acc}'_predicted.csv'
    """
}

process scpred_get_labels {
    conda "${baseDir}/envs/scpred.yaml" 

    input:
        set val(acc), file(scpred_output_tbl) from SCPRED_OUTPUT

    output:
        file("${acc}_final_labs.tsv") into SCPRED_LABELS

    """
    scpred_get_std_output.R\
            --predictions-file ${scpred_output_tbl}\
            --get-scores\
            --output-table ${acc}_labs.tsv

    # add metadata fields to output table
    echo "# dataset ${acc}" > ${acc}_final_labs.tsv
    echo "# tool scpred" >> ${acc}_final_labs.tsv
    cat ${acc}_labs.tsv >> ${acc}_final_labs.tsv
    """
}

// combine output files into single directory 
process scpred_combine_labels {
    conda "${baseDir}/envs/scpred.yaml" 

    input:
        file(labels) from SCPRED_LABELS.collect()
    output:
        file('scpred_labs') into SCPRED_LABELS_DIR

    """
    mkdir -p scpred_labs
    for file in ${labels}
    do
        mv \$file scpred_labs
    done
    """
}

process select_top_labs {
    conda "${baseDir}/envs/cell_types_analysis.yaml" 
    publishDir "${params.results_dir}", mode: 'copy'

    input:
        file(labels_dir) from SCPRED_LABELS_DIR
    output:
        file("scpred_combined_output.txt") into SCPRED_TOP_LABS

    """
    combine_tool_outputs.R\
        --input-dir ${labels_dir}\
        --top-labels-num 2\
        --scores\
        --output-table scpred_combined_output.txt
    """
}

