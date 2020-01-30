#!/usr/bin/env nextflow 


// specify query data channels 
QUERY_MAT = Channel.fromPath(params.query_mat)
QUERY_BARCODES = Channel.fromPath(params.query_barcodes)
QUERY_GENES = Channel.fromPath(params.query_genes)
QUERY_MARKERS = Channel.fromPath(params.query_markers)

// download query data 
process download_data{

  input:
    file(query_mat) from QUERY_MAT
    file(query_barcodes) from QUERY_BARCODES
    file(query_genes) from QUERY_GENES

  output:
    file("${params.query_10x_dir}") into QUERY_10X_DIR


  """
  get_query_data.R 

  """

}

// read query data from 10x directory into SCE object 
process read_query_sce {
  conda "${baseDir}/envs/dropletutils.yaml"

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


SCPRED_MODELS = Channel.fromPath(params.scpred_models_dir)
process scpred_run {
  conda "${baseDir}/envs/scpred.yaml" 
  // run scpred prediction on each classifier emitted from MODELS channel
  // need pre-processed counts matrix built from query data 
  input:
    file(query_mat) from SCPRED_QUERY_MAT
    file(model) from SCPRED_MODELS

  output:
    set val(), file() into 


  """
  name=$(echo a.txt | cut -d . -f 1)

  scpred_predict.R\
          --input-object ${eval_trained_model}\
          --pred-data ${query_mat}\
          --output-path ${}
  """



}

process scpred_get_labels {}

SCMAP_CELL_MODELS = Channel.fromPath(params.scmap_cell_models_dir)
process run_scmap_cell {}

process get_scmap_cell_labels {}

SCMAP_CLUST_MODELS = Channel.fromPath(params.scmap_clust_models_dir)
process run_scmap_clust {}

process get_scmap_clust_labels {}











