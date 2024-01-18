library(targets)
library(tarchetypes)
library(moiraine)

## add any package you need to use in the pipeline here
tar_option_set(
  packages = c(
    "moiraine",
    "MOFA2",
    "mixOmics",
    "readr",
    "dplyr",
    "ggplot2",
    "patchwork"
  )
)

## List of targets
list(

  ##==============##
  ## Data loading ----
  ##==============##

  ## Data import using a target factory
  import_dataset_csv_factory(
    files = c(
      system.file("extdata/genomics_dataset.csv", package = "moiraine"),
      system.file("extdata/transcriptomics_dataset.csv", package = "moiraine"),
      system.file("extdata/metabolomics_dataset.csv", package = "moiraine")
    ),
    col_ids = c("marker", "gene_id", "sample_id"),
    features_as_rowss = c(TRUE, TRUE, FALSE),
    target_name_suffixes = c("geno", "transcripto", "metabo")
  ),

  ## Features metadata import
  tar_target(
    fmetadata_file_geno,
    system.file("extdata/genomics_features_info.csv", package = "moiraine"),
    format = "file"
  ),

  tar_target(
    fmetadata_geno,
    import_fmetadata_csv(
      fmetadata_file_geno,
      col_id = "marker",
      col_types = c("chromosome" = "c")
    )
  ),

  import_fmetadata_csv_factory(
    files = c(
      system.file("extdata/metabolomics_features_info.csv", package = "moiraine")
    ),
    col_ids = c("feature_id"),
    target_name_suffixes = c("metabo")
  ),

  import_fmetadata_gff_factory(
    files = system.file("extdata/bos_taurus_gene_model.gff3", package = "moiraine"),
    feature_types = "genes",
    add_fieldss = c("Name", "description"),
    target_name_suffixes = "transcripto"
  ),

  ## Samples metadata import
  import_smetadata_csv_factory(
    files = system.file("extdata/samples_info.csv", package = "moiraine"),
    col_ids = "animal_id",
    target_name_suffixes = "all"
  ),

  ## Creating omics sets for each dataset
  create_omics_set_factory(
    datasets = c(data_geno, data_transcripto, data_metabo),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    features_metadatas = c(fmetadata_geno, fmetadata_transcripto, fmetadata_metabo),
    samples_metadatas = c(smetadata_all, smetadata_all, smetadata_all)
  ),

  ## Creating the MultiDataSet object
  tar_target(
    mo_set,
    create_multiomics_set(
      list(set_geno,
           set_transcripto,
           set_metabo)
    )
  ),

  ## Example with names
  tar_target(
    mo_set_with_names,
    create_multiomics_set(
      list(set_geno,
           set_transcripto,
           set_metabo),
      datasets_names = c("CaptureSeq", "RNAseq", "LCMS")
    )

  ),
  ## RNAseq differential expression results file
  tar_target(
    rnaseq_de_res_file,
    system.file(
      "extdata/transcriptomics_de_results.csv",
      package = "moiraine"
    ),
    format = "file"
  ),

  ## Reading the RNAseq differential expression results
  tar_target(
    rnaseq_de_res_df,
    read_csv(rnaseq_de_res_file) |>
      rename(feature_id = gene_id) |>
      mutate(dataset = "rnaseq")
  ),

  ## Adding the differential expression results to the MultiDataSet object
  tar_target(
    mo_set_de,
    add_features_metadata(mo_set, rnaseq_de_res_df)
  ),

  ## Applying transformations to the datasets
  transformation_datasets_factory(
    mo_set_de,
    c("rnaseq" = "vst-deseq2",
      "metabolome" = "best-normalize-manual"),
    methods = c("metabolome" = "log_x"),
    a = 0.01,
    b = 2,
    standardize = FALSE,
    transformed_data_name = "mo_set_transformed"
  ),

  ## Running a PCA on each dataset
  pca_complete_data_factory(
    mo_set_transformed,
    complete_data_name = "mo_set_complete"
  )

)
