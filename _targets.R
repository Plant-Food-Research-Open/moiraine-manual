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
    "tibble",
    "tidyr",
    "dplyr",
    "ggplot2",
    "patchwork"
  )
)

test_multilevel <- function(smeta, multilevel) {
  mo_set <- c("A", "B") |>
    rlang::set_names(c("phenomics", "metabolomics")) |>
    purrr::map(
      ~ matrix(
        data = rnorm(nrow(smeta) * 30),
        ncol = nrow(smeta),
        dimnames = list(paste0(.x, "feature", 1:30), rownames(smeta))
      )
    ) |>
    purrr::imap(
      ~ create_omics_set(.x, .y, samples_metadata = smeta)
    ) |>
    create_multiomics_set()

  res <- get_input_spls(mo_set, mode = "canonical", multilevel = multilevel)

  attr(res, "multilevel")
}

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

  ##===================================##
  ## Adding information about features ----
  ##===================================##

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

  ##=========================##
  ## Datasets transformation ----
  ##=========================##

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


  ##============================================##
  ## Individual PCA and missing data imputation ----
  ##============================================##

  ## Running a PCA on each dataset
  pca_complete_data_factory(
    mo_set_transformed,
    complete_data_name = "mo_set_complete"
  ),


  ##===============##
  ## Pre-filtering ----
  ##===============##

  ## Unsupervised feature selection based on MAD score
  feature_preselection_mad_factory(
    mo_set_complete,
    to_keep_ns = c("snps" = 1000, "rnaseq" = 1000),
    with_ties = TRUE,
    filtered_set_target_name = "mo_presel_unsupervised"
  ),

  ## Supervised feature selection based on disease status
  feature_preselection_splsda_factory(
    mo_set_complete,
    group = "status",
    to_keep_ns = c("snps" = 1000, "rnaseq" = 1000),
    filtered_set_target_name = "mo_presel_supervised"
  ),

  ##=================##
  ## DIABLO pipeline ----
  ##=================##

  ## Creating the DIABLO input
  tar_target(
    diablo_input,
    get_input_mixomics_supervised(
      mo_presel_supervised,
      group = "status"
    )
  ),

  ## Running sPLS on each dataset to construct the design matrix
  diablo_pairwise_pls_factory(diablo_input),

  ## Initial DIABLO run with no feature selection and large number of components
  tar_target(
    diablo_novarsel,
    diablo_run(
      diablo_input,
      diablo_design_matrix,
      ncomp = 7
    )
  ),

  ## Cross-validation for number of components
  tar_target(
    diablo_perf_res,
    mixOmics::perf(
      diablo_novarsel,
      validation = "Mfold",
      folds = 10,
      nrepeat = 10,
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of components)
  tar_target(
    diablo_perf_plot,
    diablo_plot_perf(diablo_perf_res)
  ),

  ## Selected value for ncomp
  tar_target(
    diablo_optim_ncomp,
    diablo_get_optim_ncomp(diablo_perf_res)
  ),

  ## Cross-validation for number of features to retain
  tar_target(
    diablo_tune_res,
    diablo_tune(
      diablo_input,
      diablo_design_matrix,
      ncomp = diablo_optim_ncomp,
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      dist = "centroids.dist",
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of features)
  tar_target(
    diablo_tune_plot,
    diablo_plot_tune(diablo_tune_res)
  ),

  ## Final DIABLO run
  tar_target(
    diablo_final_run,
    diablo_run(
      diablo_input,
      diablo_design_matrix,
      ncomp = diablo_optim_ncomp,
      keepX = diablo_tune_res$choice.keepX
    )
  ),

  ##===============##
  ## sPLS pipeline ----
  ##===============##

  ## Showcasing multilevel with one factor
  tar_target(
    spls_smeta1,
    tibble(
      sample_id = paste0("sample_", 1:10),
      id = sample_id,
      plant_id = paste0("plant_", rep(1:5, each = 2)),
      treatment = rep(LETTERS[1:2], 5)
    ) |>
      column_to_rownames("sample_id") |>
      as.data.frame()
  ),

  tar_target(
    spls_multilevel1,
    test_multilevel(spls_smeta1, multilevel = "plant_id")
  ),

  ## Showcasing multilevel with two factors
  tar_target(
    spls_smeta2,
    expand_grid(
      plant_id = paste0("plant_", 1:2),
      treatment = LETTERS[1:2],
      time = 1:3
    ) |>
      mutate(
        sample_id = paste0("sample_", 1:n()),
        id = sample_id
      ) |>
      relocate(id) |>
      column_to_rownames("sample_id") |>
      as.data.frame()
  ),

  tar_target(
    spls_multilevel2,
    test_multilevel(spls_smeta2, multilevel = c("plant_id", "treatment", "time"))
  ),

  ## Creating sPLS input
  tar_target(
    spls_input,
    get_input_spls(
      mo_presel_supervised,
      mode = "canonical",
      datasets = c("rnaseq", "metabolome")
    )
  ),

  ## Initial PLS run with no feature selection and large number of components
  tar_target(
    spls_novarsel,
    spls_run(
      spls_input,
      ncomp = 4
    )
  ),

  ## Cross-validation for number of components
  tar_target(
    spls_perf_res,
    mixOmics::perf(
      spls_novarsel,
      validation = "Mfold",
      folds = 10,
      nrepeat = 10,
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of components)
  ## Can try criterion = 'Q2.total', 'cor.tpred', 'cor.upred', 'RSS.tpred',
  ## 'RSS.upred' (but avoid 'RSS' and 'PRESS')
  tar_target(
    spls_perf_plot,
    plot(spls_perf_res, criterion = "Q2.total")
  ),

  ## Selected value for ncomp
  tar_target(
    spls_optim_ncomp,
    spls_get_optim_ncomp(spls_perf_res, min_ncomp = 2)
  ),

  ## Cross-validation for number of features to retain
  tar_target(
    spls_tune_res,
    spls_tune(
      spls_input,
      ncomp = spls_optim_ncomp,
      keepX = seq(10, 100, 10),
      keepY = seq(10, 100, 10),
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      measure = "cor",
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of features)
  tar_target(
    spls_tune_plot,
    spls_plot_tune(spls_tune_res)
  ),

  ## Final sPLS run
  tar_target(
    spls_final_run,
    spls_run(
      spls_input,
      ncomp = spls_optim_ncomp,
      keepX = spls_tune_res$choice.keepX,
      keepY = spls_tune_res$choice.keepY
    )
  ),

  ##=================##
  ## sO2PLS pipeline ----
  ##=================##

  ## Creating sO2PLS input
  tar_target(
    omicspls_input,
    get_input_omicspls(
      mo_presel_supervised,
      datasets = c("rnaseq", "metabolome")
    )
  ),

  ## Adjusted cross-validation for number of components
  tar_target(
    so2pls_cv_adj,
    so2pls_crossval_o2m_adjR2(
      omicspls_input,
      a = 1:5,
      ax = seq(0, 10, by = 2),
      ay = seq(0, 10, by = 2),
      nr_folds = 10,
      nr_cores = 6,
      seed = 127
    )
  ),
  tar_target(
    so2pls_cv_adj_res,
    so2pls_get_optim_ncomp_adj(so2pls_cv_adj)
  ),

  ## Plotting adjusted cross-validation results
  tar_target(
    so2pls_cv_adj_plot,
    so2pls_plot_cv_adj(so2pls_cv_adj)
  ),

  ## Standard cross-validation for number of components
  tar_target(
    so2pls_cv,
    so2pls_crossval_o2m(
      omicspls_input,
      so2pls_cv_adj,
      nr_folds = 10,
      nr_cores = 6,
      seed = 356
    )
  ),
  tar_target(
    so2pls_cv_res,
    so2pls_get_optim_ncomp(so2pls_cv)
  ),

  ## Plotting standard cross-validation results
  tar_target(
    so2pls_cv_plot,
    so2pls_plot_cv(so2pls_cv)
  ),

  ## Cross-validation for sparsity parameters
  tar_target(
    so2pls_cv_sparsity,
    so2pls_crossval_sparsity(
      omicspls_input,
      n = so2pls_cv_res["n"],
      nx = so2pls_cv_res["nx"],
      ny = so2pls_cv_res["ny"],
      nr_folds = 10,
      keepx_seq = c(seq(5, 30, 5), seq(40, 100, 10)),
      keepy_seq = c(seq(5, 40, 5))
    )
  ),
  tar_target(
    so2pls_cv_sparsity_res,
    so2pls_get_optim_keep(so2pls_cv_sparsity)
  ),

  ## Plotting the results of the cross-validation for the number of features
  ## to retain from each dataset for the different joint components
  tar_target(
    so2pls_cv_sparsity_plot,
    so2pls_plot_cv_sparsity(so2pls_cv_sparsity)
  ),

  ## Extracting sparsity results in table format
  tar_target(
    so2pls_cv_sparsity_table,
    so2pls_print_cv_sparsity(so2pls_cv_sparsity_res)
  ),

  ## Final sO2PLS run
  tar_target(
    so2pls_final_run,
    so2pls_o2m(
      omicspls_input,
      so2pls_cv_res,
      so2pls_cv_sparsity_res
    )
  ),

  ##===============##
  ## MOFA pipeline ----
  ##===============##

  ## Creating MOFA input
  tar_target(
    mofa_input,
    get_input_mofa(
      mo_presel_supervised,
      options_list = list(
        data_options = list(scale_views = TRUE),
        model_options = list( likelihoods = c(
          "snps" = "poisson",
          "rnaseq" = "gaussian",
          "metabolome" = "gaussian"
        )
        ),
        training_options = list(seed = 43)
      ),
      only_common_samples = FALSE
    )
  ),

  ## Training MOFA model
  tar_target(
    mofa_trained,
    run_mofa(
      mofa_input,
      save_data = TRUE,
      use_basilisk = TRUE
    )
  ),

  ##========================##
  ## Results interpretation ----
  ##========================##

  ## Formatting outputs
  tar_target(
    spls_output,
    get_output(spls_final_run)
  ),

  tar_target(
    so2pls_output,
    get_output(so2pls_final_run)
  ),

  tar_target(
    mofa_output,
    get_output(mofa_trained)
  ),

  tar_target(
    diablo_output,
    get_output(diablo_final_run)
  ),

  ## Formatting output - individual latent dimensions
  tar_target(
    spls_output_no_average,
    get_output(spls_final_run, use_average_dimensions = FALSE)
  ),

  tar_target(
    so2pls_output_no_average,
    get_output(so2pls_final_run, use_average_dimensions = FALSE)
  ),

  tar_target(
    diablo_output_no_average,
    get_output(diablo_final_run, use_average_dimensions = FALSE)
  ),

  ##====================##
  ## Results evaluation ----
  ##====================##

  ## Creating features sets from features metadata
  tar_target(
    sets_single_omics,
    make_feature_sets_from_fm(
      mo_set_complete,
      col_names = list(
        "snps" = "qtl_type",
        "rnaseq" = "de_status",
        "metabolome" = "de_status"
      )
    )
  ),

  tar_target(
    sets_single_omics_merged,
    make_feature_sets_from_fm(
      mo_set_complete,
      col_names = list(
        "snps" = "qtl_type",
        "rnaseq" = "de_status",
        "metabolome" = "de_status"
      ),
      combine_omics_sets = TRUE
    )
  ),

  ## Reading GO annotation file
  tar_target(
    rnaseq_go_terms_file,
    system.file(
      "extdata/transcriptomics_go_annotation.csv",
      package = "moiraine"
    ),
    format = "file"
  ),

  tar_target(
    rnaseq_go_df,
    read_csv(rnaseq_go_terms_file) |>
      filter(go_domain == "Biological process")
  ),

  ## Making GO terms sets
  tar_target(
    go_sets,
    make_feature_sets_from_df(
      rnaseq_go_df,
      col_id = "gene_id",
      col_set = "go_id"
    )
  ),

  ## Filtering GO term sets against measured features
  tar_target(
    go_sets_filtered,
    reduce_feature_sets_data(go_sets, mo_set_complete)
  ),

  ## Checking genes GO term sets against datasets
  tar_target(
    go_sets_check,
    check_feature_sets(
      go_sets_filtered,
      mo_set_complete,
      datasets = "rnaseq"
    )
  ),

  ## Table of information about GO terms
  tar_target(
    go_sets_info,
    rnaseq_go_df |>
      dplyr::select(go_id, go_name) |>
      dplyr::distinct()
  ),

  ## DIABLO latent components enrichment analysis
  tar_target(
    diablo_enrichment_results,
    evaluate_method_enrichment(
      diablo_output,
      go_sets_filtered,
      datasets = "rnaseq",
      use_abs = TRUE,
      min_set_size = 10,
      add_missing_features = TRUE,
      mo_data = mo_set_complete,
      sets_info_df = go_sets_info,
      col_set = "go_id"
    )
  ),

  tar_target(
    mofa_enrichment_results,
    evaluate_method_enrichment(
      mofa_output,
      go_sets_filtered,
      datasets = "rnaseq",
      latent_dimensions = paste("Factor", 1:3),
      use_abs = TRUE,
      min_set_size = 10,
      add_missing_features = TRUE,
      mo_data = mo_set_complete,
      sets_info_df = go_sets_info,
      col_set = "go_id"
    )
  ),

  ## Assessing samples clustering
  tar_target(
    diablo_silhouette,
    compute_samples_silhouette(
      diablo_output,
      mo_set_complete,
      "status"
    )
  ),

  tar_target(
    mofa_silhouette,
    compute_samples_silhouette(
      mofa_output,
      mo_set_complete,
      "status",
      latent_dimensions = paste("Factor", 1:3)
    )
  ),

  ##====================##
  ## Results comparison ----
  ##====================##

  ## Creating the input object for the MOFA pipeline
  ## using the unsupervised preselection results
  tar_target(
    mofa_unsupervised_input,
    get_input_mofa(
      mo_presel_unsupervised,
      options_list = list(
        data_options = list(scale_views = TRUE),
        model_options = list(
          likelihoods = c(
            "snps" = "poisson",
            "rnaseq" = "gaussian",
            "metabolome" = "gaussian")
        ),
        training_options = list(seed = 72)
      ),
      only_common_samples = FALSE
    )
  ),

  ## Training the model with the MOFA algorithm
  tar_target(
    mofa_unsupervised_trained,
    run_mofa(
      mofa_unsupervised_input,
      save_data = TRUE,
      use_basilisk = TRUE
    )
  ),

  ## Formatting MOFA output
  tar_target(
    mofa_unsupervised_output,
    get_output(mofa_unsupervised_trained)
  ),

  ## List of formatted output
  tar_target(
    output_list,
    list(spls_output, so2pls_output, mofa_output, diablo_output)
  ),

  tar_target(
    output_list_mofa,
    list(
      "MOFA (supervised pref.)" = mofa_output,
      "MOFA (unsupervised pref.)" = mofa_unsupervised_output
    )
  )
)
