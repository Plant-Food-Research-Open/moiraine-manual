library(targets)
library(tarchetypes)
library(moiraine)

## add any package you need to use in the pipeline here
tar_option_set(
  packages = c(
    "moiraine"
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
  )

)
