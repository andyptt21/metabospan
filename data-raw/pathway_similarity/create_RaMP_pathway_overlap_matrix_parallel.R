library(tidyverse)
library(parallel)
library(RaMP)

pkg.globals <- setConnectionToRaMP(
  dbname = "ramp2", username = "root", conpass = "mysql123",
  host = "localhost", socket = paste0(
    "/lscratch/",
    Sys.getenv("SLURM_JOB_ID"),
    "/mysql/mysql.sock"
  )
)

detectBatchCPUs <- function() {
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  if (is.na(ncores)) {
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
  }
  if (is.na(ncores)) {
    return(4) # for helix
  }
  return(ncores)
}

ncpus <- detectBatchCPUs()
options(mc.cores = ncpus)

query <- "SELECT analytehasontology.rampOntologyId, analytehasontology.rampCompoundId, ontology.commonName FROM analytehasontology  LEFT JOIN ontology on analytehasontology.rampOntologyId = ontology.rampOntologyId;"

con <- connectToRaMP()

ontology_key <- DBI::dbGetQuery(con, query)
DBI::dbDisconnect(con)

generate_long_overlap_matrix <- function(biospecimen = NULL, filename,
                                         old_ramp = FALSE) {
  if (!is.null(biospecimen)) {
    biospec_metabs <- ontology_key %>%
      dplyr::filter(commonName == biospecimen) %>%
      pull(rampCompoundId)
    if (old_ramp) {
      metabHasPathway <- readRDS("metabHasPathway_old.Rds")
    } else {
      metabHasPathway <- readRDS("metabHasPathway.Rds")
    }
    metabHasPathway_filtered <- metabHasPathway %>%
      dplyr::filter(rampId %in% biospec_metabs)
    pathway_list <- split(
      metabHasPathway_filtered$pathwayRampId,
      metabHasPathway_filtered$rampId
    )
  } else {
    pathway_list <- split(
      metabHasPathway$pathwayRampId,
      metabHasPathway$rampId
    )
  }

  pathway_overlap_list <- mclapply(1:length(pathway_list), function(i) {
    overlap_v <- c()
    colnames_v <- c()
    for (j in 1:length(pathway_list)) {
      if (length(pathway_list[[i]]) > 1 & length(pathway_list[[j]]) > 1 &
        length(intersect(pathway_list[[i]], pathway_list[[j]])) != 0) {
        overlap <-
          length(
            intersect(
              pathway_list[[i]],
              pathway_list[[j]]
            )
          ) /
            min(
              length(pathway_list[[i]]),
              length(pathway_list[[j]])
            )
        colnames <- names(pathway_list)[j]
        overlap_v <- c(overlap_v, overlap)
        colnames_v <- c(colnames_v, colnames)
      }
    }
    if (i %% 10 == 0) {
      print(i)
    }
    return(data.frame(colnames_v, overlap_v))
  }, mc.cores = ncpus)

  names(pathway_overlap_list) <- names(pathway_list)

  pathway_overlap_df <- bind_rows(pathway_overlap_list,
    .id = "rownames"
  )

  pathway_overlap_df_check <- pivot_wider(
    pathway_overlap_df,
    names_from = c(colnames_v),
    values_from = overlap_v,
    values_fill = 0
  ) %>%
    column_to_rownames("rownames") %>%
    as.matrix()

  pathway_overlap_df_check <-
    pathway_overlap_df_check[colnames(
      pathway_overlap_df_check
    ), ]
  diag(pathway_overlap_df_check) <- 1

  isSymmetric(pathway_overlap_df_check)


  saveRDS(pathway_overlap_df, filename)
}

generate_long_overlap_matrix(NULL, "RaMP_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Blood", "Blood_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix(
  "Adipose tissue",
  "Adipose_pathway_overlap_matrix.Rds"
)
generate_long_overlap_matrix("Heart", "Heart_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Urine", "Urine_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Brain", "Brain_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Liver", "Liver_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Kidney", "Kidney_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Saliva", "Saliva_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix("Feces", "Feces_pathway_overlap_matrix.Rds")
generate_long_overlap_matrix(NULL, "KEGG_pathway_overlap_matrix.Rds", TRUE)
