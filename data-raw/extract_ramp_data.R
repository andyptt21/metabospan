library(RMySQL)
library(tidyverse)
library(RaMP)

pkg.globals <- setConnectionToRaMP(
  dbname = "ramp2", username = "root", conpass = "mysql123",
  host = "localhost", socket = paste0(
    "/lscratch/",
    Sys.getenv("SLURM_JOB_ID"),
    "/mysql/mysql.sock"
  )
)

con <- connectToRaMP()

query <- "select * from analytehaspathway"

metabHasPathway <- DBI::dbGetQuery(con, query) %>%
  dplyr::filter(grepl("C", rampId))
DBI::dbDisconnect(con)

acceptable_size_pathways <- names(which(table(metabHasPathway$pathwayRampId) < 200))

metabHasPathway <- metabHasPathway %>%
  dplyr::filter(pathwayRampId %in% acceptable_size_pathways)

saveRDS(metabHasPathway, "pathway_similarity/metabHasPathway.Rds")

pkg.globals <- setConnectionToRaMP(
  dbname = "old_ramp", username = "root", conpass = "mysql123",
  host = "localhost", socket = paste0(
    "/lscratch/",
    Sys.getenv("SLURM_JOB_ID"),
    "/mysql/mysql.sock"
  )
)

con <- connectToRaMP()

query <- "select * from analytehaspathway"

metabHasPathway <- DBI::dbGetQuery(con, query) %>%
  dplyr::filter(grepl("C", rampId)) %>%
  dplyr::filter(pathwaySource == "kegg")
DBI::dbDisconnect(con)

acceptable_size_pathways <- names(which(table(metabHasPathway$pathwayRampId) < 200))

metabHasPathway <- metabHasPathway %>%
  dplyr::filter(pathwayRampId %in% acceptable_size_pathways)

saveRDS(metabHasPathway, "pathway_similarity/metabHasPathway_old.Rds")



## Chemical Similarity

## pathway_list<-split(metabHasPathway$pathwayRampId,
##                     metabHasPathway$rampId)

## metabs_of_interest<-names(pathway_list)[which(lapply(pathway_list,length) > 2)]

pkg.globals <- setConnectionToRaMP(
  dbname = "ramp2", username = "root", conpass = "mysql123",
  host = "localhost", socket = paste0(
    "/lscratch/",
    Sys.getenv("SLURM_JOB_ID"),
    "/mysql/mysql.sock"
  )
)

query <- "select * from chem_props"
con <- connectToRaMP()

chem_props <- DBI::dbGetQuery(con, query) ##  %>%
## dplyr::filter(ramp_id %in% metabs_of_interest)
DBI::dbDisconnect(con)

repeats <- chem_props$ramp_id[which(
  duplicated(chem_props$ramp_id)
)]

chem_props_duplicated <- chem_props %>%
  dplyr::filter(ramp_id %in% repeats)

`%notin%` <- Negate(`%in%`)

chem_props <- chem_props %>%
  dplyr::filter(ramp_id %notin% repeats)

for (i in unique(chem_props_duplicated$ramp_id)) {
  temp <- chem_props_duplicated %>%
    dplyr::filter(ramp_id == i)
  if ("hmdb" %in% temp$chem_data_source) {
    chem_props <- rbind(
      chem_props,
      temp %>%
        dplyr::filter(chem_data_source == "hmdb") %>%
        slice(1)
    )
  } else {
    chem_props <- rbind(
      chem_props,
      temp %>%
        dplyr::filter(chem_data_source == "chebi") %>%
        slice(1)
    )
  }
}

saveRDS(chem_props, "chemical_similarity/metabHasStructure.Rds")
