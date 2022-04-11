library(Matrix)
library(tidyverse)
library(RaMP)

pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="mysql123",
                                   host = "localhost",socket = paste0("/lscratch/",
                       Sys.getenv("SLURM_JOB_ID"),
                       "/mysql/mysql.sock"))

########################
## Pathways
########################

convert_to_similarity_matrix<-function(rds_file){
    long_matrix<-readRDS(rds_file)
    colnames(long_matrix) <-
    c("rownames","colnames","overlap")
    consensus_matrix_long<-long_matrix %>%
        rename(value=overlap) %>%
        distinct(rownames,colnames, .keep_all=TRUE)

    consensus_matrix<-pivot_wider(
        consensus_matrix_long,
        names_from=c(colnames),
        values_from=value,
        values_fill = 0)  %>%
        column_to_rownames("rownames") %>%
        as.matrix()

    consensus_matrix =
        consensus_matrix[colnames(
            consensus_matrix),]

    diag(consensus_matrix)=1
    
    if(isSymmetric(consensus_matrix)){
        return(consensus_matrix)
    }else{
        stop("Output is not symmetrical")
    }
}

RaMP_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/RaMP_pathway_overlap_matrix.Rds")
Blood_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Blood_pathway_overlap_matrix.Rds")
Adipose_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Adipose_pathway_overlap_matrix.Rds")
Heart_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Heart_pathway_overlap_matrix.Rds")
Urine_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Urine_pathway_overlap_matrix.Rds")
Brain_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Brain_pathway_overlap_matrix.Rds")
Liver_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Liver_pathway_overlap_matrix.Rds")
Kidney_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Kidney_pathway_overlap_matrix.Rds")
Saliva_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Saliva_pathway_overlap_matrix.Rds")
Feces_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/Feces_pathway_overlap_matrix.Rds")
KEGG_pathway_similarity_matrix<-convert_to_similarity_matrix("pathway_similarity/KEGG_pathway_overlap_matrix.Rds")

########################
## Structures
########################
full_chemical_similarity_matrix <- readRDS("chemical_similarity/RaMP_chemical_similarity_matrix.Rds")

convert_to_context_similarity_matrix<-function(ontology){
    full_matrix<-full_chemical_similarity_matrix

    query = paste0("SELECT analytehasontology.*, ontology.* from analytehasontology, ontology where ontology.commonName in ('",
                   ontology,
                   "') and ontology.rampOntologyId = analytehasontology.rampOntologyId")
    con <- connectToRaMP()
    ontology_metabs <- DBI::dbGetQuery(con,query)
    
    DBI::dbDisconnect(con)

    long_matrix<-full_matrix %>%
        dplyr::filter(ref_id %in% ontology_metabs$rampCompoundId &
                      ids_keep %in% ontology_metabs$rampCompoundId) %>%
        dplyr::relocate(ids_keep,ref_id,weights_keep)
    
    colnames(long_matrix) <-
        c("rownames","colnames","overlap")
    consensus_matrix_long<-long_matrix %>%
        rename(value=overlap) %>%
        distinct(rownames,colnames, .keep_all=TRUE)

    consensus_matrix<-pivot_wider(
        consensus_matrix_long,
        names_from=c(colnames),
        values_from=value,
        values_fill = 0)  %>%
        column_to_rownames("rownames")
    ## browser()
    dropped = setdiff(colnames(consensus_matrix),rownames(consensus_matrix))

    dropped_rows <- data.frame(matrix(0,ncol = ncol(consensus_matrix),
                                  nrow = length(dropped)))
    rownames(dropped_rows) <- dropped
    colnames(dropped_rows) <-colnames(consensus_matrix)

    consensus_matrix<-rbind(consensus_matrix,dropped_rows) %>% as.matrix
    
    consensus_matrix =
        consensus_matrix[colnames(
            consensus_matrix),]

    diag(consensus_matrix)=1
    
    if(isSymmetric(consensus_matrix)){
        return(consensus_matrix)
    }else{
        stop("Output is not symmetrical")
    }
}

Blood_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Blood")
Adipose_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Adipose tissue")
Heart_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Heart")
Urine_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Urine")
Brain_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Brain")
Liver_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Liver")
Kidney_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Kidney")
Saliva_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Saliva")
Feces_chemical_similarity_matrix<-convert_to_context_similarity_matrix("Feces")

########################
## Save
########################

usethis::use_data(RaMP_pathway_similarity_matrix,
                  Blood_pathway_similarity_matrix,
                  Adipose_pathway_similarity_matrix,
                  Heart_pathway_similarity_matrix,
                  Urine_pathway_similarity_matrix,
                  Brain_pathway_similarity_matrix,
                  Liver_pathway_similarity_matrix,
                  Kidney_pathway_similarity_matrix,
                  Saliva_pathway_similarity_matrix,
                  Feces_pathway_similarity_matrix,
                  KEGG_pathway_similarity_matrix,
                  Blood_chemical_similarity_matrix,
                  Adipose_chemical_similarity_matrix,
                  Heart_chemical_similarity_matrix,
                  Urine_chemical_similarity_matrix,
                  Brain_chemical_similarity_matrix,
                  Liver_chemical_similarity_matrix,
                  Kidney_chemical_similarity_matrix,
                  Saliva_chemical_similarity_matrix,
                  Feces_chemical_similarity_matrix,
                  overwrite = TRUE,
                  internal = TRUE)
