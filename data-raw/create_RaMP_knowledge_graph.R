library(Matrix)
library(tidyverse)

pathway_overlap_matrix_long<-
    readRDS("pathway_similarity/RaMP_pathway_overlap_matrix.Rds")
colnames(pathway_overlap_matrix_long)<-
    c("rownames","colnames","overlap")
chemical_similarity_matrix_long<-
    readRDS("chemical_similarity/RaMP_chemical_similarity_matrix.Rds")
colnames(chemical_similarity_matrix_long)<-
    c("rownames","similarity","colnames")

consensus_matrix_long<-pathway_overlap_matrix_long %>%
    left_join(chemical_similarity_matrix_long,
              by=c("rownames","colnames")) %>%
    replace_na(list(overlap=0,similarity=0)) %>%
    mutate(value=overlap+similarity) %>%
    select(-overlap,-similarity) %>%
    distinct(rownames,colnames, .keep_all=TRUE)

just_pathways=FALSE
if(just_pathways){
    consensus_matrix_long<-pathway_overlap_matrix_long %>%
        rename(value=overlap) %>%
        distinct(rownames,colnames, .keep_all=TRUE)
}


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
diag(consensus_matrix)=2
if(just_pathways){
    diag(consensus_matrix)=1
}

isSymmetric(consensus_matrix)
## Matrix is not symmetric at this stage, must induce symmetry

## for(i in 1:nrow(consensus_matrix)){
##     for(j in i:ncol(consensus_matrix)){
##         if(consensus_matrix[i,j]!=consensus_matrix[j,i]){
##             print(c(i,j))
##             break
##             if(consensus_matrix[i,j]==0){
##                 ##consensus_matrix[i,j]<-consensus_matrix[j,i]
##             }else if(consensus_matrix[j,i]==0){
##                 ##consensus_matrix[j,i]<-consensus_matrix[i,j]
##             }else{
##                 print(paste0("Conflicting information found at ",
##                              i,", ",j))
##             }
##         }
##     }
## }

## isSymmetric(consensus_matrix)

## disconnected_nodes<-which(apply(consensus_matrix,1,sum)==2)
## consensus_matrix<-consensus_matrix[-disconnected_nodes,
##                                      -disconnected_nodes]

saveRDS(as.matrix(consensus_matrix),"RaMP_consensus_matrix.Rds")

