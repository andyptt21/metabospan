metabHasOntology<-readRDS("metabHasOntology.Rds")

ontology_list<-split(metabHasOntology$rampOntologyIdLocation,
                    metabHasOntology$rampCompoundId)
overlap=c()
row=c()
column=c()
rownames=c()
colnames=c()
for(i in 751:1350){
    for(j in i:length(ontology_list)){
        if(length(ontology_list[[i]])>1 & length(ontology_list[[j]])>1 & length(intersect(ontology_list[[i]],ontology_list[[j]]))!=0){
            overlap=c(overlap,length(intersect(ontology_list[[i]],ontology_list[[j]]))/min(length(ontology_list[[i]]),
                                                                                         length(ontology_list[[j]])))
            ## row=c(row,i)
            rownames=c(rownames,names(ontology_list)[i])
            ## column=c(column,j)
            colnames=c(colnames,names(ontology_list)[j])
        }
    }
    if(i%%10==0){
        print(i)
    }
}
## ontology_overlap_matrix<-sparseMatrix(i=row,j=column,x=overlap,
##                                      dimnames=list(rownames,colnames),symmetric=TRUE)
ontology_overlap_df<-data.frame(rownames,colnames,overlap)
saveRDS(ontology_overlap_df,"RaMP_ontology_overlap_matrix_5.Rds")
