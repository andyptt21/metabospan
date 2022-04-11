##' Extract proximal metabolites from RaMP knowledge graph
##'
##' Uses output from "measure_metabolite_distance" to assemble vector of proximal metabolites to the set of interest
##' 
##' @param metabolites A list of database identifiers associated with metabolites of interest. Database of origin must be prepended to identifier (e.g. hmdb:HMDB0000001)
##' @param crosstalkers_full_set The output from "measure_metabolite_distance": a vector of distances between each node in RaMP and the seed set
##' @param cutoff Distance cutoff for extracted metabolites. "visualize_metabolite_distances" is helpful for determining a cutoff
##' @param percentile_cutoff Extracts top n percentile metabolites by lowest distance
##' @param absolute_cutoff Extracts top n metabolites by lowest distance
##' @param pathwaysOrChemicals Test distance on a pathway similarity or chemical similarity matrix. Choices are "pathways" or "chemicals"
##' @param biospecimen If NULL, test all metabolites in RaMP. Else, test networks built for specific biospecimens. Choices are "Blood", "Adipose", "Heart", "Urine", "Brain", "Liver", "Kidney", "Saliva", and "Feces"
##' @param keep_seeds Boolean for whether or not seed distances from input set should be calculated (almost always FALSE)
##' @return Vector of metabolites proximal to the seed set that met specified cutoff criterion
##' @importFrom magrittr %>%
##' @export
##' @author Andrew Christopher Patt
extract_proximal_metabolites<-function(metabolites,
                                       crosstalkers_full_set,
                                       cutoff=NULL,
                                       percentile_cutoff=NULL,
                                       absolute_cutoff=NULL,
                                       biospecimen = NULL,
                                       pathwaysOrChemicals = "pathways",
                                       keep_seeds=FALSE){
    target_ids<-rampFindSourceRampId(metabolites) %>%
        dplyr::select("rampId")
    knowledge_matrix<-load_knowledge_matrix(biospecimen,pathwaysOrChemicals)
    seeds<-which(colnames(knowledge_matrix) %in% as.matrix(target_ids))
    seeds=ifelse(colnames(knowledge_matrix) %in% as.matrix(target_ids),
                 1,0)
    if(!keep_seeds){
        db_metabolites<-rownames(knowledge_matrix)[-which(seeds==1)]
    }
    ##crosstalkers_trimmed_set_scores<-crosstalkers_full_set[which(crosstalkers_full_set>cutoff)]
    if(!is.null(cutoff)){
        crosstalkers_trimmed_set<-db_metabolites[which(crosstalkers_full_set>cutoff)]
    }else if(!is.null(absolute_cutoff)){
        crosstalkers_trimmed_set<-
            names(crosstalkers_full_set)[order(crosstalkers_full_set,decreasing=TRUE)[1:absolute_cutoff]]
    }else if(!is.null(percentile_cutoff)){
        percentile_cutoff_threshold<-stats::quantile(crosstalkers_full_set,percentile_cutoff)
        crosstalkers_trimmed_set<-
            db_metabolites[which(crosstalkers_full_set>percentile_cutoff_threshold)]
        
    }
    crosstalkers_trimmed_set<-data.frame(rampId=crosstalkers_trimmed_set)
    if(nrow(crosstalkers_trimmed_set)>1){
        crosstalkers_trimmed_set<-
            rampFindSourceFromId(
                       crosstalkers_trimmed_set %>% as.matrix() %>% as.vector()) ## %>%
        ## group_by("rampId") %>%
        ## filter(row_number() == 1) %>% 
        ## ungroup() %>%
        crosstalkers_trimmed_set <-
            crosstalkers_trimmed_set[!duplicated(crosstalkers_trimmed_set$rampId),] %>%
            dplyr::select("sourceId")
    }
    print(paste0("Extracting ",nrow(crosstalkers_trimmed_set)," metabolite(s)"))
    return(crosstalkers_trimmed_set %>% as.matrix() %>% as.vector())    
}
