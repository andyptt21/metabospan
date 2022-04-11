##' Compute distance to seed set for each node in RaMP knowledge graph
##'
##' Uses Markov random walk with restarts algorithm to generate a by-node distance measurement from each node in the RaMP knowledge graph to each node associated with the input metabolites.
##'
##' @param metabolites A list of database identifiers associated with metabolites of interest. Database of origin must be prepended to identifier (e.g. hmdb:HMDB0000001)
##' @param pathwaysOrChemicals Test distance on a pathway similarity or chemical similarity matrix. Choices are "pathways" or "chemicals"
##' @param biospecimen If NULL, test all metabolites in RaMP. Else, test networks built for specific biospecimens. Choices are "Blood", "Adipose", "Heart", "Urine", "Brain", "Liver", "Kidney", "Saliva", and "Feces"
##' @param keep_seeds Boolean for whether or not seed distances from input set should be calculated (almost always FALSE)
##' @return Vector of distances from each node in RaMP knowledge network to seed set
##' @importFrom magrittr %>%
##' @export
##' @author Andrew Christopher Patt
measure_metabolite_distance <- function(metabolites,
                                        biospecimen = NULL,
                                        pathwaysOrChemicals = "pathways",
                                        keep_seeds = FALSE) {
  target_ids <- rampFindSourceRampId(metabolites) %>%
    dplyr::select("rampId")

  knowledge_matrix <- load_knowledge_matrix(biospecimen, pathwaysOrChemicals)

  seeds <- which(colnames(knowledge_matrix) %in% as.matrix(target_ids))
  if (length(seeds) < 3) {
    stop("Less than three input metabolites could be mapped to the specified network. Three or more must be mappable for valid MetaboSPAN results")
  }
  seeds <- ifelse(colnames(knowledge_matrix) %in% as.matrix(target_ids),
    1, 0
  )
  pt <- diffusr::random.walk(seeds, knowledge_matrix,
    do.analytical = TRUE, correct.for.hubs = TRUE
  )
  crosstalkers_full_set <- pt$p.inf
  names(crosstalkers_full_set) <- colnames(knowledge_matrix)
  if (!keep_seeds) {
    crosstalkers_full_set <- crosstalkers_full_set[-which(seeds == 1)]
  }
  return(crosstalkers_full_set)
}


load_knowledge_matrix <- function(biospecimen, pathOrChem) {
  if (pathOrChem == "pathways") {
    if (is.null(biospecimen)) {
      return(RaMP_pathway_similarity_matrix)
    } else if (biospecimen == "Blood") {
      return(Blood_pathway_similarity_matrix)
    } else if (biospecimen == "Adipose") {
      return(Adipose_pathway_similarity_matrix)
    } else if (biospecimen == "Heart") {
      return(Heart_pathway_similarity_matrix)
    } else if (biospecimen == "Urine") {
      return(Urine_pathway_similarity_matrix)
    } else if (biospecimen == "Brain") {
      return(Brain_pathway_similarity_matrix)
    } else if (biospecimen == "Liver") {
      return(Liver_pathway_similarity_matrix)
    } else if (biospecimen == "Kidney") {
      return(Kidney_pathway_similarity_matrix)
    } else if (biospecimen == "Saliva") {
      return(Saliva_pathway_similarity_matrix)
    } else if (biospecimen == "Feces") {
      return(Feces_pathway_similarity_matrix)
    } else if (biospecimen == "KEGG") {
      return(KEGG_pathway_similarity_matrix)
    } else {
      stop("Biospecimen not recognized")
    }
  } else if (pathOrChem == "chemicals") {
    if (biospecimen == "Blood") {
      return(Blood_chemical_similarity_matrix)
    } else if (biospecimen == "Adipose") {
      return(Adipose_chemical_similarity_matrix)
    } else if (biospecimen == "Heart") {
      return(Heart_chemical_similarity_matrix)
    } else if (biospecimen == "Urine") {
      return(Urine_chemical_similarity_matrix)
    } else if (biospecimen == "Brain") {
      return(Brain_chemical_similarity_matrix)
    } else if (biospecimen == "Liver") {
      return(Liver_chemical_similarity_matrix)
    } else if (biospecimen == "Kidney") {
      return(Kidney_chemical_similarity_matrix)
    } else if (biospecimen == "Saliva") {
      return(Saliva_chemical_similarity_matrix)
    } else if (biospecimen == "Feces") {
      return(Feces_chemical_similarity_matrix)
    } else {
      stop("Must choose a specific biospecimen for chemical analysis")
    }
  } else {
    stop("pathwaysOrChemicals must be 'pathways' or 'chemicals'")
  }
}
