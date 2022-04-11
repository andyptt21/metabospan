##' Visualize the distance distribution between nodes in the RaMP network and metabolites of interest
##'
##' Generate a histogram of distance measurements from each node in the RaMP knowledge network to the metabolite set of interest (uses the output from measure_metabolite_distance)
##'
##' @param crosstalkers_full_set The output from "measure_metabolite_distance": a vector of distances between each node in RaMP and the seed set
##' @return Histogram of distance distribution (plot)
##' @export
##' @author Andrew Christopher Patt
visualize_metabolite_distances <- function(crosstalkers_full_set) {
  graphics::hist(crosstalkers_full_set,
    breaks = 200,
    main = paste0(
      "Metabolite proximity distribution"
    )
  )
}


##' Internal function for visualizing metabolite network
##'
##' Shows seed metabolites, extracted metabolites, and their connectivity in the RaMP knowledge graph
##'
##' @param metabolites A list of database identifiers associated with metabolites of interest. Database of origin must be prepended to identifier (e.g. hmdb:HMDB0000001)
##' @param crosstalkers Vector of metabolites returned by "extract_proximal_metabolites"
##' @param pathwaysOrChemicals Test distance on a pathway similarity or chemical similarity matrix. Choices are "pathways" or "chemicals"
##' @param biospecimen If NULL, test all metabolites in RaMP. Else, test networks built for specific biospecimens. Choices are "Blood", "Adipose", "Heart", "Urine", "Brain", "Liver", "Kidney", "Saliva", and "Feces"
##' @param label Boolean specifying if nodes in visual should be labeled with metabolite common names
##' @param filterIslands Boolean to display "Islands" (nodes with no direct connections to other nodes in visual)
##' @return Visual showing connectivity of seeds and extracted nodes from data
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
##' @author Andrew Christopher Patt
visualize_metabolite_network <- function(metabolites, crosstalkers,
                                         biospecimen = NULL,
                                         pathwaysOrChemicals = "pathways",
                                         label = TRUE, filterIslands = TRUE) {
  browser()
  knowledge_matrix <- load_knowledge_matrix(biospecimen, pathwaysOrChemicals)
  metabolites <- rampFindSourceRampId(metabolites) %>%
      dplyr::pull("rampId") %>%
      unique()
  crosstalkers <- rampFindSourceRampId(crosstalkers) %>%
      dplyr::pull("rampId") %>%
      unique()
  metabolites <- metabolites[!is.na(metabolites)]
  crosstalkers <- crosstalkers[!is.na(crosstalkers)]
  metabolites <- intersect(
    metabolites,
    colnames(knowledge_matrix)
  )
  metabolites_source <- rampFindSourceFromId(
    metabolites
  ) %>%
    dplyr::group_by(.data$rampId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select("sourceId")
  full_graph <-
      igraph::graph.adjacency(knowledge_matrix, mode = "undirected",
                              weighted  = TRUE, diag = F)
  ## Find shortest paths between all metabolites in set of interest + crosstalkers
  combined_metabolites <- c(metabolites,crosstalkers)
  path_metabolites <- lapply(1:length(combined_metabolites), function(x) {
    lapply(x:length(combined_metabolites), function(y) {
      if (x == y) {
        return(NULL)
      } else {
          shortestPath <- igraph::shortest_paths(full_graph,combined_metabolites[x],
                                                 combined_metabolites[y])
        return(shortestPath$vpath[[1]])
      }
    })
  })

  path_metabolites <- names(unlist(path_metabolites)) %>% unique %>%
    setdiff(c(metabolites, crosstalkers))

  crosstalkers_matrix <- knowledge_matrix[
    c(metabolites, crosstalkers, path_metabolites),
    c(metabolites, crosstalkers, path_metabolites)
  ]
  igraph <- igraph::graph.adjacency(crosstalkers_matrix,
    mode = "undirected", weighted = TRUE, diag = F
  )
  ## V(igraph)$color<-cluster_louvain(igraph)$membership
  igraph <- igraph %>%
      igraph::set_vertex_attr("size", value =
                                          c(rep(4,times = length(c(metabolites,crosstalkers))),
                                            rep(1,times = length(path_metabolites))
                              ))
  labels <- c(source_to_common(metabolites_source) %>%
    as.matrix() %>% as.vector(), rep(
    "",
    length(c(crosstalkers, path_metabolites))
  ))
  igraph <- igraph %>%
    igraph::set_vertex_attr("label", value = labels)

  custom_pal <- RColorBrewer::brewer.pal(12, "Set3")
  ## g<-network::as.network.matrix(as.data.frame(crosstalkers_matrix))
  ## g %v% "color" = igraph::cluster_louvain(igraph)$membership

  colors <- igraph::cluster_louvain(igraph)$membership
  ## g %v% "color" = c(c(g %v% "color")[1:length(metabolites)]*2,
  ##                   c(g %v% "color")[(length(metabolites)+1):ncol(crosstalkers_matrix)]*2-1)
  colors <- custom_pal[colors]
  colors[1:length(metabolites)] <-
    colorspace::darken(colors[1:length(metabolites)],
      amount = 0.3
    )
  colors[(length(c(metabolites, crosstalkers, path_metabolites)) - (length(path_metabolites))+1):length(c(metabolites, crosstalkers, path_metabolites))] <- "black"

  shapes <- c(rep("square", length(metabolites)), rep("circle", length(c(
    crosstalkers,
    path_metabolites
  ))))

  igraph <- igraph %>%
    igraph::set_vertex_attr("color", value = colors)
  ## V(igraph)$color <- g %v% "color"
  ## V(igraph)$color[1:length(metabolites)]=
  ##             colorspace::darken(V(igraph)$color[1:length(metabolites)],
  ##                                amount=0.3)

  igraph <- igraph %>%
    igraph::set_vertex_attr("label.color",
      value = colorspace::darken(colors,
        amount = 0.3
      )
    )
  ## V(igraph)$label.color=
  ##             colorspace::darken(V(igraph)$color,amount=0.3)
  ## V(igraph)$label.dist=1
  igraph <- igraph %>%
    igraph::set_vertex_attr("label.dist", value = 1)
  igraph <- igraph %>%
    igraph::set_vertex_attr("shape", value = shapes)
  igraph <- igraph %>%
    igraph::set_edge_attr("color",
      value = "gray80"
    )
  igraph <- igraph %>%
    igraph::set_vertex_attr("label.dist", value = 1)
  ## V(igraph)$shape="circle"
  ## V(igraph)$shape[1:length(metabolites)]=
  ## rep("square",length(metabolites))
  ## E(igraph)$color="gray80"
  if (filterIslands) {
    isolated <- which(igraph::degree(igraph) == 0)
    igraph <- igraph::delete.vertices(igraph, isolated)
  }
  if (label) {
    igraph::plot.igraph(igraph, margin = c(0, 0, 0, 0))
  } else {
    igraph::plot.igraph(igraph, margin = c(0, 0, 0, 0), vertex.label = NA)
  }
}

#' Find rampId from given source ID
#' The rampId can be plugged in other functions to continue query
#' @param sourceId a data frame or string separated by comma or string
#' separated by new line
#' @return data.frame that has sourceId and rampId and source as columns
rampFindSourceRampId <- function(sourceId) {
  if (is.character(sourceId)) {
    if (grepl("\n", sourceId)[1]) {
      list_metabolite <- strsplit(sourceId, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else if (grepl(",", sourceId)[1]) {
      list_metabolite <- strsplit(sourceId, ",")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- sourceId
    }
  } else if (is.data.frame(sourceId)) {
    list_metabolite <- as.character(unlist(sourceId$sourceId))
  } else {
    message("Wrong Format of argument")
    return(NULL)
  }
  list_metabolite <- sapply(list_metabolite, shQuote)
  list_metabolite <- paste(list_metabolite, collapse = ",")
  con <- RaMP::connectToRaMP()
  query <- paste0("select sourceId,IDtype as analytesource, rampId from source where sourceId in (", list_metabolite, ");")
  df <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)
  return(df)
}

#' Find all source from given list of RaMP Ids
#' @param rampId could be a data frame return by rampFindSynonymFromSynonym
#' containing all information related to synonym. Or can be a list of
#' rampIds
#' @param full return whole searching result or not (TRUE/FALSE)
#' @return a data frame that has all source Id in the column or the source table that has metaoblites entry
rampFindSourceFromId <- function(rampId = NULL, full = TRUE) {
  if (is.data.frame(rampId)) {
    list_id <- rampId$rampId
  } else if (is.character(rampId)) {
    if (grepl("\n", rampId)[1]) {
      list_id <- strsplit(rampId, "\n")
      list_id <- unlist(list_id)
    } else if (grepl(",", rampId)[1]) {
      list_id <- strsplit(rampId, ",")
      list_id <- unlist(list_id)
    } else {
      list_id <- rampId
    }
  } else {
    message("Wrong format of input")
    return(NULL)
  }
  list_id <- unique(list_id)
  list_id <- sapply(list_id, shQuote)
  list_id <- paste(list_id, collapse = ",")
  query <- paste0("select * from source where rampId in (", list_id, ");")

  con <- RaMP::connectToRaMP()
  df <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)
  if (full) {
    return(df)
  } else {
    return(df[, 1])
  }
}

#' Find rampId from given common name vector
#' @param source_ids a vector of metabolite common name strings
#' @return data.frame that has rampId as column
source_to_common <- function(source_ids) {
  query <- "select * from source"
  con <- RaMP::connectToRaMP()
  source_table <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

  source_ids <- data.frame(sourceId = source_ids)
  common_names <- source_ids %>%
    dplyr::left_join(source_table, by = "sourceId")
  common_names <- common_names[!duplicated(common_names$rampId), ] %>%
    dplyr::select("commonName")
  return(common_names)
}
