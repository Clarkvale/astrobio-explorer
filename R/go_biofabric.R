#'GO_Fabric_Builder
#'
#'Shiny module which creates a reactive biofabric graph from a gprofiler2 fgse 
#'results

library(shiny)
library(gprofiler2)
library(igraph)
library(RBioFabric)
library(dplyr)
library(htmlwidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(shinyBS)

validate_empty_set <- function(df, term_label){
    if(nrow(df) == 0){
        return(paste0("No enriched gene functional groups of type ", term_label, "."))
    }else{
        return(NULL)
    }
}

#'Builds adjacency list from gpprofiler2:: gost() output
build_adjacency_list <- function(gost_out, dataset, term_type, upordown){
    
    gost_out <- gost_out$result %>% dplyr::filter(source == term_type) %>% dplyr::filter(query == upordown)
    
    validate(validate_empty_set(gost_out, term_label =  term_type))
    intersections <- as.list((sapply(gost_out$intersection,FUN = entrez2name, 
                            dataset = dataset, USE.NAMES = T, simplify = F)))
    
    names(intersections) <- gost_out$term_name
    
    out <- c()
    for(i in 1:length(intersections)){
        a <- intersections[[i]]
        names(a) <- rep(names(intersections)[i], length(a))
        out <- append(out, a)
    }


    
    adj_list <- unlist(out)
    return(adj_list)
    
    
}
#' helper function which converts entrez ID to gene name using a limma derived dataset
entrez2name <- function(go_entry, dataset){
    return(unique(dataset$Gene.symbol[which(dataset$Entrez.ID %in% .format_intersection(go_entry))]))
}


format_term_id <- function(intersection, term){
    
    term_rep <- rep(term, length(intersection))
    names(intersection) <- term_rep
    return(intersection)
    
}

#'Helper function for formatting mappings between term ids and genes associated 
#'with them.
.format_intersection <- function(intersection_string){
    return(as.numeric(strsplit(intersection_string, split = ",")[[1]]))
}

go_fabric_ui  <- function(id, study_name = NULL, color = NULL, label = NULL){
    ns <- NS(id)
    
    #enabling draggable & resizable box using. 
   div(tags$style(HTML(paste0( "#", id, "-", study_name, "_box .box-header.with-border", "{ background-color:", color, " !important;}"))), #Need to specify header as a handle by using the class id
   jqui_draggable( options = list(handle = ".box-header.with-border", containment = "#plot_space"), 
        box(title = paste0(study_name, " BioFabric Gene Group Enrichment Plot"), 
        width = 12,
        id = ns(paste0(study_name, "_box")),
        collapsible = F,
        closable = F,
        sidebar = boxSidebar(
            id = ns("mycardsidebar"),
            width = 40,
            selectInput(
                ns("upordown"),
                label = "Fold Change Direction",
                choices = list("Up" = "up-regulated", "Down" = "down-regulated"),
                selected = "up-regulated"
                
            ),
            selectInput(
                ns("term_type"),
                label = "Enrichment Type",
                choices = list("Biological Process" = "GO:BP", 
                               "Cell Component" = "GO:CC",
                               "Molecular Function" = "GO:MF",
                               "Transcription Factor Motif" = "TF",
                               "KEGG" = "KEGG",
                               "Reactome"= "REAC",
                               "Human Phenotype Ontology" = "HP",
                               "miRNA Target Annotations" = "MIRNA",
                               "Protein Complex Annotations (CORUM)" = "CORUM",
                               "Wikipathway Annotations" = "WP"),
                selected = "GO:BP"
                
            ),
            shinyBS::bsTooltip(ns("term_type"), "Subset genes by their inclusion in an enriched functional group category."),
            checkboxInput(ns("shadowLink"), "Display Shadow Links"),
            shinyBS::bsTooltip(ns("shadowLink"), "Show duplicate edges.")
        ),
    
    
       bioFabric_htmlwidgetOutput(ns("fabric"))
    ))
    
   )
    
    
}

go_fabric_server <- function(id, gprofiler, dataset, subset = NULL, label = NULL){
    moduleServer(id, 
                 function(input, output, session){
                     
                   if(length(subset) != 0){
          
                     gprofiler$result <- gprofiler$result[subset,]
      
                     
                   }
                    else{
                     gp <- gprofiler
                     
                     }
                   d <- dataset
                     
                     
                     network <- reactive({ 
                                            adj_list <- build_adjacency_list(gprofiler, 
                                                               dataset, 
                                                               term_type  = input$term_type, 
                                                               upordown = input$upordown)
                     
                                         nodes <- data.frame( id = unique(unname(c(adj_list, names(adj_list)))))
                                         edges <- data.frame(from = names(adj_list), to = adj_list)
                                         igraph::graph_from_data_frame(d=edges, directed=F)})
                     
                     
                     output$fabric <- renderBioFabric_htmlwidget(
                         bioFabric_htmlwidget(bioFabric(network(), shadowLinks = input$shadowLink), 
                                              zoomMin = 0.5, zoomMax = 180))
                    
                     
                     
                 })
}

ui <- dashboardPage( title = "fabric-test",
     #wrapping everything in a dashBoard page makes the box work                
    dashboardHeader(),
    dashboardSidebar(),
        dashboardBody(
            tags$style("body { background-color: ghostwhite}"),
            go_fabric_ui("plot", study_name = "GSE123", color = "blue")
         )
        
)


server <- function(input, output) {
    #example data    
    data <- read.csv("C:\\Users\\Benja\\repos\\de-expression-igem\\datasets\\GSE4136_Scer\\GSE4136_25thGen.csv")
    fdata <- data %>% dplyr::filter(adj.P.Val <= 0.1) %>% dplyr::filter(logFC >= 1 || logFC <= -1)
    
    up <- fdata %>% filter(logFC >= 1) %>% na.omit 
    down <- fdata %>% filter(logFC <= -1) %>% na.omit
    
    gp <- gost(list("up-regulated" =  as.numeric(up$Entrez.ID), "down-regulated" = as.numeric(down$Entrez.ID)), 
               organism = "scerevisiae", numeric_ns = "ENTREZGENE_ACC"
               ,ordered_query = T, multi_query = F, evcodes = T)
    
    
    
    
    
    go_fabric_server("plot", gp, data)
}

# Run the application 
shinyApp(ui = ui, server = server)
