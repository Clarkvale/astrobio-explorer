#'volcano.R 
#'
#'Shiny Module for building an interactive plotly volcano plot in a shiny 
#'dashboard box.
library(plotly)
library(shiny)
library(manhattanly)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(shinycssloaders)
library(shinyBS)



#'Helper function for formatting mappings between term ids and genes associated 
#'with them.
.format_intersection <- function(intersection_string){
    return(as.numeric(strsplit(intersection_string, split = ",")[[1]]))
}

volcano_ui  <- function(id, label = NULL, color = NULL, study_name = NULL){
    ns <- NS(id)
    #browser()
    div(tags$style(HTML(paste0( "#", id,"_box", " .box-header.with-border", "{ background-color:", color, " !important;}"))),
    
    #enabling draggable & resizable box using. 
    #Need to specify header as a handle by using the class id
    jqui_draggable(options = list(handle = ".box-header.with-border", containment = "#plot_space"), 
      
        box(id = paste0(id, "_box"),
            title = paste0(study_name, " Gene Expression Volcano Plot"),
            width = 12,
            collapsible = F,
            closable = F,
            sidebar = boxSidebar(
                id = ns("myboxsidebar"),
                width = 40,
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
                                   "Wikipathway Annotations" = "WP",
                                   "All" = "all"),
                    selected = "all"),
                shinyBS::bsTooltip(ns("term_type"), "Subset genes by their inclusion in an enriched functional group category.")
                ),
            plotlyOutput(ns("plot")) %>% withSpinner(type = 4)
            )
        
    )
)  
    
   
    
}

volcano_server <- function(id, gprofiler, dataset,  subset = NULL, label = NULL){
    moduleServer(id, 
                 function(input, output, session){
                    
                     ns <- session$ns
                     
                     #if I ask for a parameter only inside a reactive context it seems like they won't update properly
                     gp <- gprofiler
                     d <- dataset
                     s <- subset
                     #browser()
                     
                     #getting variable gene ids
                     if(!any(grepl("Entrez\\.ID", colnames(dataset), ignore.case = T))){
                         if(any(grepl("Platform\\.ORF", colnames(dataset), ignore.case = T))){
                              gene_ids <- grep("Platform\\.ORF", colnames(dataset), value = T, ignore.case = T) 
                         } else {
                             gene_ids <- grep("Gene\\.symbol", colnames(dataset), ignore.case = T, value = T)
                         }
                     } else {
                         gene_ids <- grep("Entrez\\.ID", colnames(dataset), ignore.case = T, value = T)[[1]]
                     }
                     
                     #apparently this is a problem
                     dataset$adj.P.Val <- as.numeric(dataset$adj.P.Val)
                     dataset$logFC <- as.numeric(dataset$logFC)
                     
                     filtered_data <- reactive({
                                                if(input$term_type != "all"){
                         
                                                    genes <- dplyr::filter(gprofiler$result, 
                                                                       source == input$term_type)$intersection
                                                    genes <- unique(unlist(sapply(genes, .format_intersection, USE.NAMES = F)))
                                                
                                                    d_out <- dataset[which(dataset[[gene_ids]] %in% genes),] 
                                                    } else {
                                                        d_out <- dataset 
                                                    }
                                                  if(!is.null(subset) && length(subset) != 0){
                                                    d_out <-  d_out[which(dataset[[gene_ids]] %in% subset),] 
                                                  }
                       
                                                  d_out
                         
                                                })
                      
                    
                     
                     output$plot <- renderPlotly(
                         volcanoly(filtered_data(), gene = gene_ids, 
                                   effect_size = "logFC", 
                                   p = "adj.P.Val", 
                                   genomewideline = -log10(1e-2), effect_size_line = FALSE, 
                                   snp = "Gene.symbol", annotation1 = "adj.P.Val", title = NULL))
                    
                     
                     #placeholder return Value to make sure shiny knows its there when rendered
                     return(TRUE)
                    })
}

#testing shiny app
ui <- dashboardPage( title = "volcano-test",
                     #wrapping everything in a dashBoard page makes the box work                
                     dashboardHeader(),
                     dashboardSidebar(),
                     dashboardBody(
                         tags$style("body { background-color: ghostwhite}"),
                         volcano_ui("plot")
                     )
                     
)


server <- function(input, output) {
    #example data    
    data <- read.csv("test_data\\GSE4136_25thGen.csv")
    fdata <- data %>% dplyr::filter(adj.P.Val <= 0.1) %>% dplyr::filter(logFC >= 1 || logFC <= -1)
    
    up <- fdata %>% filter(logFC >= 1) %>% na.omit 
    down <- fdata %>% filter(logFC <= -1) %>% na.omit
    
    gp <- gprofiler2::gost(list("up-regulated" =  as.numeric(up$Entrez.ID), "down-regulated" = as.numeric(down$Entrez.ID)), 
               organism = "scerevisiae", numeric_ns = "ENTREZGENE_ACC"
               ,ordered_query = T, multi_query = F, evcodes = T)
    subset <- c("853158","855581","853043")
    
    
    
    volcano_server("plot", gp, data, subset = subset)
}

# Run the application 
shinyApp(ui = ui, server = server)