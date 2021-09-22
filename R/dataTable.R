#' Shiny dashboardbox for a DT datatable 
#' 


library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(gprofiler2)
library(dplyr)
library(DT)

dataTable_ui  <- function(id, label = NULL, study_name = NULL, color = NULL){
  ns <- NS(id)
  
  div(tags$style(HTML(paste0("#", id, "-", study_name, "_box .box-header.with-border", "{ background-color:", color, " !important; }",
                             "#", study_name, "_table-", study_name,"_box { width: 100%;}"))),
  
  
    box(id = ns(paste0(study_name, "_box")),
        #TODO: Add variable title
        #paste0("Gene Set Functional Groups ", "(", "GSEBLALA", " ",  gp$meta$query_metadata$organism, ")")
        title = paste0(study_name, " ", "Enriched Gene Functional Groups"),
        width = 12,
        collapsible = T,
        solidHeader = T,
  
        
        div(DT::DTOutput(ns("plot")), style = "font-size: 100%; width: inherit;")
      )
  )
  
  
  
  
  
  
  
}

dataTable_server <- function(id, gprofiler,  dataset, label = NULL){
  moduleServer(id, 
               function(input, output, session){
                 
                 ns <- session$ns
                 
                 slim_gp <- gprofiler$result %>% 
                   dplyr::select(query, term_id, source, term_name, intersection_size, p_value) %>% dplyr::rename(direction = query)
                 
                 slim_dt <- dataset %>% dplyr::select(Gene.symbol, Gene.title, 
                                                   Entrez.ID, logFC, P.Value, 
                                                   adj.P.Val) %>% 
                                                   filter(Entrez.ID != "")
                 output$plot <- renderDT(
                   datatable(slim_gp, options =  list(pageLength = 15,  
                                                      lengthMenu = c(5, 10, 15, 20),
                                                      autoWidth = TRUE,
                                                      columnDefs = list(
                                                          list(
                                                            targets = c(4))
                                                            # render = JS("function(data, type, row, meta) { 
                                                            #         return type === 'display' && data.length > 20 ?
                                                            #         '<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;
                                                            #         }")
                                                          
                                                          
                                                      )
                                                  )) %>% DT::formatStyle(names(slim_gp),lineHeight='70%') 
                                  )
                 # I dont think this actually does something anymore
                 reactive({
                      as.numeric(strsplit(gprofiler$result$intersection[which(gprofiler$result$term_id == slim_gp[input$plot_rows_selected,]$term_id)] , ",")[[1]])
            
                      })
               })
}

ui <- dashboardPage( title = "table-test",
                     #wrapping everything in a dashBoard page makes the box work                
                     dashboardHeader(),
                     dashboardSidebar(),
                     dashboardBody(
                       tags$style("body { background-color: ghostwhite}"),
                       dataTable_ui("plot", study_name = "GSE123", color = "aliceblue")
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
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
