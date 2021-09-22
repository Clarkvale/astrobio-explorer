#' Shiny dashboardbox for a gprofiler2 result, gene set expression analysis 
#' 


library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(shinycssloaders)
library(gprofiler2)
library(plotly)


#'Helper function for formatting mappings between term ids and genes associated 
#'with them.
.format_intersection <- function(intersection_string){
  return(as.numeric(strsplit(intersection_string, split = ",")[[1]]))
}


manhat_ui  <- function(id, label = NULL, color = NULL, study_name = NULL){
  ns <- NS(id)
  div(tags$style(HTML(paste0( "#", id, "-", study_name, "_box .box-header.with-border", "{ background-color:", color, " !important;}"))),
  #enabling draggable & resizable box using. 
  #Need to specify header as a handle by using the class id
    jqui_draggable(options = list(handle = ".box-header.with-border",  containment = "#plot_space"), 
                   jqui_resizable(options =  list(minWidth = 300, maxWidth = 1000, maxHeight = 461, minHeight = 461 ),
                     box(id = ns(paste0(study_name, "_box")),
                         #TODO: Add variable title
                         #paste0("Gene Set Functional Groups ", "(", "GSEBLALA", " ",  gp$meta$query_metadata$organism, ")")
                         title = paste0(study_name, " Gene Functional Enrichment"),
                         width = 12,
                         collapsible = F,
                         closable = F,
                         plotlyOutput(ns("plot")) %>% withSpinner(type = 4)
                     )
                   )
    )
  
  )
  
  
}

manhat_server <- function(id, gprofiler, label = NULL){
  moduleServer(id, 
               function(input, output, session){
                 
                 ns <- session$ns
                 gp <- gprofiler
                 
                 output$plot <- renderPlotly({
                   g <- gostplot(gprofiler)
                   g$x$source <- id
                   
                   event_register(g, "plotly_click")
                   event_register(g, "plotly_selected")
                   g
                   
                 })
                 
                
                 output <- reactiveVal()
                 
                 # observeEvent(event_data(event = "plotly_click", source = id), 
                 #              ignoreInit = T, ignoreNULL = T,
                 #              { 
                 #                click  <- event_data(event = "plotly_click", source = id)$key
                 #                req(click)
                 #                  
                 #                out <-  .format_intersection(gp$result$intersection[which(gp$result$term_id %in% click)])
                 #                output(out)
                 #                  
                 #                print(output)}
                 #              )
                 #              
                 # 
                 # 
                 # observeEvent(event_data(event = "plotly_selected", source = id), 
                 #              ignoreInit = T, ignoreNULL = T,
                 #              { 
                 #                
                 #                sel  <- event_data(event = "plotly_selected", source = id)$key
                 #                req(sel)
                 #                  
                 #                out <-  .format_intersection(gp$result$intersection[which(gp$result$term_id %in% sel)])
                 #                output(out)
                 #                
                 #                
                 #                print(output)}
                 # )
                 # 
                 # 
                 # return(output)
                 
               })
}

ui <- dashboardPage( title = "manhat-test",
                     #wrapping everything in a dashBoard page makes the box work                
                     dashboardHeader(),
                     dashboardSidebar(),
                     dashboardBody(
                       tags$style("body { background-color: ghostwhite}"),
                       manhat_ui("plot")
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
  
  
  manhat_server("plot", gp, data)
}

# Run the application 
shinyApp(ui = ui, server = server)
