#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("R/query_row.R")




query_builder <- function(id, label = NULL){
    ns <- NS(id)
    
   
    fluidRow(
        
            column(2,
                 actionButton(inputId = ns("add_term"), label = " + Search Term", style = "top : 172.5px") ,
                 actionButton(inputId = ns("rem_term"), label = " - Search Term", style = "margin-top : 0px; top : 247.5px")
                ),
            

        tags$div(class = "box_child",
            column(12, 
                               
                    tags$div(id = ns("new_row")),
                    row_query_ui_no_op("row0")
               
           
            
            )
        )
    )
   
}

query_builder_server <- function(id){
    moduleServer(id,
                 function(input, output, session){
                     ns <- session$ns
                     
                     row_count <- reactiveVal(0)       
                     
                     observeEvent(input$rem_term, {
                         if(row_count() >= 1 ){
                         newValue <- row_count() - 1     
                         row_count(newValue)}
                         #print(row_count)
                     })
                     
                     observeEvent(input$add_term, {
                         if(row_count() < 3){
                         newValue <- row_count() + 1     
                         row_count(newValue)}
                         #print(row_count)
                         
                     })
                     
                     inserted <- c()
                     
                     observeEvent(input$add_term, {
                         if(length(inserted) < 3){
                             id <- paste0("row",row_count())
                             
                             insertUI(selector = paste0("#", ns("new_row")),
                                       ui = 
                                    tags$div(id = id,
                                          
                                      row_query_ui(ns(id))
                                )
                             )
                             inserted <<- c(id, inserted)
                             #print(inserted)
                         }
                     })
                     
                     observeEvent(input$rem_term, {
                         if(length(inserted) > 0 ){
                             
                             removeUI(selector = paste0("#", inserted[length(inserted)])
                                     
                             )
                             inserted <<- inserted[-length(inserted)]
                             #print(inserted)
                         }
                     })
                     
                     
                     ##TODO:
                     #add server logic for receiving outputs from inner module
                     
                   
                     
                     

                     
                 })
}


# Define UI for application that draws a histogram
ui <- fluidPage(
    query_builder("q_build")
    

    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    query_builder_server("q_build")
    
}

# Run the application 
shinyApp(ui = ui, server = server)
