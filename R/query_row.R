#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


row_query_ui <- function(id, label = NULL){
    ns <- NS(id)
    
        tags$div(
                       
                       fluidRow(
                           column(4, style = "padding-right : 0px;",
                                       selectInput(inputId = ns("org_query"), label = NULL,
                                                   choices = c("E.coli"= "ecol_term" , "Homo sapiens" = "hsap_term")
                                                   
                                       )
                       ),
                       
                       column(4, style = "padding-left : 0px; padding-right: 0px;",
                              selectInput(
                                  inputId = ns("plat_query"), label = NULL,
                                  choices = c("RNAseq" = "rna_term", "Microarray" = "micro_term", "Either"= "either_term")
                                  
                              )
                              
                       ),
                       
                       column(4, style = "padding-left : 0px; padding-right: 0px;",
                              selectInput(
                                  inputId = ns("grav_query"), label = NULL,
                                  choices = c("Spaceflown" = "space_term", "HARV" = "harv_term", "RPM" = "rpm_term", "RCCS" = "rccs_term")
                              )
                       )),
                       
                       
                       
                      
                      fluidRow(
                        column(3, offset = 4,
                               selectInput(
                                   inputId = ns("op_query"), label = NULL,
                                   choices = c("AND" = "and_op", "OR" = "or_op"))
                               ))
            )
                       
    
}

row_query_ui_no_op <- function(id, label = NULL){
    ns <- NS(id)
    
    tags$div(
        
        
        fluidRow(
            column(4, style = "padding-right : 0px;",
                   selectInput(inputId = ns("org_query"), label = NULL,
                               choices = c("E.coli"= "ecol_term" , "Homo sapiens" = "hsap_term")
                               
                   )
            ),
            
            column(4, style = "padding-left : 0px; padding-right: 0px;",
                   selectInput(
                       inputId = ns("plat_query"), label = NULL,
                       choices = c("RNAseq" = "rna_term", "Microarray" = "micro_term", "Either"= "either_term")
                       
                   )
                   
            ),
            
            column(4, style = "padding-left : 0px; padding-right: 0px;",
                   selectInput(
                       inputId = ns("grav_query"), label = NULL,
                       choices = c("Spaceflown" = "space_term", "HARV" = "harv_term", "RPM" = "rpm_term", "RCCS" = "rccs_term")
                   )
            ))
    )
}



row_query_ui_server <- function(id){
    moduleServer(id, 
                 function(input, output, session){
                    return(list( org = reactive(input$org_query),
                                 plat = reactive(input$plat_query),
                                 grav = reactive(input$grav_query),
                                 op = reactive(input$op_query)
                        
                    ))
                     
                 })
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    row_query_ui("row1")

   
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    row_query_ui_server("row1")

    
}

# Run the application 
shinyApp(ui = ui, server = server)
