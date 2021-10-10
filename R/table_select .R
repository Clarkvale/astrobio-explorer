#Shiny module for selecting gses via table


library(shiny)
library(DT)
library(RSQLite)
library(shinyFeedback)
library(shinyWidgets)
library(dplyr)
library(spsComps)
library(shinyjs)

size_validate <- function(df, in_p, in_fc){
  f_rows <- nrow(dplyr::filter(df, logFC >= in_fc) %>% dplyr::filter(adj.P.Val <= in_p))
  
  if(f_rows > 5000){
    return("Too many genes filtered for FGSE analysis. Try lowering the p-value or increasing the fold change requirement.")
  }
  else if(f_rows <= 2){
    return("Too few genes filtered for FGSE analysis. Try increasing the p-value or decreasing the fold change requirement.")
  }
  else{
    return(NULL)
  }
}

study_validate <- function(codes, df){
  if(length(codes) > 3){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

make_org_name <- function(string){
  prefix <- tolower(substr(strsplit(string, " ")[[1]][1], 1, 1))
  suffix <- tolower((strsplit(string, " ")[[1]][2]))
  return(paste0(prefix, suffix))
}

process_data <- function(dataset, GSEcode, Organism, pval_filter = 0.1, logfc_filter = 1, exclude = F, perform_gea = T){
  require(gprofiler2)
  #browser()
  if(perform_gea){
    filtered <- dataset %>% dplyr::filter(adj.P.Val <= pval_filter) %>% 
      dplyr::filter(logFC >= logfc_filter | logFC <= -logfc_filter)
    up <- filtered %>% dplyr::filter(logFC >= logfc_filter)  
    down <- filtered %>% dplyr::filter(logFC <= -logfc_filter) 
    
    shinyCatch( 
      gp <- gprofiler2::gost(list("up-regulated" =  as.numeric(up$Entrez.ID), 
                                    "down-regulated" = as.numeric(down$Entrez.ID)), 
                               organism = Organism, numeric_ns = "ENTREZGENE_ACC"
                               ,ordered_query = T, multi_query = F, evcodes = T, exclude_iea = exclude),
      blocking_level = "none"
    
  
    )
  }
  #here I'm making sure we got a result from gprofiler 
  gp <- tryCatch(
    expr = gp,
    
    error = function(e){
      e
      gp <- NULL
    }
    
  )
  
  outname <- paste0(GSEcode, ": ", Organism)
  return(list( "data" = list("gene_expression" = dataset,"gp_result" = gp), "name" = outname))
  
  
}



table_select_ui <- function(id, label = NULL){
  ns <- NS(id)
  
  
  
    fluidPage(
      useShinyjs(),
      useShinyFeedback(),
      tags$div(class = "box_parent", 
        fluidRow(
       
       
         column(10,DTOutput(ns("gseTable"))),
         
         column(2,textAreaInput(
           inputId = ns("gse_codes"), placeholder = "The resulting GSE codes will appear here. Each code links to the GEO page where the study was found. 
                                                                  Remove codes to omit from analysis.", 
           label = "Selected Studies", rows = 3, resize = "vertical"),
           tags$h5(tags$b("Links to Source Database")),
           htmlOutput(ns("links")),
           
           checkboxInput(ns("perform"), label = "Perform Gene Enrichment Analysis"),
           
           tags$div(class = "options",
            
            tags$h5(tags$b("Enrichment Analysis Options")),
            numericInput(inputId = ns("pval"), value = 0.05, label = "Filter values by P value (less than)", min = 0, max = 1, step = 0.01),
            numericInput(inputId = ns("logfc"), value = 1, label = "Filter values by log2 fold change (greater than)", min = 0, max = 10, step = 1),
            checkboxInput(ns("exclude"), label = "Exclude In silico GO Annotations (IEA)")
           ),
           actionBttn(inputId = ns("import"), "Import!", style = "fill", color = "primary")
           )
         )
      ),
     tags$br(),
     tags$div(class = "summary_box",
      
     fluidRow(column(2, offset = 6),
              htmlOutput(ns("title"))),
     fluidRow(
       
       
         column(5, tags$h5(tags$b("Experiment Summary:")), 
                  htmlOutput(ns("summary"))),
         
         column(7, tags$h5(tags$b("Differential Gene Expression Contrast Model:")),
                  uiOutput(ns("metaTable"))))
       )
    )
  
  
}


table_select_server <- function(id){
  moduleServer(id,
               function(input, output, session){
                 ns <- session$ns
                 
                 #defining output from module
                 codes <- reactiveValues(gses = NULL,
                                         btn = NULL,
                                         org = NULL)
                 
                 #connecting to metaDB, sqlite DB which stores study metadata
                 conn <- dbConnect(RSQLite::SQLite(), dbname = "AstroMeta.db")
                 gse_table <- as.data.frame(tbl(conn, "STUDIES"))
                 dbDisconnect(conn)
                 
                 

            
                 
                 
                 output$gseTable = renderDT(gse_table, options = list(lengthChange = FALSE,
                                                                      rowCallback = JS(sprintf('function(row, data) {
                                                                                               $(row).mouseenter(function(){
                                                                                                   var hover_index = $(this)[0]._DT_RowIndex
                                                                                                   /* console.log(hover_index); */
                                                                                                   Shiny.onInputChange("%s", hover_index);
                                                                                              });}', ns("hoverIndexJS"))
                                                                        ),
                                                                      pageLength = 15))
                 
                 #adding proxy to edit table on client side
                 proxy <-  dataTableProxy("gseTable")
                 
                 observe({
                   #waiting for a row input from the interactive table
                   selected_rows <- input$gseTable_rows_selected
                   #displays rows based on table selection
                   gses <- gse_table$acc[selected_rows]
                   #print(gses)
                   
                   
                   updateTextAreaInput(session, "gse_codes", 
                                       placeholder = "The resulting GSE codes will appear here. Each code links to the GEO page where the study was found. 
                                                                Remove codes to omit from analysis.", 
                                       value =  paste0(gses, sep = "\n", collapse = ""))
                   
                   
                   
                 })
                 
                 
                 #generating contrast tables on row hover
                 contrast_table_DB <- list.files(path = "datasets", pattern = "meta.*csv", recursive = T)
                 meta_info_DB <- list.files(path = "datasets", pattern = "meta.*txt", recursive = T)
                 observeEvent(input$hoverIndexJS,{
         
                   hovered_gse <- gse_table$acc[input$hoverIndexJS + 1] #JS is 0 indexed 
                   #browser()
                   query <- sapply(hovered_gse, paste0, "_meta", USE.NAMES = F)
                   meta_paths <- sapply(query, FUN = function(q){paste0("datasets/", grep(pattern = q, contrast_table_DB, value = T))}, USE.NAMES = F)
                   
                   tables <-  tagList(sapply(meta_paths, FUN = function(path){renderTable(read.csv(path))}))
                   
                   
                   #getting study summaries from .txt files
                   summary_path <- paste0("datasets/", grep(pattern =  paste0(hovered_gse, "_meta"), meta_info_DB, value = T))
                   meta.txt <-  readr::read_delim((summary_path), delim = "\n", show_col_types = F)
                   pos <- grep(meta.txt[[1]], pattern = "summary:")
                   summary_end <- T
                   summary <- c()
                   #browser()
                   while(summary_end){
                     pos <- pos + 1
                     if(grepl(meta.txt[[pos,1]], pattern = "supplementary_file:") | grepl(meta.txt[[pos,1]], pattern = "title:")){
                       summary_end <- F
                     }else{
                       summary <- append(summary, meta.txt[[pos,1]])
                     }
                   }
                   output$title <- renderUI(tags$h3(tags$b(paste0(hovered_gse, " Info"))))
                   output$summary <- renderUI(paste0(summary, collapse = "\n"))
                   output$metaTable <- renderUI(tables)
                   
                 })
                 
                 #generating links to GEO database
                 output$links <- renderUI({
                   
                   selected<- input$gseTable_rows_selected
                   req(selected)
                   
                   #displays rows based on table selection
                   gses <- gse_table$acc[selected]
                   pref <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
                   links <- sapply(gses, function(gse){return(paste0("<a href=",pref, gse,">", gse, "</a>"))}, USE.NAMES = F)
                   #print(links)
                   HTML(paste(links, sep = "<br/>"))
                 })
                 
                 #using observeEvent to only trigger the update when textAreaBox is changed by the user.
                 #This block also does the work of both validating gse codes and returning the codes to the rest of the app
                 observe({
                   
                   text_box_acc <- na.omit(as.vector(strsplit(input$gse_codes, "\n")))[[1]]
                   
                   
                   written_rows <- which(gse_table$acc %in% text_box_acc)
                   #print(written_rows)
                   
                   
                   #Showing feedback on textarea, showFeedBackSuccess() has some weird behavior so I need to wrap my conditions into a function
                   #TODO: Update with more conditions
                   hideFeedback("gse_codes")
                   req(written_rows != 0)
                 
                   if(!study_validate(written_rows, NULL)){
                     showFeedbackDanger(
                       inputId = "gse_codes", 
                       text ="Too many studies selected for analysis"
                     )
                    
                   }
                   
                   if(any(which(gse_table$GSE_Analysis_Compatible[written_rows] == "no"))){
                     shinyjs::disable(id = "perform")
                     updateCheckboxInput(inputId = "perform", value = FALSE)
                     showFeedbackWarning(inputId = "gse_codes", text = "Study not supported for functional gene set enrichment analysis.")
                   }
                   
                   else{
                     
                     shinyjs::enable(id = "perform")
                     showFeedbackSuccess("gse_codes")
                     
                     #selection the rows based off the textAreaInput
                     #proxy %>% selectRows(written_rows)
                     
                     org_name <- gse_table$org[written_rows]
                     
                     #finally storing the gse codes  and organisms
                     codes$gses <- text_box_acc
                     codes$org <- org_name
                     
                     
                     
                   }
                 })
                 
                 
                 # observeEvent(input$import,{
                 #   codes$btn <- input$import
                 # })
                 
                 
                 
                 #processing data with gprofiler2, this needs to be ..//datasets when running the module independently
                 myDB <- list.files(path = "datasets", pattern = "GSE\\d*\\.csv", recursive = T)
                 #Select studies page
                 
                
                 #listening for an actionbutton press in the select table module
                 paths <- reactiveVal()
                 observeEvent(input$import, {
                  
                   path_num <-  sapply(codes$gses, grep,  x = myDB)
                   paths(myDB[path_num])
                   
                 })
                 
                 
                 
                 
               
                
                 output_data <- reactiveVal()
                 pre_output <- list()
                 
                 observeEvent(input$import,  ignoreNULL = T, ignoreInit = T,{
                   
                   pre_output <- list() 
                  
                  tryCatch(
                    #browser(),
                    #this needs to be ..//datasets// when running the module independently
                    selected_data <- lapply(paste0("datasets//", paths()), read.csv),
                    
                    error = function(cond) {
                      showNotification(message(cond), type = "error")
                      message("Problem with accessing data: ")
                      message(cond)
                      validate(need(NULL, message = paste(cond)))
                      
                    })
                   
                   
                   #validating size
                  
                    for(i in 1:length(selected_data)){
                      erstring <- NULL
                      #browser()
                      if(input$perform){
                        erstring <- size_validate(selected_data[[i]], in_p = input$pval, in_fc = input$logfc)
                      }
                      
                      if(!is.null(erstring)){
                        showNotification(erstring, type = "error")
                        validate(need(NULL, message = erstring))
                        }
                      }
                    
                  #TODO: change message based on input
                  withProgress(message = "Performing Gene Functional Enrichment Analysis",value = 0, max = length(codes), 
                   {
                     
                     for(i in 1:length(selected_data)){
                       incProgress(amount = 1, detail = paste0("Working on ", codes$gses[i]))
                      
                         pre_output[[codes$gses[i]]] <- process_data(dataset =  selected_data[[i]],
                                                                              GSEcode = codes$gses[i],
                                                                              Organism =  make_org_name(codes$org[i]),
                                                                              pval_filter = input$pval,
                                                                              logfc_filter = input$logfc,
                                                                              exclude = input$exclude,
                                                                              perform_gea = input$perform)
                       
                         

                     }
                  
                      
                                
                  })
                  
                  
                  output_data(pre_output)
                  
                  if(any(sapply(output_data(), is.null))){
                    validate(need(NULL, message = "Error"))
                  }
                  else{
                    #browser() 
                    show_alert(title = "Complete!", type = "success")
                  }
                 })
                 
                 
                 observeEvent(input$perform, {
                   shinyjs::toggle(selector = ".options" )
                 })
                 
               
                output_data
                 
                 }
               )
  
          
}




#' TESTING SHINY APP 
#'______________________________________________________________________________



ui <- fluidPage(
  table_select_ui("table")
  
  
)

server <- function(input, output) {
  
  
  pdata <- table_select_server("table")
  
  observeEvent(pdata(), {print(pdata())})
  
 
}

# Run the application 
shinyApp(ui = ui, server = server)
