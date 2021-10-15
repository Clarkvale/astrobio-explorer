#License

# Â© Copyright 2021 iGEM Concordia, Benjamin Clark
# This file is part of AstroBio Explorer.
# AstroBio Explorer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# AstroBio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with AstroBio Explorer.  If not, see <https://www.gnu.org/licenses/>.

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(gprofiler2)
library(shinyWidgets)
library(dplyr)
library(dbplyr)
library(shinyjs)
library(DT)
library(shinyjqui)
library(purrr)



#'building dynamic sidebar menu
#'
build_menu <- function(input_data){
    tags_out <- list()

    for(i in 1:length(names(input_data))){
        study <- names(input_data)[i]
        if(!is.null(input_data[[study]]$data$gp_result)){
          tags_out[[i]] <-  menuItem(text = study, icon = icon("server"),
                                             menuSubItem(tabName =  paste0(study,"_gexp_tab"), text = paste0("Gene Expression"), icon = icon("line-chart")),
                                             menuItem(text = paste0("Gene Functional Groups"), 
                                                      menuSubItem(tabName =  paste0(study, "_manhat_tab"), text  = "Manhattan Plot", icon = icon("line-chart")),
                                                      menuSubItem(tabName =  paste0(study, "_biofab_tab"), text = "BioFabric Plot", icon = icon("line-chart")),
                                                      menuSubItem(tabName =  paste0(study, "_results_tab"), text = "Gene Group Table", icon = icon("table")))
            )
        } else{
          tags_out[[i]] <-  menuItem(text = study, icon = icon("server"),
                                     menuSubItem(tabName =  paste0(study,"_gexp_tab"), text = paste0("Gene Expression"), icon = icon("line-chart"))
                                     
                                    
          )
          
        }
                                  
    }
    
    return(tags_out)
}

#' Creates an organism name suitable for gprofiler2 to use 
make_org_name <- function(string){
  prefix <- tolower(substr(strsplit(string)[[1]][1], 1, 1))
  suffix <- tolower((strsplit(string)[[1]][2]))
  return(paste0(prefix, suffix))
}

#' Formats selected genes into a string describing the groups size.
make_gene_readout_string <- function(gene_list){
  out_vec <- c()
  for(i in 1:length(gene_list)){
    out_vec <- append(out_vec, paste0(names(gene_list[i]), ": ", length(gene_list[[i]]), " genes selected"))
  }
  return(paste(out_vec, collapse =  " "))
}

#Converting a standard display name into a id
name2idPrefix <- function(name){
  return(paste0(strsplit(name, split = ":")[[1]][1]))
}






#UI DEFINITION



ui <- dashboardPage(title = "AstroBio Explorer",
    
    dashboardHeader(title = "AstroBio Explorer"),
    dashboardSidebar(collapsed = F,
      sidebarMenu(id = "import_tab", menuItem(text = "Import From Database", tabName = "import", 
                                          icon = icon("database"))),
      useShinyjs(),
        
        
        
          tags$br(),
          tags$h5(class = "sidebar_head", "Study Views"),
          sidebarMenu(id = "tabs_menu", menuItemOutput("menuItem")),
          tags$h5(class = "sidebar_head", "Gene Selections"), 
          textOutput("gene_select"),
          actionButton("gene_show", label = "Reset\\Show Genes on Plot"), 
          actionButton("reset", label = "Reset Selections"),
          tags$br(),
          downloadButton(outputId = "download", label = "Download Results")
        ),
       
      
    dashboardBody(
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
      useShinyjs(),
        tabItem(tabName = "import_tab", table_select_ui("import")),
        uiOutput(outputId = "tabs")
  )

   
)




#SERVER DEFINITION
server <- function(input, output) {
  
  
  #SETTING UP DB CONNECTION AND USER SELECTION PROMPT TAB
  
  input_data <- table_select_server("import")
 
    
    
    #DEFINING COLOR PALETTE AND ALL POSSIBLE PLOT IDS
  studies_palette <- reactiveVal()
  observeEvent(input_data(), ignoreInit = T,{
    colors <- list()
    colors<- suppressWarnings(RColorBrewer::brewer.pal(n = length(input_data()), name = "Pastel2"))
    names(colors) <- names(input_data())
    
    studies_palette(colors)
    
    
  
    #BUILDING SIDEBAR MENU
    output$menuItem <- renderMenu({sidebarMenu(id = "studies_tab", .list =  build_menu(input_data()))})
    
  })
    
    
    
    #rendering gene selection readout
    output$gene_select <- renderText({
      validate(need(input_data(), "No input data"))
      
      n <- sel_genes()
      make_gene_readout_string(n)
    })
    
    #datatable row selections
    sel_row <- reactive({
      
      #req(input_data()$data$gp_result)
      gses <- sapply(input_data(), function(x){return(name2idPrefix(x$name))})
      l_out <- lapply(input_data(), function(x){
        
        return(input[[paste0(name2idPrefix(x$name), "_results-plot_rows_selected")]])
      })
      names(l_out) <- gses
      l_out
      
    })
    
  
    
    #hiding elements when sidebar is collapsed
    observe({
      if(input$sidebarCollapsed){
        shinyjs::hide(selector = "#data_sel", anim = T, animType = "slide", time = 0.1)
        shinyjs::hide(selector = "#gene_select", anim = T, animType = "slide", time = 0.1)
        shinyjs::hide(selector = ".sidebar.shiny-bound-input .btn", anim = T, animType = "slide", time = 0.1)
        shinyjs::hide(selector = ".sidebar_head", anim = T, animType = "slide", time = 0.1)
       
      }
      else{
        shinyjs::show(selector = "#data_sel", anim = T, animType = "slide", time = 0.1)
        shinyjs::show(selector = "#gene_select", anim = T, animType = "slide", time = 0.1)
        shinyjs::show(selector = ".sidebar.shiny-bound-input .btn", anim = T, animType = "slide", time = 0.1)
        shinyjs::show(selector = ".sidebar_head", anim = T, animType = "slide", time = 0.1)
      
      }
    })
    
    #_____________________________________

    #GENERAL GENE SELECTION LOGIC
    
    #Passing row selections to gene selection in gp table
    sel_genes <- reactive({
      genes_list <- list()
      
      for(i in 1:length(sel_row())){
        if(!is.null(input_data()[[names(sel_row())[i]]]$data$gp_result)){  
          gp_table <- input_data()[[names(sel_row())[i]]]$data$gp_result$result
          genes_list[[names(sel_row())[i]]] <- as.numeric(unlist(strsplit(gp_table$intersection[sel_row()[[i]]], ",")))
        }
        else{
          
        }
        
        }
      # #TODO: Add plot selections here
      # for(i in 1:length(input_data)){
      #   manh_out <- module_outputs[[paste0(study, "_manhat")]]
      #   print(manh_out())
      # 
      # }
      
      genes_list
    })
    
    #resetting table rows on action button
    observeEvent(input$reset, {
      validate(need(input_data(), label = "Missing input data"))

      for(i in 1:length(input_data())){
        if(!is.null(input_data()[[names(input_data())[i]]]$data$gp_result)){
        #clearing tables
          proxy <-  dataTableProxy(paste0(names(input_data())[i],"_results-plot"))
          proxy %>% selectRows(NULL)
        }
        #redrawing manhat plots
        # js$resetClick(paste0(names(input_data)[i]))
        # 
        # # man_proxy <- plotlyProxy(outputId = paste0(names(input_data)[i], "_manhat-plot"))
        # # plotlyProxyInvoke(man_proxy, method = "react", )
        # 
      }
    })
    
    
    #redrawing graphs based on table selections
    observeEvent(input$gene_show, ignoreInit = T, 
      {
        validate(need(input_data(), label = "Missing input data"))
        tab_type <- strsplit(input$tabs_menu, "_")[[1]][2]
        study <-  strsplit(input$tabs_menu, "_")[[1]][1]     
        switch(tab_type,
               
               gexp = {module_outputs[[paste0(study,"_gexp")]] <- volcano_server(id = paste0(study,"_gexp"), 
                                                                                   gprofiler = input_data()[[study]]$data$gp_result,
                                                                                   dataset = input_data()[[study]]$data$gene_expression,
                                                                                   subset = sel_genes()[[study]])
                      },
               
               biofab = {module_outputs[[paste0(study,"_biofab")]] <- go_fabric_server(id = paste0(study,"_biofab"),
                                                                                       gprofiler = input_data()[[study]]$data$gp_result,
                                                                                       dataset = input_data()[[study]]$data$gene_expression,
                                                                                       subset = sel_row()[[study]])
                      }
               
               )
      
    })
   
    
    
    #_______________________________________________
    
    
    #BUILDING TABSETS
    #all_tabs <- reactiveValues()
    all_tabs <- list()
    observeEvent(input_data(), ignoreInit = T,{
      
      all_tabs <<- list()
      
      output$tabs <- renderUI({
          
          
          for(i in 1:length(input_data())){
              study <- names(input_data())[i]
              #'volcano plot definitions
              
              all_tabs[[paste0(study,"_gexp_tab")]] <<- tabItem(tabName = paste0(study, "_gexp_tab"),
                                                              volcano_ui(id = paste0(study, "_gexp"), 
                                                                         color = studies_palette()[[study]], 
                                                                         study_name = study))
              if(!is.null(input_data()[[study]]$data$gp_result)){
                
              
                all_tabs[[paste0(study,"_manhat_tab")]] <<- tabItem(tabName = paste0(study, "_manhat_tab"),
                                                                  manhat_ui(id = paste0(study, "_manhat"),
                                                                             color = studies_palette()[[study]],
                                                                             study_name = study))
    
                all_tabs[[paste0(study,"_biofab_tab")]] <<- tabItem(tabName = paste0(study, "_biofab_tab"),
                                                                  go_fabric_ui(id = paste0(study, "_biofab"),
                                                                            color = studies_palette()[[study]],
                                                                            study_name = study))
                
                all_tabs[[paste0(study, "_results_tab")]] <<- tabItem(tabName = paste0(study, "_results_tab"),
                                                                      dataTable_ui(id = paste0(study, "_results"),
                                                                                   color = studies_palette()[[study]],
                                                                                   study_name = study)
                                                                      )
              }
          }
          return(tags$div(class = "tab-content", all_tabs))
      })
    })
    
    
    #HIDING TABSETS BASED ON CLICKED ON MENU ITEMS
    observeEvent(input$studies_tab, ignoreNULL = F, ignoreInit = T,{
      #browser()
     if(!is.null(input$studies_tab) && input$studies_tab != "import"){
       
       shinyjs::hide(selector =  "#shiny-tab-import_tab", anim = T, animType = "slide")
       shinyjs::show(selector = "#tabs")
     } else{
       shinyjs::show(selector = "#shiny-tab-import_tab")
       shinyjs::hide(selector =  "#tabs")
     }
      
    }
  )
    
    
    
    #MODULE DEFINITIONS
    module_outputs <- reactiveValues()
    observeEvent(input_data(), ignoreInit = T, {
      validate(need(input_data(), label = "Missing input data"))
      #resetting just in case
      module_outputs <- reactiveValues()
      
      for(i in 1:length(input_data())){
          study <- names(input_data())[i]
          module_outputs[[paste0(study,"_gexp")]] <- volcano_server(id = paste0(study,"_gexp"), 
                                                                    gprofiler = input_data()[[study]]$data$gp_result,
                                                                    dataset = input_data()[[study]]$data$gene_expression
                                                                    )
          if(!is.null(input_data()[[study]]$data$gp_result)){
          
            module_outputs[[paste0(study,"_manhat")]] <- manhat_server(id = paste0(study,"_manhat"),
                                                                   gprofiler = input_data()[[study]]$data$gp_result
                                                                      
                                                                    )
            module_outputs[[paste0(study,"_biofab")]] <- go_fabric_server(id = paste0(study,"_biofab"),
                                                                     gprofiler = input_data()[[study]]$data$gp_result,
                                                                     dataset = input_data()[[study]]$data$gene_expression
                                                                   )
            module_outputs[[paste0(study, "_results")]] <- dataTable_server(id = paste0(study,"_results"),
                                                                            gprofiler = input_data()[[study]]$data$gp_result,
                                                                            dataset = input_data()[[study]]$data$gene_expression
                                                                            )
          }
          
      }
  })
    
    #DOWNLOAD HANDLER
    output$download <- downloadHandler(
      filename = function() {
        paste('AstroExplorer-', Sys.Date(), '.xlsx', sep='')
      },
      content = function(file){
        #browser()
        flat_list <- unlist(unlist(input_data(), recursive = F), recursive = F)
        
        enrichment_data <- lapply(flat_list[grep(pattern = "gp_result", names(flat_list))], function(x) return(x$result))
        expression_data <- flat_list[grep(pattern = "gene_expression", names(flat_list))]
        out <- append(expression_data, enrichment_data)
        return(openxlsx::write.xlsx(out, file))
      }
  
    )

}

# Run the application 
shinyApp(ui = ui, server = server)
