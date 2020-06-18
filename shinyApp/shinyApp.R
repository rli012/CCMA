
setwd('C:\\Users/rli3/Documents/CCMA/')

library(shiny)
library(shinydashboard)
library(shinyjs)
library(readxl)
library(ggplot2)
library(Matrix)
library(stringr)
library(Biobase)
library(survival)
library(survminer)
library(limma)
library(edgeR)

library(shinydashboardPlus)

#source('shinyApp/shiny_functions.R')

##### ui


header=dashboardHeader(title = 'CCMA')

#header=dashboardHeaderPlus(title = 'PCa Transcriptomes', fixed = TRUE)

sidebar=dashboardSidebar(
  
  ### header color
  # https://stackoverflow.com/questions/50467136/how-to-create-custom-shinydashboard-skin
  tags$head(tags$style(HTML('.logo {
                              background-color: #4285F4 !important;
                              }
                              .navbar {
                              background-color: #4285F4 !important;
                              }
                              '))),
  
  sidebarMenu(
    style = 'position:fixed; overflow: visible', # 
    
    menuItem("Home", tabName = 'tab_home', icon = icon("home")),
    menuItem("miRNA Browser", tabName = 'browser', icon = icon("search")),
    menuItem("miRNA Datasets", tabName = 'analysis', icon = icon("database")),
    menuItem("Download", tabName = 'download', icon = icon("download")),
    menuItem("Pipelines", tabName = 'pipelines', icon = icon("dna"))
    
  )
)

tab_home <- fluidRow(
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    width = 12,
    
    h2(strong("Welcome to CCMA (Cancer Circulating miRNA Atlas)")),
    
    hr(),
    
    valueBox(value = '58', color = 'orange', width = 3,
             subtitle = tags$p(strong("Datasets"), style = "font-size: 200%;"), icon = icon("database")),
    valueBox(value = '10,072', color = 'orange', width = 3,
             subtitle = tags$p(strong("Samples"), style = "font-size: 200%;"),  icon = icon("user-circle")),
    valueBox(value = '30', color = 'orange', width = 3,
             subtitle = tags$p(strong("Signatures"), style = "font-size: 200%;"),  icon = icon("dna"))
    
  ),
  
  box(
    title = NULL, solidHeader = TRUE, collapsible = FALSE,
    width = 12, # solidHeader=TRUE can remove the top boarder
    
    h2(strong("Introduction")),
    h3(strong("What is circulating miRNA?")),
    tags$p("miRNAs can function as potential oncogenes or tumor suppressors. 
           Altered expression of these molecules was correlated with the 
           occurrence of many cancer diseases and therefore they 
           are considered a molecular tool for non-invasive cancer diagnosis and prognosis", 
           style = "font-size: 150%;"),
    
    #column(width = 600) {
      
    #}
    tags$img(src='mirna.jpg', width=550), # in www
    tags$img(src='fluid.jpg', width=550),
    
    
    
    br(),
    h3(strong("About CCMA")),
    tags$p('We collected 80 public circulating miRNA expression datasets in cancer', style = "font-size: 150%;"),
    
    br(),
    h3(strong("Cite")),
    tags$p('Please cite the following publication:
           Li, R., et al., CCMA: a human cancer circulating miRNA atlas', style = "font-size: 150%;")
  )
  
)




body=dashboardBody(
  
  useShinyjs(),
  #selectizeInput(inputId = "gene", label='Gene', choices = gene.annotation, selected = gene.default, multiple=FALSE, 
  #               options = list(
  #                 placeholder = 'Select a gene', maxOptions = 10,
  #                 selectOnTab=TRUE)),
  
  #tags$head(tags$style(HTML('.box {margin: 1px;}'))), # distance between box
  
  #https://stackoverflow.com/questions/45706670/shiny-dashboadpage-lock-dashboardheader-on-top
  tags$script(HTML("$('body').addClass('fixed');")), # fix header & sidebar
  
  tabItems(
    tabItem(tabName="tab_home", tab_home)#,
    
    #tabItem(tabName="tab_boxplot", tab_boxplot),
    #tabItem(tabName="tab_kmplot",tab_kmplot),
    #tabItem(tabName="tab_dataset",tab_dataset)
  )
  
)


ui <- dashboardPage(title='CCMA', header, sidebar, body) # skin = 'blue', 


######## Server

server <- function(input, output, session) { 
  
}





shinyApp(
  ui = ui,
  server = server
)

