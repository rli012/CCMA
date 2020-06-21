
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
library(DT)
library(dplyr)
library(pROC)
library(ROCR)


library(dashboardthemes)
library(shinydashboardPlus)

source('shinyApp/shiny_functions.R')

##### ui


############################################################################################

################### Theme

logo_blue_gradient <- shinyDashboardLogoDIY(
  
  boldText = "CCMA",
  mainText = "",
  textSize = 30,
  badgeText = "Beta v1.0",
  badgeTextColor = "white",
  badgeTextSize = 2,
  badgeBackColor = "#40E0D0", # 
  badgeBorderRadius = 10
  
)


### creating custom theme object
theme_blue_gradient <- shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(0,0,0)"
  ,primaryFontColor = "rgb(0,0,0)"
  ,infoFontColor = "rgb(0,0,0)"
  ,successFontColor = "rgb(0,0,0)"
  ,warningFontColor = "rgb(0,0,0)"
  ,dangerFontColor = "rgb(0,0,0)"
  ,bodyBackColor = "rgb(248,248,248)"
  
  ### header
  ,logoBackColor = "rgb(23,103,124)"
  
  ,headerButtonBackColor = "rgb(238,238,238)"
  ,headerButtonIconColor = "rgb(75,75,75)"
  ,headerButtonBackColorHover = "rgb(210,210,210)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(238,238,238)"
  ,headerBoxShadowColor = "#aaaaaa"
  ,headerBoxShadowSize = "2px 2px 2px"
  
  ### sidebar
  ,sidebarBackColor = cssGradientThreeColors(
    direction = "down"
    ,colorStart = "rgb(20,97,117)"
    ,colorMiddle = "rgb(56,161,187)"
    ,colorEnd = "rgb(3,22,56)"
    ,colorStartPos = 0
    ,colorMiddlePos = 50
    ,colorEndPos = 100
  )
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = "3px 5px 5px"
  ,sidebarShadowColor = "#aaaaaa"
  
  ,sidebarUserTextColor = "rgb(255,255,255)"
  
  ,sidebarSearchBackColor = "rgb(55,72,80)"
  ,sidebarSearchIconColor = "rgb(153,153,153)"
  ,sidebarSearchBorderColor = "rgb(55,72,80)"
  
  ,sidebarTabTextColor = "rgb(255,255,255)"
  ,sidebarTabTextSize = 16
  ,sidebarTabBorderStyle = "none none solid none"
  ,sidebarTabBorderColor = "rgb(35,106,135)"
  ,sidebarTabBorderWidth = 1
  
  ,sidebarTabBackColorSelected = cssGradientThreeColors(
    direction = "right"
    ,colorStart = "rgba(44,222,235,1)"
    ,colorMiddle = "rgba(44,222,235,1)"
    ,colorEnd = "rgba(0,255,213,1)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorSelected = "rgb(0,0,0)"
  ,sidebarTabRadiusSelected = "0px 20px 20px 0px"
  
  ,sidebarTabBackColorHover = cssGradientThreeColors(
    direction = "right"
    ,colorStart = "rgba(44,222,235,1)"
    ,colorMiddle = "rgba(44,222,235,1)"
    ,colorEnd = "rgba(0,255,213,1)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorHover = "rgb(50,50,50)"
  ,sidebarTabBorderStyleHover = "none none solid none"
  ,sidebarTabBorderColorHover = "rgb(75,126,151)"
  ,sidebarTabBorderWidthHover = 1
  ,sidebarTabRadiusHover = "0px 20px 20px 0px"
  
  ### boxes
  ,boxBackColor = "rgb(255,255,255)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "0px 1px 1px"
  ,boxShadowColor = "rgba(0,0,0,.1)"
  ,boxTitleSize = 16
  ,boxDefaultColor = "rgb(210,214,220)"
  ,boxPrimaryColor = "rgba(44,222,235,1)"
  ,boxInfoColor = "rgb(210,214,220)"
  ,boxSuccessColor = "rgba(0,255,213,1)"
  ,boxWarningColor = "rgb(244,156,104)"
  ,boxDangerColor = "rgb(255,88,55)"
  
  ,tabBoxTabColor = "rgb(255,255,255)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(0,0,0)"
  ,tabBoxTabTextColorSelected = "rgb(0,0,0)"
  ,tabBoxBackColor = "rgb(255,255,255)"
  ,tabBoxHighlightColor = "rgba(44,222,235,1)"
  ,tabBoxBorderRadius = 5
  
  ### inputs
  ,buttonBackColor = "rgb(245,245,245)"
  ,buttonTextColor = "rgb(0,0,0)"
  ,buttonBorderColor = "rgb(200,200,200)"
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(235,235,235)"
  ,buttonTextColorHover = "rgb(100,100,100)"
  ,buttonBorderColorHover = "rgb(200,200,200)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(200,200,200)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(200,200,200)"
  
  ### tables
  ,tableBackColor = "rgb(255,255,255)"
  ,tableBorderColor = "rgb(240,240,240)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
)


header=dashboardHeader(title = logo_blue_gradient)

sidebar=dashboardSidebar(
  
  sidebarMenu(
    style = 'position:fixed; overflow: visible', # 
    
    menuItem("Home", tabName = 'tab_home', icon = icon("home")),
    menuItem("miRNA Browser", tabName = 'tab_browser', icon = icon("search")),
    menuItem("miRNA Datasets", tabName = 'tab_datasets', icon = icon("database")),
    menuItem("Download", tabName = 'tab_download', icon = icon("download")),
    menuItem("Pipelines", tabName = 'tab_pipelines', icon = icon("dna"))
    
  )
)

tab_home <- fluidRow(
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    width = 12, 
    
    h2(strong("Welcome to CCMA (Cancer Circulating miRNA Atlas)")),
    
    hr(),
    
    valueBox(value = '21', color = 'teal', width = 3,
             subtitle = tags$p(strong("Cancer types"), style = "font-size: 200%;"),  icon = icon("dna")),
    valueBox(value = '80', color = 'teal', width = 3,
             subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
    valueBox(value = '26,072', color = 'teal', width = 3,
             subtitle = tags$p(strong("Samples"), style = "font-size: 200%;"),  icon = icon("user-circle"))
    
    
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


#################################################################################

mir.default <- 'MIMAT0000062' # hsa-let-7a-5p

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')
#colnames(mir.annotation) <- c('mir_id', 'mir_name', 'mir_seq')

mir.id <- selectizeInput(inputId = "mir.id", label=strong('miRNA'),# h4(strong('miRNA'))
                         choices = NULL, selected = mir.default, 
                         multiple = FALSE, width = 300,
                         options = list(placeholder = 'Select a miRNA',
                                        server = TRUE, selectOnTab=TRUE,
                                        searchField = c('Name', 'ID', 'Previous_ID'),
                                        labelField = "Name",
                                        valueField = "ID",
                                        #maxOptions = 5,
                                        render = I("{option: function(item, escape) 
                                                      {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                                         gene = gene + '<li>' + 'Previous IDs:' + item.Previous_ID + '</li>';
                                                         gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                                         return gene
                                                      }
                                                    }")
                          ))





project.default <- 'TCGA-CHOL' # hsa-let-7a-5p

project.id <- selectizeInput(inputId = "project.id", label=NULL,# h4(strong('miRNA'))
                         choices = NULL, selected = project.default, 
                         multiple = FALSE, width = 150,
                         options = list(placeholder = 'Select a project',
                                        server = TRUE, selectOnTab=TRUE#,
                                        #searchField = c('Name', 'ID', 'Previous_ID'),
                                        #labelField = "Name",
                                        #valueField = "ID",
                                        #maxOptions = 5,
                                        #render = I("{option: function(item, escape) 
                                        #           {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                        #           gene = gene + '<li>' + 'Previous IDs:' + item.Previous_ID + '</li>';
                                        #           gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                        #           return gene
                                        #           }
                                        #           }")
                          ))


datasets <- readRDS('shinyApp/data/CCMA_Datasets.RDS')

eSet <- readRDS('shinyApp/data/GSE73002_GPL18941_eSet.RDS')
expr <- exprs(eSet)
meta <- pData(eSet)

colnames(meta)[9] <- 'Sample.Type'

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
expr.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

projects.tcga <- meta.tcga %>% group_by(project_id) %>% 
  summarise(group=length(unique(sample_type)))

projects.tcga <- sort(projects.tcga[projects.tcga$group==2,]$project_id)


tab_browser <- fluidRow(
  box(
    title = 'Select a miRNA', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12, 
    
    mir.id,
    hr(),

    strong(uiOutput("mir.name")),
    strong(textOutput("mir.preid")),
    strong(textOutput("mir.info")),
    strong(textOutput("mir.seq")),
    strong(uiOutput("mir.targets")),
    
    hr(),
    
    column(8,
           br(),
      plotOutput('tcga_boxplot',width = 800, height = 350)
    ),
    
    column(4, 
           project.id, 
           plotOutput('tcga_rocplot',width = 400, height = 300)
    )
    
    #strong("Targets: ", a("ENCORI", href = textOutput('mir.encori'), style = "font-size: 100%;"))

    ),
  
  #box(
  #  title = 'miRNA expression in TCGA', status = "primary", solidHeader = FALSE, collapsible = FALSE,
  #  width = 12, 
    
  #  plotOutput('tcga_boxplot',width = 1000, height = 400)
  #),
  
  
  
  box(
    title = 'Select a dataset', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12, 
    
    DT::dataTableOutput("browser_datasets"),
    hr(),
    plotOutput('mir_boxplot',width = 400, height = 400)
  )
  
  #box(title = NULL,
  #    status = "primary", solidHeader = FALSE, collapsible = TRUE,
  #    width = 4,
  #    height = 400,
  #    plotOutput('mir_boxplot')
  #),
  
  
  )


tab_datasets <- fluidRow(
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    width = 12, 
    
    DT::dataTableOutput("datasets")
    
  ),
  
  
  box(
    title = 'Dataset-level Analysis', status = "primary", solidHeader = TRUE, collapsible = TRUE,
    width = 12,
    
    tabBox(width = 12, id = 'degbox',
           
           tabPanel("Summary", 
                    #column(12, 
                    #),
                    
                    #box(title = 'Summary',  
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #     width = 12,
                    #     #height = 400,
                    #     textOutput("dataset_summary"),
                    #     
                    #     div("div creates segments of text with a similar style. This division of text is all blue because I passed the argument 'style = color:blue' to div", style = "color:blue"),
                    #     br(),
                    #     p("span does the same thing as div, but it works with",
                    #       span("groups of words", style = "color:blue"),
                    #       "that appear inside a paragraph.")
                    # ),
                    
                    div(h4(strong(textOutput("dataset_summary"))), style = "color:black", align='center'),
                    
                    #div(h4(strong(dataset_summary)), style = "color:black", align='center'),
                    br(),
                    #p(span("Experimental design:", style = 'color:blue'), overall_design),
                    
                    
                    #conditionalPanel(condition = "input.datasets_rows_selected==1 || input.datasets_rows_selected==3 || input.datasets_rows_selected==6 || 
                    #                 input.datasets_rows_selected==8 || input.datasets_rows_selected==9",
                    #                 box(title = 'Sample Type',  
                    #                     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #                     width = 4,
                    #                     height = 400,
                    #                     plotOutput('pie_sample_type')
                    #                 )
                    #),
                    
                    box(title = 'Sample Type',
                        status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        width = 4,
                        height = 400,
                        plotOutput('pie_sample_type')
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected==1 || input.datasets_rows_selected==5 || input.datasets_rows_selected==12",
                                     box(title = 'Pathological T Stage', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('pie_pstage')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected==1 || input.datasets_rows_selected==2 || input.datasets_rows_selected==3",
                                     box(title = 'Clinical T Stage', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('pie_cstage')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected!=3 && input.datasets_rows_selected<10",
                                     box(title = 'Preoperative PSA', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('bar_psa')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected<12 && input.datasets_rows_selected!=2 && input.datasets_rows_selected!=10 && 
                                     input.datasets_rows_selected!=11",
                                     box(title = 'Gleason Score', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('bar_gleason')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected==1",
                                     box(title = 'Overall Survival Status', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('pie_os_status')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected==1",
                                     box(title = 'Overall Survival', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('km_os_time')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected<=12",
                                     box(title = 'Relapse-free Survival Status', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('pie_bcr_status')
                                     )
                    ),
                    
                    
                    conditionalPanel(condition = "input.datasets_rows_selected<=10",
                                     box(title = 'Relapse-free Survival', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('km_bcr_time')
                                     )
                    ),
                    
                    
                    conditionalPanel(condition = "input.datasets_rows_selected==10 || input.datasets_rows_selected==12 || input.datasets_rows_selected==13 || 
                                     input.datasets_rows_selected==14",
                                     box(title = 'Metastasis-free Survival Status', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('pie_metastasis_status')
                                     )
                    ),
                    
                    conditionalPanel(condition = "input.datasets_rows_selected==10",
                                     box(title = 'Metastasis-free Survival', 
                                         status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         width = 4,
                                         height = 400,
                                         plotOutput('km_metastasis_time')
                                     )
                    )#,
                    
                    
                    #box(title = 'Preop PSA', 
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #     width = 4,
                    #     #height = 400,
                    #     plotOutput('pie_psa')
                    # ),
                    
                    
                    #htmlOutput("gse")
                    
                    #splitLayout(cellWidths = c("50%", "50%"), 
                    #             plotOutput('pie_sample_type'), 
                    #             plotOutput('pie_gleason')),
  ),
  
  tabPanel(title = "Differential Expression Analysis", value = 'sample_type',
           checkboxGroupInput(inputId = "control_sample_type", label = "Control group",
                              choices = c('Normal', 
                                          'Primary Tumor' = 'Primary',
                                          'Metastasis'),
                              selected = 'Normal',
                              inline = TRUE),
           checkboxGroupInput(inputId = "case_sample_type", label = "Case group",
                              choices = c('Normal',
                                          'Primary Tumor' = 'Primary',
                                          'Metastasis'),
                              selected = 'Primary',
                              inline = TRUE),
           
           
           sliderInput(inputId = "foldchange", label = h5(strong('Fold Change')), 
                       min = 0, max = 3,  step = 0.1, value = 2, width = 300),
           
           sliderInput(inputId = "fdr", label = h5(strong('BH Adjusted P Value')), 
                       min = 0, max = 0.1,  step = 0.01, value = 0.01, width = 300),
           
           
           column(4, br(), plotOutput('volcano_sample_type', height = 500)),
           column(6, br(), DT::dataTableOutput("table_sample_type"))
           
           
           
  ),
  
  tabPanel("Highly Expressed miRNAs", 
           radioButtons(inputId = "gleason_control", label = "Control group",
                        c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                        inline = TRUE),
           radioButtons(inputId = "gleason_case", label = "Case group",
                        c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                        inline = TRUE)
  ),
  
  tabPanel("Hierarchical Clustering", 
           
           sliderInput("psa", 
                       "Preoperative PSA Value", 
                       min = 0,
                       max = 50, 
                       value = 4,
                       width = 500)
  ),
  
  tabPanel("PC Analysis", 
           
           sliderInput("psa", 
                       "Preoperative PSA Value", 
                       min = 0,
                       max = 50, 
                       value = 4,
                       width = 500)
  ),
  
  
  
  tabPanel("ROC Analysis", 
           radioButtons(inputId = "stage_control", label = "Control group",
                        c('Stage I', 'Stage II','Stage III','Stage IV'),
                        inline = TRUE),
           radioButtons(inputId = "stage_case", label = "Case group",
                        c('Stage I', 'Stage II','Stage III','Stage IV'),
                        inline = TRUE)
  )
  
  )
  
  )
  
)



tab_download <- fluidRow(
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    width = 12, 
    
    DT::dataTableOutput("download")
    
  )
  
)


tab_pipelines <- fluidRow(
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    width = 12#, 
    
    
    
  )
  
)



body=dashboardBody(
  
  useShinyjs(),
  
  tags$script(HTML("$('body').addClass('fixed');")), # fix header & sidebar
  
  theme_blue_gradient,
  
  tabItems(
    tabItem(tabName="tab_home", tab_home),
    
    tabItem(tabName="tab_browser", tab_browser),
    tabItem(tabName="tab_datasets",tab_datasets),
    tabItem(tabName="tab_download",tab_download)
    
  )
  
)


ui <- dashboardPage(title='CCMA', header, sidebar, body) # skin = 'blue', 


######## Server

server <- function(input, output, session) { 
  updateSelectizeInput(session, 'mir.id', choices = mir.annotation, selected = mir.default, server = TRUE)
  updateSelectizeInput(session, 'project.id', choices = projects.tcga, selected = project.default, server = TRUE)
  
  
  output$mir.name <- renderUI({ 
    mir.id <- input$mir.id
    mir.name <- mir.annotation[mir.id, 'Name']
    mir.url <- paste0('http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=', mir.id)
    mir.url <- a(mir.name, href = mir.url, style = "font-size: 150%;")
    tagList(mir.url)
  })
  
  output$mir.preid <- renderText({ 
    mir.id <- input$mir.id
    mir.preid <- mir.annotation[mir.id, 'Previous_ID']
    mir.preid <- paste0('Previous IDs: ', mir.preid)
    mir.preid
  })
  
  output$mir.info <- renderText({ 
    mir.id <- input$mir.id
    mir.name <- mir.annotation[mir.id, 'Name']
    mir.info <- paste0('Accession: ', mir.id)
    mir.info
  })
  
  output$mir.seq <- renderText({ 
    mir.id <- input$mir.id
    mir.seq <- mir.annotation[mir.id, 'Sequence']
    mir.seq <- paste0('Sequence: ', mir.seq)
    mir.seq
  })
  
  
  output$mir.targets <- renderUI({ 
    mir.id <- input$mir.id
    mir.name <- mir.annotation[mir.id, 'Name']
    mir.encori <- paste0('http://starbase.sysu.edu.cn/agoClipRNA.php?source=mRNA&flag=miRNA&clade=mammal&genome=human&assembly=hg19&miRNA=',
                         mir.name, '&clipNum=&deNum=&panNum=&proNum=&program=&target=')
    mir.encori <- a('ENCORI', href = mir.encori, style = "font-size: 100%;")
    
    mir.mirdb <- paste0('http://mirdb.org/cgi-bin/search.cgi?searchType=miRNA&full=mirbase&searchBox=',mir.id)
    mir.mirdb <- a('miRDB', href = mir.mirdb, style = "font-size: 100%;")
    
    mir.mirtarbase <- paste0('http://mirtarbase.cuhk.edu.cn/php/search.php?opt=b_mirna&org=hsa&bname=',mir.name)
    mir.mirtarbase <- a('miRTarBase', href = mir.mirtarbase, style = "font-size: 100%;")
    
    mir.targetscan <- paste0('http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi?mirg=',mir.name)
    mir.targetscan <- a('TargetScan', href = mir.targetscan, style = "font-size: 100%;")
    
    mir.dianatarbase <- paste0('http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=',mir.name)
    mir.dianatarbase <- a('Diana-TarBase', href = mir.dianatarbase, style = "font-size: 100%;")
    
    tagList("Targets:", mir.encori, mir.mirdb, mir.mirtarbase, mir.targetscan, mir.dianatarbase)
  })
  
  output$tcga_boxplot <- renderPlot({
    
    mir.id <- input$mir.id
    mir.name <- mir.annotation[mir.id, 'Name']
    #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
    
    #grp <- group.expression()
    group <- meta.tcga[,'sample_type']
    expr <- expr.tcga[mir.id,]
    project <- meta.tcga[,'project_id']
    
    dataForBoxPlot <- data.frame(expr, group, project, mir.name)
    
    p <- tcgaboxplotFun(dataForBoxPlot)
    p
  })
  
  
  output$tcga_rocplot <- renderPlot({
    
    mir.id <- input$mir.id
    mir.name <- mir.annotation[mir.id, 'Name']
    #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
    
    #grp <- group.expression()
    project.id <- input$project.id
    idx <- meta.tcga$project_id==project.id
    
    group <- meta.tcga[idx,'sample_type']
    expr <- expr.tcga[mir.id,idx]
    
    dataForROCPlot <- data.frame(expr, group)
    dataForROCPlot$group <- factor(dataForROCPlot$group, levels = c('Normal','Tumor'))
    
    p <- rocplotFun(dataForROCPlot)
    p
  })
  
  
  
  
  output$browser_datasets <- DT::renderDataTable({datasets},
                                                 options = list(pageLength = 5),
                                                 selection = list(mode='multiple', selected=1)
  )
  
  
  output$mir_boxplot <- renderPlot({
    
    mir.id <- input$mir.id
    #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
    
    #grp <- group.expression()
    group <- meta[,'Sample.Type']
    expr <- expr[mir.id,]
    dataset <- 'GSE73002'
    
    dataForBoxPlot <- data.frame(expr, group, dataset)

    p <- boxplotFun(dataForBoxPlot)
    p
  })
  
  
  output$datasets <- DT::renderDataTable({datasets},
                                         options = list(pageLength = 5),
                                         selection = list(mode='single', selected=1)
  )
  
  output$pie_sample_type <- renderPlot({
    #idx <- input$dataset_rows_selected
    #accession <- as.character(datasets[idx,'Accession'])
    
    sample.freq <- table(meta$Sample.Type)
    dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
    
    p <- pieplotFun(dataForPiePlot)
    
    p
  })
  
  
  output$dataset_summary <- renderText({ 
    idx <- input$datasets_rows_selected
    dataset_summary <- as.character(paste0(datasets[idx,'Accession'], ': ', datasets[idx,'Title']))
    dataset_summary
  })
  
  
  
  
  
  
  output$download <- DT::renderDataTable({datasets},
                                        options = list(pageLength = 100, fixedHeader=TRUE,
                                                       #autoWidth = TRUE,
                                                       #columnDefs = list(list(width = '50px', targets = c(2,3,4,6))), # "_all"
                                                       initComplete = JS(
                                                         "function(settings, json) {",
                                                         "$(this.api().table().header()).css({'background-color': '#20B2AA', 'color': '#fff'});",
                                                         "}")),
                                        selection = list(mode='none')
  )
  
}


shinyApp(
  ui = ui,
  server = server
)


