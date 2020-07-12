
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
library(digest)
library(pheatmap)


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
    valueBox(value = '88', color = 'teal', width = 3,
             subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
    valueBox(value = '31,933', color = 'teal', width = 3,
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


ccma.datasets <- readRDS('shinyApp/data/CCMA_Datasets.RDS')
ccma.primary <- readRDS('shinyApp/data/CCMA_Datasets_Primary.RDS')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
expr.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

projects.tcga <- meta.tcga %>% group_by(project_id) %>% 
  summarise(group=length(unique(sample_type)))

projects.tcga <- sort(projects.tcga[projects.tcga$group==2,]$project_id)


expr.ccma <- readRDS('shinyApp/data/CCMA_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/CCMA_Metadata.RDS')


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
    
    DT::dataTableOutput("analysis_datasets")
    
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
                    #                     plotOutput('pie_disease_status')
                    #                 )
                    #),
                    
                    box(title = 'Disease Status',
                        status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_disease_status')
                    ),
                    
                    box(title = 'Subgroups',
                        status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_group')
                    ),
                    
                    #box(title = 'Preop PSA', 
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #     width = 4,
                    #     #height = 400,
                    #     plotOutput('pie_psa')
                    # ),
                    
                    box(title = NULL,
                        status = "primary", solidHeader = FALSE, collapsible = TRUE,
                        width = 12,
                        height = 600,
                        htmlOutput("gse")
                    )
                    
                    #htmlOutput("gse")
                    
                    
                    #splitLayout(cellWidths = c("50%", "50%"), 
                    #             plotOutput('pie_disease_status'), 
                    #             plotOutput('pie_gleason')),
  ),
  
  tabPanel("Highly Expressed miRNAs", 
           radioButtons(inputId = "gleason_control", label = "Control group",
                        c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                        inline = TRUE),
           radioButtons(inputId = "gleason_case", label = "Case group",
                        c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                        inline = TRUE)
  ),
  
  
  tabPanel(title = "Differential Expression Analysis", value = 'sample_type',
           
           column(6, br(), DT::dataTableOutput("groups")),
           
           # checkboxGroupInput(inputId = "control_sample_type", label = "Control group",
           #                    choices = c('Normal', 
           #                                'Primary Tumor' = 'Primary',
           #                                'Metastasis'),
           #                    selected = 'Normal',
           #                    inline = TRUE),
           # checkboxGroupInput(inputId = "case_sample_type", label = "Case group",
           #                    choices = c('Normal',
           #                                'Primary Tumor' = 'Primary',
           #                                'Metastasis'),
           #                    selected = 'Primary',
           #                    inline = TRUE),
           
           column(1, br()),
           
           column(4, br(), 
                  selectInput("deg.test", "Method", width = 300,
                              c("T test" = "t",
                                "Wilcoxon test" = "wilcox",
                                "Limma" = "limma")),
                  
                  sliderInput(inputId = "foldchange", label = h5(strong('Fold Change')), 
                              min = 0, max = 3,  step = 0.1, value = 2, width = 300),
                  
                  sliderInput(inputId = "fdr", label = h5(strong('BH Adjusted P Value')), 
                              min = 0, max = 0.1,  step = 0.01, value = 0.01, width = 300),
                  
                  br(),
                  
                  actionButton(inputId = 'deg.submit', label = strong('Submit'), icon=icon("check"), width = 300)
                  
                  ),
           
           
           column(12, hr(),
                  column(4, br(), plotOutput('volcano_sample_type', height = 500)),
                  column(6, br(), DT::dataTableOutput("table_sample_type"))
                  )
           
  ),
  
  tabPanel("ROC Analysis", 
           radioButtons(inputId = "stage_control", label = "Control group",
                        c('Stage I', 'Stage II','Stage III','Stage IV'),
                        inline = TRUE),
           radioButtons(inputId = "stage_case", label = "Case group",
                        c('Stage I', 'Stage II','Stage III','Stage IV'),
                        inline = TRUE)
  ),
  
  tabPanel("Feature Selection", 
           radioButtons(inputId = "stage_control", label = "Control group",
                        c('Stage I', 'Stage II','Stage III','Stage IV'),
                        inline = TRUE),
           radioButtons(inputId = "stage_case", label = "Case group",
                        c('Stage I', 'Stage II','Stage III','Stage IV'),
                        inline = TRUE)
  ),
  
  tabPanel("Hierarchical Clustering", 
           
           column(3, #offset = 1,
                  checkboxGroupInput(inputId = "cluster", label =  h5(strong('Cluster')),
                                     choices = c('Row', 'Column'), selected = 'Column',
                                     inline = TRUE),
                  checkboxGroupInput(inputId = "names", label = h5(strong('Name')),
                                     choices = c('Row', 'Column'), selected = NULL,
                                     inline = TRUE),
                  
                  sliderInput(inputId = "font.row", label = h5(strong('Font size (Row)')), 
                              min = 1, max = 20,  value = 10, width = 200),
                  
                  sliderInput(inputId = "font.column", label = h5(strong('Font size (Column)')), 
                              min = 1, max = 20,  value = 14, width = 200),
                  
                  radioButtons(inputId = "angle.column", label = "Angle (Column)",
                               c(0, 45, 90), selected = 45,
                               inline = TRUE)
           ),
           
           
           column(9,
                  plotOutput("heatmap")
           )
           
  ),
  
  tabPanel("PC Analysis", 
           
           column(3
                  
           ),
           
           column(9,
                  plotOutput("pca")
           )
           
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
  
  #ns <- session$ns
  
  seed <- reactiveVal()
  
  ################################################################
  ######################## Information ###########################
  
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
  
  
  #########################################################
  ######################## TCGA ###########################
  
  observeEvent(input$mir.id, {
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
    
    observeEvent(input$project.id, {
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
      
    })
    
  })
    
  
  ########################################################################
  ######################## GENE-LEVEL ANALYSIS ###########################
  
  output$browser_datasets <- DT::renderDataTable({ccma.primary},
                                                 options = list(pageLength = 5),
                                                 selection = list(mode='multiple', selected=1)
  )
  
  
  observeEvent(input$browser_datasets_rows_selected, {
    output$mir_boxplot <- renderPlot({
      
      mir.id <- input$mir.id
      #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
      
      #grp <- group.expression()
      
      idx <- input$browser_datasets_rows_selected
      datasets <- as.character(ccma.datasets[idx,'Dataset'])
      
      group <- c()
      expr <- c()
      dataset <- c()
      for (dt in datasets) {
        group <- c(group, meta.ccma[[dt]][,'Disease.Status'])
        expr <- c(expr, expr.ccma[[dt]][mir.id,])
        dataset <- c(dataset, rep(dt, nrow(meta.ccma[[dt]])))
        
      }
      
      dataForBoxPlot <- data.frame(expr=unlist(expr), group=unlist(group), dataset, 
                                   stringsAsFactors = F)
      
      p <- boxplotFun(dataForBoxPlot)
      p
    })
    
  })
  
  
  
  ###########################################################################
  ######################## DATASET-LEVEL ANALYSIS ###########################
  
  output$analysis_datasets <- DT::renderDataTable({ccma.primary},
                                         options = list(pageLength = 5),
                                         selection = list(mode='single', selected=1)
  )
  
  observeEvent(input$analysis_datasets_rows_selected, {
    
    idx <- input$analysis_datasets_rows_selected
    req(idx)
    
    dataset <- as.character(ccma.datasets[idx,'Dataset'])
    
    meta <- meta.ccma[[dataset]]
    expr <- expr.ccma[[dataset]]
    
    seed(sample(1e9,1))
    
    groups <- meta$Group
    
    group.levels <- unlist(sapply(ccma.datasets[idx,'Group'], 
                                  function(x) strsplit(x, '; ')[[1]]))
    groups <- factor(groups, levels = group.levels)

    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    groups <- cbind(t(sapply(groups[,1], function(x) sprintf('<input type="radio" name="%s" value="%s"/>', digest(x,algo='murmur32',seed=seed()), 1:2))), groups, stringsAsFactors=F)
    colnames(groups) <- c('Case','Control','Groups','N')
    
    #groups <- groups[order(groups$N,decreasing=T),]
    #print(groups)
    #TODO too many groups need different approach. (issue is too much overhead to client side and it's not intuitive
    #groups <- groups[seq(min(nrow(groups),100)),]
    #print(dim(groups))
    #shownGroups(groups)
    output$groups <- DT::renderDataTable(groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                         options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                      initComplete = JS("
                                                                        function(setting, json) {
                                                                        $(this.api().table().container())
                                                                        .find('div.dataTables_paginate')
                                                                        .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                        }"))
                                         )
          # drawCallback = JS("
          #                   function(settings) {
          #                   Shiny.unbindAll(this.api().table().node());
          #                   Shiny.bindAll(this.api().table().node());
          #                   }")
          # ),
          # callback = JS("
          #               table.rows().every(function(i, tab, row) {
          #               var $this = $(this.node());
          #               //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
          #               //console.log($this.children()[0]);
          #               //console.log($.parseHTML(this.data()[0])[0].name);
          #               $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
          #               //$this.attr('id', this.data()[2]); //Group Name
          #               $this.addClass('shiny-input-radiogroup');
          #               //console.log($this.prop('checked'));
          #               });
          #               Shiny.unbindAll(table.table().node());
          #               Shiny.bindAll(table.table().node());
          #               "
          #               )
          # 
          # 
          # )
    
    
    output$dataset_summary <- renderText({ 
      dataset_summary <- as.character(paste0(ccma.datasets[idx,'Accession'], ': ', ccma.datasets[idx,'Title']))
      dataset_summary
    })
    
    
    output$gse <- renderUI({
      link <- as.character(ccma.datasets[idx,'Links'])
      tags$iframe(src=link, seamless="seamless", width='100%', height='600')
    })
    
    
    output$pie_disease_status <- renderPlot({
      sample.freq <- table(meta$Disease.Status)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    output$pie_group <- renderPlot({
      sample.freq <- table(meta$Group)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    
    observeEvent(input$deg.submit, {
      
      deg.group <- meta$Group
      
      # group[group %in% input$control_group] <- 'Control'
      # group[group %in% input$case_group] <- 'Case'
      
      deg.group <- ifelse(deg.group=='Healthy', 'Control', 'Case')
      
      deg.group <- factor(deg.group)
      
      design <- model.matrix(~0+deg.group)
      colnames(design) <- levels(deg.group)
      
      contrast.matrix <- makeContrasts(contrasts='Case - Control',
                                       levels=design)
      contrast.matrix
      
      ### Differential gene expression analysis (limma)
      
      fit <- lmFit(expr, design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      
      dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
      
      dataForVolcanoPlot <- dgeTable
      
      logFcThreshold <- log2(input$foldchange)
      adjPvalThreshold <- input$fdr
      
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                           logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                           logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                           logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'
      
      
      output$volcano_sample_type <- renderPlot({
        
        p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
          #xlim(-2,2) +
          labs(x=expression('log'[2]*'(Fold Change)'), 
               y=(expression('-log'[10]*'(FDR)')), 
               title=NULL) +
          geom_point(aes(color=Significance), alpha=1, size=2) +
          geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
                     color='darkgreen', linetype='dashed') +
          geom_hline(yintercept = -log10(adjPvalThreshold), 
                     color='darkgreen',linetype='dashed')+
          #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
          #scale_y_continuous(expand = c(0.3, 0)) +
          #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
          scale_color_manual(values = c(google.green,"gray", google.red)) +
          #facet_wrap(~Comparison, ncol = 2) +
          #geom_text_repel(data = subset(dataForVolcanoPlot, 
          #                              adj.P.Val < adjPvalThreshold & logFC > logFcThreshold), 
          #                segment.alpha = 0.4, aes(label = Symbol), 
          #                size = 3.5, color='red', segment.color = 'black') +
          #geom_text_repel(data = subset(dataForVolcanoPlot, 
          #                              adj.P.Val < adjPvalThreshold & logFC < logFcThreshold*-1), 
          #                segment.alpha = 0.4, aes(label = Symbol), 
          #                size = 3.5, color='green3', segment.color = 'black') +
          
          theme_bw() +
          theme(axis.line = element_blank(),
                #panel.grid.major = element_blank(),
                #panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) +
          theme(legend.position="none") +
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=16),
                strip.text = element_text(size=14, face='bold')) +
          theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
        
        
        p
        
        
      })
      
      
      output$table_sample_type <- DT::renderDataTable({
        
        dataForVolcanoPlot[,-ncol(dataForVolcanoPlot)] <- apply(dataForVolcanoPlot[,-ncol(dataForVolcanoPlot)], 2, 
                                                                function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
        
      })
      
    })
    
    output$heatmap <- renderPlot({
      #req(input$file.upload)
      
      idx <- 1:50
      
      dataForHeatmap <- t(scale(t(expr[1:50,])))
      
      sample.annotation <- data.frame(Group=meta$Group, Disease.Status=meta$Disease.Status,
                                      row.names = colnames(dataForHeatmap), stringsAsFactors = F)
      sample.annotation$Group <- factor(sample.annotation$Group, levels = group.levels)
      
      mx <- max(dataForHeatmap, na.rm = T)
      
      cluster.row <- cluster.col <- name.row <- name.col <- F
      
      if ('Row' %in% input$cluster) {
        cluster.row = T
      }
      
      if ('Column' %in% input$cluster) {
        cluster.col = T
      }
      
      if ('Row' %in% input$names) {
        name.row = T
      }
      
      if ('Column' %in% input$names) {
        name.col = T
      }
      
      p <- pheatmap(dataForHeatmap,
                    scale = 'none',
                    cluster_cols = cluster.col,
                    border_color = NA,
                    cluster_rows = cluster.row,
                    #treeheight_row = 0,
                    show_rownames = name.row,
                    show_colnames = name.col,
                    fontsize_row = input$font.row, 
                    fontsize_col = input$font.column,
                    angle_col = input$angle.column, 
                    annotation_col = sample.annotation,
                    #annotation_legend = T,
                    breaks = c(seq(-1*mx,mx, 2*mx/100)),
                    color=col_fun
      )
      
      p
      
    }) # }, height = 700, width = 700)
    
    
    output$pca <- renderPlot({
      #req(input$file.upload)
      
      filter <- which(rowSums(is.na(expr))>0)

      dataForPCA <- t(scale(t(expr[-filter,])))
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      
      pcDf <- data.frame(PC1=pcaResults$rotation[,1],
                         PC2=pcaResults$rotation[,2],
                         Sample=rownames(pcaResults$rotation),
                         Group=factor(meta$Group, levels=group.levels),
                         Disease.Status=meta$Disease.Status,
                         stringsAsFactors = F)
      
      p <- ggplot(pcDf, aes(PC1, PC2, shape = Disease.Status, color = Group)) + #, shape = "Group"
        theme_bw() +
        #scale_alpha_manual(values = c(0.4, 1)) +
        #scale_color_manual(values = google.colors) +
        #scale_shape_manual(values = shapeScale) +
        #scale_size_manual(values = c(3,7)) +
        geom_point(size = 5) +
        labs(x=paste0('PC1 (', pc1, '%)'), y=paste0('PC2 (', pc2, '%)')) +
        # geom_text_repel(aes(label = extractSampleName(Sample)), size = 2) +
        guides(color = guide_legend(order = 1, override.aes = list(size = 2)),
               shape = guide_legend(order = 2, override.aes = list(size = 4))) +
        #       shape = guide_legend(order = 3, override.aes = list(size = 3)),
        #       size = guide_legend(order = 4)) +
        theme(legend.title = element_text(size = 15, face = 'bold'),
              legend.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 15))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    
    
  })
  
  
  
  
  
  
  
  
  
  
  
  #############################################################
  ######################## DOWNLOAD ###########################
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  download_datatable <- reactive(data.frame(ccma.primary,
                                            ExpressionSet=shinyInput(downloadButton, 88,
                                                                     'button_',
                                                                     label='Download'#,
                                                                     #onclick = sprintf("Shiny.setInputValue('%s', this.id)","select_button")
                                                                     )
                                            )
  )
  
  
  output$download <- DT::renderDataTable({download_datatable()},
                                        options = list(pageLength = 100, fixedHeader=TRUE,
                                                       #autoWidth = TRUE,
                                                       #columnDefs = list(list(width = '50px', targets = c(2,3,4,6))), # "_all"
                                                       initComplete = JS(
                                                         "function(settings, json) {",
                                                         "$(this.api().table().header()).css({'background-color': '#20B2AA', 'color': '#fff'});",
                                                         "}")),
                                        escape = FALSE,
                                        selection = list(mode='none')
  )
  

  lapply(1:88, function(i){
    output[[paste0("button_",i)]] <- downloadHandler(
      filename = function() {
        paste0(ccma.datasets$Accession[i], '_', ccma.datasets$Annotation[i], '_ExpressionSet.RDS')
      },
      content = function(file) {
        exprData <- expr.ccma[[ccma.datasets$Dataset[i]]]
        metaData <- meta.ccma[[ccma.datasets$Dataset[i]]]
        platform <- ccma.datasets$Annotation[i]
        eSet <- ExpressionSet(assayData = as.matrix(exprData),
                              phenoData = AnnotatedDataFrame(metaData),
                              #featureData = annoData,
                              annotation = platform)
        saveRDS(eSet, file)
      }
    )
  })

  
}


shinyApp(
  ui = ui,
  server = server
)


