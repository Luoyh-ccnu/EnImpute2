# library(reticulate)
# use_condaenv("D:/anaconda/envs/enimpute2_py37/python.exe")
# py_config()
# ui.R
library(shiny)
library(shinydashboard)
library(shinythemes)
library(lemon)
library(patchwork)
library(ggplot2)
library(EnImpute2)
library(Seurat)
source("shiny_tools.R")
ui = navbarPage(
  "EnImpute2 Imputation",
  theme = shinytheme("cerulean"),
  tabPanel("Home",
           fluidPage(
             titlePanel(h1("Welcome to EnImpute2 shiny platform", align = 'center')),
             sidebarLayout(
               sidebarPanel(
                 h4('Welcome to EnImpute2 platform where you can use to have a quick experience with EnImpute2.'),
                 br(),
                 'The right picture briefly summaries the workflow of EnImpute2, it is based on the ensemble idea of EnImpute which uses trimmed mean as an ensembling strategy.',
                 br(),
                 br(),
                 'Furthermore,EnImpute2 borrows the change point detection idea from Time Series Analysis to select the most suitable to impute the corresponding raw matrix.',
                 'For more details,please refer to User Guide section',
                 width = 6
               ),
               mainPanel(img(
                 src = "mindmapv2.jpg",
                 width = 540,
                 height = 567
               ),width=6)
             )
           )),
  tabPanel("Imputation",
           sidebarLayout(
             sidebarPanel(
               fileInput("count",
                         label = "Count matrix",
                         multiple = FALSE),
               
               numericInput(
                 "scale.factor",
                 label = "scale.factor",
                 value = 10000,
                 min =  1,
                 step = 1000
               ),
               numericInput(
                 "trim",
                 label = "trim",
                 value = 0.3,
                 min =  0,
                 max = 0.5,
                 step = 0.1
               ),
               numericInput(
                 "threshold",
                 label = "threshold",
                 value = 0.7,
                 min =  0,
                 max = 1,
                 step = 0.1
               ),
               
               checkboxInput("ALRA", "ALRA", value = TRUE),
               checkboxInput("DCA", "DCA", value = TRUE),
               checkboxInput("DrImpute", "DrImpute", value = TRUE),
               checkboxInput("knn_smooth", "knn_smooth", value = TRUE),
               checkboxInput("MAGIC", "MAGIC", value = TRUE),
               checkboxInput("SAVER", "SAVER", value = TRUE),
               checkboxInput("scImpute", "scImpute", value = TRUE),
               checkboxInput("scNPF", "scNPF", value = TRUE),
               checkboxInput("SCRABBLE", "SCRABBLE", value = TRUE),
               checkboxInput("scRMD", "scRMD", value = TRUE),
               checkboxInput("scTSSR", "scTSSR", value = TRUE),
               checkboxInput("scTSSR2", "scTSSR2", value = TRUE),
               checkboxInput("SDImpute", "SDImpute", value = TRUE),
               checkboxInput("VIPER", "VIPER", value = TRUE),
               checkboxInput("zinbwave", "zinbwave", value = TRUE),
               
               actionButton("Run_EnImpute2",
                            label = "Run EnImpute2"),
               
               width = 2
             ),
             
             mainPanel(fluidRow(
               column(9,
                      textOutput("summary1"),
                      br(),
                      textOutput("summary2")),
               column(
                 3,
                 radioButtons(
                   "scale",
                   label = "Output File Scale",
                   choices = c(
                     "log (logarithmic scale)" = "exp",
                     "exp (Exponential  scale)" = "log"
                   ),
                   selected = "exp"
                 ),
                 #设置单选按钮
                 radioButtons(
                   "fileformat",
                   label = "Output File Format",
                   choices = c(
                     ".txt (tab-delimited text)" = "txt",
                     ".csv (comma-separated values)" = "csv"
                   ),
                   selected = "csv"
                 ),
                 textInput("dlname",
                           label = "Output File Name (Do not include file extension)"),
                 downloadButton("download",
                                label = "Download")
               )
             ))
           )),
  navbarMenu(
    "Downstream analysis",
    tabPanel(
      "Recover gene expression",
      sidebarLayout(
        sidebarPanel(
          'This is where you can see the visulization of the recovery experiment',
          br(),
          'When you upload the raw count matrix and the imputed matrix,UMAP dimensionality reduction visulization of two matrices will be shown on the right side.',
          br(),
          'You can adjust the number and the type of individual methods to be considered in EnImpute2 by the visualization image.',
          fileInput("count1",
                    label = "Count matrix",
                    multiple = FALSE),
          fileInput("Imputed1",
                    label = "Imputed matrix",
                    multiple = FALSE),
          actionButton("Recovery_visualiazation",
                       label = "Start Recovery Visualiazation!")
        ),
        mainPanel(plotOutput("recover_plot"))
      )
    ),

    tabPanel("Cell clustering",
             sidebarLayout(
               sidebarPanel(
                 'This is where you can see the result of the clustering experiment',
                 br(),
                 'When you upload the raw count matrix and the imputed matrix, ARI and NMI of two matrices will be shown on the right side.',
                 br(),
                 'You can adjust the number and the type of individual methods to be considered in EnImpute2 by comparing the difference between two pictures.',
                 fileInput("count2",
                           label = "Count matrix",
                           multiple = FALSE),
                 fileInput("Imputed2",
                           label = "Imputed matrix",
                           multiple = FALSE),
                 actionButton("clustering_visualiazation",
                              label = "Start Cluster Visualiazation!")
               ),
               mainPanel(plotOutput("cluster_plot"))
             )),

    tabPanel(
      "Differential analysis",
      sidebarLayout(
        sidebarPanel(
          'This is where you can see the volcano plot of the differential experiment',
          br(),
          'The uploaded matric must contain three columns:
          one is all genes\' name,one is log2FC and the other one is pvalue',
          br(),
          'You can adjust the number and the type of individual methods to be considered in EnImpute2 by the volcano pictures.',
          fileInput("count3",
                    label = "Count matrix(including log2FC and pvalue)",
                    multiple = FALSE),
          fileInput("Imputed3",
                    label = "Imputed matrix(including log2FC and pvalue)",
                    multiple = FALSE),
          actionButton("Differential_visualiazation",
                       label = "Start Differential Visualiazation!")
        ),
        mainPanel(plotOutput("differential_plot"))
      )
    ),

    tabPanel(
      "Trajectory analysis",
      sidebarLayout(
        sidebarPanel(
          'This is where you can see the visulization of the Trajectory experiment',
          br(),
          'When you upload the raw count matrix and the imputed matrix,lineage construction of two matrices will be shown on the right side.',
          br(),
          'You can adjust the number and the type of individual methods to be considered in EnImpute2 by the visualization image.',
          fileInput("count4",
                    label = "Count matrix",
                    multiple = FALSE),
          fileInput("Imputed4",
                    label = "Imputed matrix",
                    multiple = FALSE),
          actionButton("Trajectory_visualiazation",
                       label = "Start Trajectory Visualiazation!")
        ),
        mainPanel(plotOutput("tra_plot"))
      )
    )
  ),

  tabPanel("User Guide", includeMarkdown("README.md")),
  
  tabPanel(
    "Contact us",
    # img(src = "pic.jpg",
    #   width = 590,
    #   height = 826,
    #   text-align='center'
    #   # margin: -70px 0 0 -70px
    # ),
    h3(
      'If you have any questions,please do not hesitate to send emails to luoyuheng@mails.ccnu.edu.cn',
      align = 'center'
    )
  ),
  
  tabPanel("Reference", includeMarkdown("enimpute2_ref.md"))
)


options(shiny.maxRequestSize=100*1024^2)
server = function(input, output, session) {
  observeEvent(input$Run_EnImpute2, {
    count = read.csv(
      file = input$count$datapath,
      header = TRUE,
      row.names = 1
    )
    count = as.matrix(count)
    
    time.used = system.time({
      Impute.count = try(EnImpute2(
        count,
        scale.factor = input$scale.factor,
        trim = input$trim,
        threshold = input$threshold,
        ALRA = input$ALRA,
        DCA = input$DCA,
        DrImpute = input$DrImpute,
        knn_smooth = input$knn_smooth,
        MAGIC = input$MAGIC,
        SAVER = input$SAVER,
        scImpute = input$scImpute,
        scNPF = input$scNPF,
        SCRABBLE = input$SCRABBLE,
        scRMD = input$scRMD,
        scTSSR = input$scTSSR,
        scTSSR2 = input$scTSSR2,
        SDImpute = input$SDImpute,
        VIPER = input$VIPER,
        zinbwave = input$zinbwave
      ))
    })
    
    output$summary1 = renderText("Imputation is finished (EnImpute2)")
    output$summary2 = renderText({
      paste("Time used (seconds):", round(as.numeric(time.used[3]), 2))
    })
    output$download = downloadHandler(
      filename = function() {
        paste(input$dlname, ".", input$fileformat, sep = "")
      },
      
      content = function(file) {
        if (input$fileformat == "txt") {
          if (input$scale == "exp")
            write.table(Impute.count$count.EnImpute2.exp,
                        file,
                        sep = "\t",
                        quote = FALSE)
          else
            write.table(Impute.count$count.EnImpute2.log,
                        file,
                        sep = "\t",
                        quote = FALSE)
        }
        else{
          if (input$scale == "exp")
            write.csv(Impute.count$count.EnImpute2.exp, file, quote = FALSE)
          else
            write.csv(Impute.count$count.EnImpute2.log, file,  quote = FALSE)
        }
      }
    )
    
  })
  
  observeEvent(input$Recovery_visualiazation, {
    count = read.csv(
      file = input$count1$datapath,
      header = TRUE,
      row.names = 1, check.names=FALSE
    )
    count = as.matrix(count)

    count.EnImpute2.exp= read.csv(
      file = input$Imputed1$datapath,
      header = TRUE,
      row.names = 1, check.names=FALSE
    )
    count.EnImpute2.exp = as.matrix(count.EnImpute2.exp)
    
    output$recover_plot = renderPlot({
      count.tsne = UMAP_cluster(count)
      count.EnImpute2.tsne = UMAP_cluster(count.EnImpute2.exp)
      p1=DimPlot(count.tsne, reduction = "umap")+ggtitle('Before EnImpute2')
      p2=DimPlot(count.EnImpute2.tsne, reduction = "umap")+ggtitle('After EnImpute2')
      p1+p2

    })
  })
  
  observeEvent(input$clustering_visualiazation, {
    count = read.csv(
      file = input$count2$datapath,
      header = TRUE,
      row.names = 1, check.names=FALSE
    )
    count = as.matrix(count)
    
    count.EnImpute2.exp= read.csv(
      file = input$Imputed2$datapath,
      header = TRUE,
      row.names = 1, check.names=FALSE
    )
    count.EnImpute2.exp = as.matrix(count.EnImpute2.exp)
    
    output$cluster_plot = renderPlot({
      count.score = clusterscore(count)
      count.EnImpute2.score = clusterscore(count.EnImpute2.exp)
      p1=ggplot(count.score,aes(x=criteria,y=value,fill=criteria))+
        geom_bar(stat = "identity")+xlab('Before EnImpute2')+ylab('')
      p2=ggplot(count.EnImpute2.score,aes(x=criteria,y=value,fill=criteria))+
        geom_bar(stat = "identity")+xlab('After EnImpute2')+ylab('')
      grid_arrange_shared_legend(p1,p2, nrow = 1,ncol = 2, position='right')
    })
  })
  
  observeEvent(input$Differential_visualiazation, {
    count = read.csv(
      file = input$count3$datapath,
      header = TRUE, check.names=FALSE,
      stringsAsFactors = F
    )
    # count = as.matrix(count)
    
    count.EnImpute2.exp= read.csv(
      file = input$Imputed3$datapath,
      header = TRUE,check.names=FALSE,
      stringsAsFactors = F
    )
    # count.EnImpute2.exp = as.matrix(count.EnImpute2.exp)
    
    output$differential_plot = renderPlot({
      count.dif = volcano_plot(count)+xlab('Before EnImpute2')
      count.EnImpute2.dif = volcano_plot(count.EnImpute2.exp)+xlab('After EnImpute2')
      count.dif+count.EnImpute2.dif
    })
  })
  
  observeEvent(input$Trajectory_visualiazation, {
    count = read.csv(
      file = input$count4$datapath,
      header = TRUE,
      row.names = 1, check.names=FALSE
    )
    count = as.matrix(count)
    
    count.EnImpute2.exp= read.csv(
      file = input$Imputed4$datapath,
      header = TRUE,
      row.names = 1, check.names=FALSE
    )
    count.EnImpute2.exp = as.matrix(count.EnImpute2.exp)
    
    output$tra_plot = renderPlot({
      cat('1')
      count.tra = trascore(count)
      cat('2')
      count.EnImpute2.tra = trascore(count.EnImpute2.exp)
      p1=ggplot(count.tra$score,aes(x=criteria,y=value,fill=criteria))+
        geom_bar(stat = "identity")+
        theme(legend.position ='none')+xlab('Before EnImpute2')
      p2=plot_cell_trajectory(count.tra$dCellData, color_by = "timepoint")+xlab('Before EnImpute2')
      p3=ggplot(count.EnImpute2.tra$score,aes(x=criteria,y=value,fill=criteria))+
        geom_bar(stat = "identity")+
        theme(legend.position ='none')+xlab('After EnImpute2')
      p4=plot_cell_trajectory(count.EnImpute2.tra$dCellData, color_by = "timepoint")+xlab('After EnImpute2')
      (p1+p2)/(p3+p4)
    })
  })
}


shinyApp(ui = ui, server = server)