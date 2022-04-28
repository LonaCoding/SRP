#loads in libraries for apps
library(shiny)
library(shinythemes)
library(Seurat)
library(uwot)

#opens and reads the rds file based on seurat analysis into rshiny app
set_1 <- readRDS("/home/fic2/Desktop/SRP_APP/pipeline2_seurat_object.rds")
data_list = list(set_1=set_1)

#Sets the page design for the rshiny app
ui <- fluidPage(
  
#main panel of rshiny app with text
                    mainPanel(
                    #buttons to display information regarding app
                    radioButtons("dataset", label = h3("Gene Cluster Identifier"),
                                 choices = list("Pipeline 2" = "set_1"),
                                 selected = "set_1"),
                    helpText("Enter gene names in CAPS.  Gene names must be exact e.g. CD24 otherwise you will recieve an Error"),
                    #allows user to input gene name for app interface to work
                    textInput("gene2", label = "Gene Name", value = "CD24"),
                    #defines cluster plots for main panel
                    plotOutput("dimPlot1", width = "75%"),
                    plotOutput("genePlot2", width = "75%"),
                  )
           )
#server functionality of app  
server <- function(input, output, session) {
  #inputs dataset for app
  datasetInput <- reactive({ df <- data_list[[input$dataset]]
  })
  #server fucntionality to dispaly both clsuter plots
  output$dimPlot1 <- renderPlot({DimPlot(datasetInput(), reduction = "umap", label =T) })
  output$genePlot2 <- renderPlot({FeaturePlot(datasetInput(), features = (input$gene2), reduction = "umap")})
  
 
}
#call to shiny app function
shinyApp(ui, server)
