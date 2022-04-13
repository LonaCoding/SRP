library(shiny)

source("seurat.R")

# User interface ----
ui <- fluidPage(

  # App title ----
  titlePanel("Pipeline2"),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "DimPlot")

    )
  )

server <- function(input, output) {

  output$DimPlot <- renderPlot({


    })

}


# Run app ----
shinyApp(ui, server)
