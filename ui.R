#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme=shinytheme("cerulean"),
    navbarPage(
        "tRuffle",
        tabPanel("About",
                 splitLayout(cellWidths = c("20%","80%"),
                     HTML('<center><img src="tRuffle.png" width="30%" height="30%" ></center>'),
                     tags$h1("Welcome to tRuffle!")
                 )
                 ),
        tabPanel("1) Prepare",
                 sidebarLayout(
                     sidebarPanel(
                         HTML('<left><img src="tRuffle.png" width="18%" height="18%" ></left>'),
                         tags$h3("Data Upload"),
                         fileInput("DEG_table", "Choose your DEG results table", accept = ".xlsx"),
                         tags$h4("Choose column for"),
                         numericInput("col_cluster","cluster ID", min = 1, max = 20,value = 2),
                         numericInput("col_gene","feature", min = 1, max = 20,value = 1),
                         numericInput("col_log","log2FC", min = 1, max = 20,value = 3),
                         numericInput("col_p","p-value", min = 1, max = 20,value = 4),
                         numericInput("col_adj_p","adjusted p-value", min = 1, max = 20,value = 5),
                         width = 2
                         ),
                     mainPanel(tags$h2("Uploaded dataset"),
                               dataTableOutput("DEG_table_out"),
                               width = 10
                     )
                 )),
        tabPanel("2) Explore",
                 sidebarLayout(
                             sidebarPanel(
                                 HTML('<left><img src="tRuffle.png" width="18%" height="18%" ></left>'),
                                 tabsetPanel(
                                     tabPanel("Thresholds",
                                              tags$h3("Choose thresholds"),
                                              sliderInput("p_threshold",label = "max p-value",min = 0.001, max = 0.5, value = 0.01,),
                                              sliderInput("p_adj_threshold",label = "max adj. p-value",min = 0.001, max = 0.5, value = 0.01),
                                              sliderInput("logFC_threshold",label = "absolute value of min log2FC",min = 0.1, max = 5, value = 1),
                                              sliderInput("uniqueness_threshold",label = "min uniqueness",min = 0, max = 1, value = c(0.7,1)),
                                              checkboxGroupInput("considered_thresholds", label = "Considered thresholds", 
                                                                 choices = c("p-value", "adjusted p-value", "log2FC", "uniqueness"),
                                                                 selected = c("p-value", "adjusted p-value", "log2FC", "uniqueness"),
                                                                 inline = TRUE)),
                                     tabPanel("Clusters",
                                              tags$h3("Choose clusters"),
                                              checkboxGroupInput("considered_clusters", label = "Considered clusters for filtering", inline = TRUE)),
                                     tabPanel("Intersection",
                                              tags$h3("Choose clusters"),
                                              checkboxGroupInput("intersection_clusters", label = "Considered clusters for intersection", inline = TRUE))
                                     ), width = 2
                             ),
                             mainPanel(
                                 tabsetPanel(
                                     tabPanel("Table",
                                              tags$h2("Filtered dataset"),
                                       dataTableOutput("DEG_table_filtered_out"),
                                       downloadButton("download_filtered_table", "Download .xlsx")
                                       #textOutput("test"),
                                       #dataTableOutput("test2")
                                       ),
                                     tabPanel("Upset",
                                              tags$h2("Cluster Intersections"),
                                              tags$h4("Please note that intersection sizes depend on your filtering settings."),
                                              tags$text(" "),
                                              plotOutput("upset")
                                              ),
                                     tabPanel("Intersection DEGs",
                                              splitLayout(cellWidths = c("50%","50%"),
                                                  fluidPage(
                                                    tags$h4("Upregulated DEGs"),
                                                    dataTableOutput("DEGs_intersection_up")),
                                                  fluidPage(
                                                    tags$h4("Downregulated DEGs"),
                                                    dataTableOutput("DEGs_intersection_down"),
                                                    downloadButton("download_DEGs", "Download Intersection as .xlsx"))
                                              ))
                         ),width = 10
                         ))))))
                
