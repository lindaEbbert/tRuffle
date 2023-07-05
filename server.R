#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=40*1024^2)
library(shiny)
library(readxl)
library(UpSetR)
library(openxlsx)
library(magrittr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    uniqueness <- function(DEG_clusters,total_clusters){
        score <- (-1/(total_clusters-1))*DEG_clusters+(total_clusters/(total_clusters-1))
        return(score)
    }
    
    DEG_table_input <- reactive({
        req(input$DEG_table)
        readxl::read_xlsx(input$DEG_table$datapath, sheet = 1)
        })


    output$DEG_table_out <- renderDataTable({
        DEG_table_input()
    })
    
    clusters <- reactive({
        vector <- as.data.frame(DEG_table_input())[,input$col_cluster] %>% unique()
        selected <- c()
        for (i in 1:length(vector)){
            selected <- c(selected, as.character(vector[i]))
        }
        if (is.null(selected)){
            selected <- character(0)
        }
        selected
    })
    
    intersection_clusters <- reactive({
        vector <- input$considered_clusters
        selected <- c()
        for (i in 1:length(vector)){
            selected <- c(selected, as.character(vector[i]))
        }
        if (is.null(selected)){
            selected <- character(0)
        }
        selected
    })
    
    observe({
        updateCheckboxGroupInput(session=getDefaultReactiveDomain(), inputId = "considered_clusters",choices = clusters(), selected = clusters(),inline = TRUE)
    })
    
    observe({
        updateCheckboxGroupInput(session=getDefaultReactiveDomain(), inputId = "intersection_clusters",choices = intersection_clusters(), selected = intersection_clusters()[c(1,2)],inline = TRUE)
    })
    
    data1 <- reactive({
        cbind(DEG_table_input(), cluster_up=NA, cluster_down=NA, uniqueness=NA)
    })
    
    totalClusters <- reactive({length(unique(data1()[,input$col_cluster]))})
    
    data2 <- reactive({
        mydata <- data1()
        for (i in 1:length(data1()[,1])){

            if (is.na(data1()$cluster_up[i])||is.na(data1()$cluster_down[i])){
                
                gene <- data1()[i,input$col_gene]
                number_up <- length(unique(data1()[data1()[,input$col_gene]==gene&data1()[,input$col_log]>0,input$col_cluster]))
                number_down <- length(unique(data1()[data1()[,input$col_gene]==gene&data1()[,input$col_log]<0,input$col_cluster]))
                mydata[i,"cluster_up"] <- number_up
                mydata[i,"cluster_down"] <- number_down
            }
            
            if (mydata[i,input$col_log]>0){
                mydata[i,"uniqueness"] <- uniqueness(DEG_clusters=number_up,total_clusters=totalClusters())
            }
            
            if (mydata[i,input$col_log]<0){
                mydata[i,"uniqueness"] <- uniqueness(DEG_clusters=number_down,total_clusters=totalClusters())
            }
        }
        mydata
    })


    # "p-value", "adjusted p-value", "log2FC", "uniqueness"
    data3 <- reactive({
        data_temp <- data2()
        if("p-value" %in% input$considered_thresholds){
            data_temp <- data_temp[as.numeric(data_temp[,input$col_p])<=input$p_threshold,]
        }

        if("adjusted p-value" %in% input$considered_thresholds){
            data_temp <- data_temp[as.numeric(data_temp[,input$col_adj_p])<=input$p_adj_threshold,]
        }

        if("log2FC" %in% input$considered_thresholds){
            data_temp <- data_temp[data_temp[,input$col_log]>=input$logFC_threshold | data_temp[,input$col_log]<=(-input$logFC_threshold),]
        }

        if("uniqueness" %in% input$considered_thresholds){
            data_temp <- data_temp[data_temp$uniqueness>=input$uniqueness_threshold,]
        }
        data_temp
    })

    output$DEG_table_filtered_out <- renderDataTable({
        req(input$DEG_table)
        table_of_considered_clusters()[,-c(length(colnames(data3()))-1,length(colnames(data3()))-2)] #print table without cluster_up und cluster_down
    })

    table_of_considered_clusters <- reactive({
        data3()[as.character(data3()[,input$col_cluster])%in%input$considered_clusters,]
    })
    
    upset_input <- reactive({
        temp <- data3()
        upsetinput <- data.frame("gene"=unique(temp[,input$col_gene]),"up_down"=rep("up",length(unique(temp[,input$col_gene]))))
        upsetinput <- rbind(upsetinput,data.frame("gene"=unique(temp[,input$col_gene]),"up_down"=rep("down",length(unique(temp[,input$col_gene])))))
        
        for (i in input$considered_clusters){
            upsetinput <- cbind(upsetinput,rep(0,length(upsetinput[,1])))
        }
        # fill in 1 for every DEG
        colnames(upsetinput) <- c("gene","up_down",input$considered_clusters)
        for (i in 1:length(temp[,1])){
            gene <- temp[i,input$col_gene]
            #cluster_name <- 
            cluster <- temp[i,input$col_cluster]
            if (temp[i,input$col_log]>0){
                upsetinput[upsetinput$gene==gene&upsetinput$up_down=="up",cluster] <- 1
            } else {
                upsetinput[upsetinput$gene==gene&upsetinput$up_down=="down",cluster] <- 1
            }
        }
        upsetinput
    })
    
    upregulated_genes_vector <- reactive({
        index <- c()
        selection <- as.data.frame(upset_input())[upset_input()$up_down=="up",c("gene","up_down",input$intersection_clusters)]
        for (i in 1:length(selection[,1])){
            if(all(selection[i,3:length(selection[1,])]==1)){
                index <- c(index,i)
            }
        }
        selection[index,"gene"]
    })
    
    downregulated_genes_vector <- reactive({
        index <- c()
        selection <- as.data.frame(upset_input())[upset_input()$up_down=="down",c("gene","up_down",input$intersection_clusters)]
        for (i in 1:length(selection[,1])){
            if(all(selection[i,3:length(selection[1,])]==1)){
                index <- c(index,i)
            }
        }
        selection[index,"gene"]
    })
    
    upregulated_genes <- reactive({
        index <- c()
        temp <- data3()[,shown_columns()]
        for (i in 1:length(temp$gene)){
            if(temp$gene[i] %in% upregulated_genes_vector()){
                if(temp[i,6] %in% input$intersection_clusters){
                    index <- c(index,i)
                }
            }
        }
        temp[index,]
    })
    
    downregulated_genes <- reactive({
        index <- c()
        temp <- data3()[,shown_columns()]
        for (i in 1:length(temp$gene)){
            if(temp$gene[i] %in% downregulated_genes_vector()){
                if(temp[i,6] %in% input$intersection_clusters){
                    index <- c(index,i)
                }
            }
        }
        temp[index,]
    })   
    
    shown_columns <- reactive({
        vector <- c("gene")
        if("p-value" %in% input$considered_thresholds){
            vector <- c(vector,colnames(data3()[input$col_p]))
        }
        
        if("adjusted p-value" %in% input$considered_thresholds){
            vector <- c(vector,colnames(data3()[input$col_adj_p]))
        }
        
        if("log2FC" %in% input$considered_thresholds){
            vector <- c(vector,colnames(data3()[input$col_log]))
        }
        
        if("uniqueness" %in% input$considered_thresholds){
            vector <- c(vector,"uniqueness")
        }
        vector <- c(vector,colnames(data3())[input$col_cluster])
        vector
    })
    
     observe({
         message(shown_columns())
    #     message(downregulated_genes_vector())
     })
    
    output$upset <- renderPlot({
        upset(upset_input(),
              nsets = length(input$considered_clusters),
              nintersects = NA, 
              sets = input$considered_clusters,
              keep.order = TRUE,
              #point.size = 2.8,
              text.scale=1.5,
              #line.size = 1,
              mb.ratio = c(0.6, 0.4))
    })
    
    output$DEG_table_clusters <- renderDataTable({
        req(input$DEG_table)
        #table_of_considered_clusters()
        upset_input()
    })
    
    output$DEGs_intersection_up <- renderDataTable(
        #req(input$DEG_table)
        upregulated_genes()
        #upset_input()
    )

    
    output$DEGs_intersection_down <- renderDataTable(
        #req(input$DEG_table)
        downregulated_genes()
    )
    
    output$download_filtered_table <- downloadHandler(
        filename = function() {
            paste0("filtered_table", ".xlsx")
        },
        content = function(file) {
            myoutput <- openxlsx::createWorkbook()
            openxlsx::addWorksheet(wb=myoutput, sheetName = "Filtered Table")
            openxlsx::writeData(wb=myoutput, sheet = "Filtered Table",x=table_of_considered_clusters())
            openxlsx::saveWorkbook(wb=myoutput,file = file, overwrite = T)
            #vroom::vroom_write(data(), file)
        }
    )
    
    # output$test <- renderText({
    #     print(totalClusters())
    # })
    # 
    # output$test2 <- renderDataTable({
    #     req(input$DEG_table)
    #     data3()
    # })    
})
