rm(list = ls())
# setwd("scripts/external_app/gallery")

# Base libraries
library(shiny)
library(argonR)
library(argonDash)
library(magrittr)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(stringr) # To capitalise the first character of each word
library(tidyverse)
library(DT) # For rendering data tables

# Libraries for visualisations
library("ggpubr") # For adding p-values to the plots
library("rstatix") # For adding p-values to the plots
library("ggplot2")
library("lattice")
library("gghalves") # Package for half-violin plots
library("gtable")
library("gridExtra")
library("grid")

# template
# browser()
source("sidebar.R")
source("footer.R")

# Tab elements in the app
source("ext_tab_about.R")
source("ext_tab_home.R")
source("ext_tab_ranking_gene_list.R")
source("ext_tab_visualisations.R")


# Calling source of functions
source("ext_functions_basic.R")



# App
shiny::shinyApp(
  ui = argonDashPage(
    title = "TCI",
    author = "Novo Nordisk A/S. This is tool was created in the area of Systems Biology and Target Discovery for Digital Science & Innovation.",
    description = "TCI",
    sidebar = argonSidebar,
    header =  argonDashHeader(
      
                    titlePanel(div("", style = "color: white")),
                    gradient = TRUE,
                    color = "purple",
                    separator = FALSE,
                    separator_color = "white",
                    bottom_padding = 1,
                    top_padding = 6,
                    background_img = NULL,
                    mask = TRUE,
                    opacity = 1,
                    height = 1000
                    
                  
      
    ),
    body = argonDashBody(
      
              argonTabItems(
                home_tab,
                ranking_gene_list_tab,
                visualisations_tab,
                about_tab
              ),
              
              # shiny::tags$script( "$(document).on('click', function(event) {
              #         Shiny.onInputChange('activeTab', $('.active').data().value);});"),
              
              # This routine activates another tab after the user click a link
              shiny::tags$head(
                shiny::tags$script(
                  'Shiny.addCustomMessageHandler("update-tabs_ranking_gene_list", function(message) {
                         
                        // hide and inactivate all not selected tabs
                        $(".active.show").removeClass("active show");
                        $(".tab-pane.active.show").removeClass("active show");
                        
                        // add active class to the current selected tab and show its content
                        $("#tab-tab_ranking_gene_list").addClass("active show");
                        $("#shiny-tab-tab_ranking_gene_list").addClass("active show");
                       });
                  
                    '
                  
                )
              )
      
    ),
    footer = argonFooter
  ),
  
  
  server = function(input, output, session) {
   
   
    
    # Move to tab of data upload
    observeEvent(input$tab2link, {
      
      session$sendCustomMessage(
        type = "update-tabs_ranking_gene_list",
        message = "" 
      )
    })
    
    
    # Select weight for the layers
    # These are 4 observe events for selecting weights for the layers
    observe({
      val1 <- input$layer1
      # Control the value, min, max, and step.
      updateSliderInput(session, "layer2", value = 0,
                        min = 0, max = 100-val1, step = 1)

    })

    observe({
      val1 <- input$layer1
      val2 <- input$layer2
      # Control the value, min, max, and step.
      updateSliderInput(session, "layer3", value = 0,
                        min = 0, max = 100-val1-val2, step = 1)

    })

    observe({
      val1 <- input$layer1
      val2 <- input$layer2
      val3 <- input$layer3
      # Control the value, min, max, and step.
      updateSliderInput(session, "layer4", value = 0,
                        min = 100-val1-val2-val3, max = 100-val1-val2-val3, step = 1)

    })


    # Computing scores  ----
    # Use a reactive function to compute scores using the weights per layer and pre-computed objects 
    compute_scores.function <- eventReactive( 
      
            eventExpr = input$buttom_compute_scores,  
            valueExpr = ({
        
        
        req(input$layer1, input$layer2, input$layer3, input$layer4)
        
        
        all.ttci.score = NULL
        
        # Load the pre-transformed data
        load(file = "dummy_data/transformed_evidence_matrix.RData")
        head(evidence.trans.matrix)
        
        # Load the metadata for datasets
        load(file = "dummy_data/table_datasets.RData")
        tdatasets$weights_total = NULL
        head(tdatasets)
        
        # Load the pre-computed completeness scores
        tconfidence.score = fread("dummy_data/table_confidence_score_atherosclerosis.txt", header = T)
        tconfidence.score = data.frame(tconfidence.score)
        head(tconfidence.score)
        
        # List the evidence layers
        levidence.layers = unique(tdatasets$evidence_layer)
        
        # Read the weight provided by the user
        lweigths.layer = c(input$layer1, input$layer2, input$layer3, input$layer4)/100
        cat(" ++ Weights selected by the user:", lweigths.layer, "\n")
        
        
        
        # Define number of permutations for computing the p-value
        n.permutations = 100
        
        
        # Initialise table of tci per layer
        ttci.score = NULL
        
        # Compute tci evidence score per evidence layer
        cat("-------------------\n Compute tci evidence score per evidence layer. \n------------------- \n")
        
        # Loop across evidence layers
        for(i in 1:length(levidence.layers)){ # i = 1
          
          # Call tci function for the datasets in this layer
          cat(" ++ Calling tci function for the datasets in layer:", levidence.layers[i], "\n")
          
          
          # Sum the scores across data sets per layer and compute the layer score (0-1) and total score.
          res.basic.tci.score = basic_tci_score( evidence.trans.matrix = evidence.trans.matrix[  , tdatasets$evidence_layer == levidence.layers[i], drop = FALSE  ] ,
                                                 basic.weights = tdatasets$weights [ tdatasets$evidence_layer == levidence.layers[i] ] ,
                                                 main = paste("TCI, layer: ", levidence.layers[i], sep = "" ),
                                                 replace.NA.0 = TRUE, 
                                                 p.value.permutation = FALSE, 
                                                 n.permutations = NULL )
          
          # Add score vector to table of tci                      
          ttci.score = cbind(  ttci.score,  
                               res.basic.tci.score$gene.tci.score )
          
          
          
        } # End loop of tci evidence score per evidence layer
        
        
        # Compute the total tci
        cat("-------------------\n Compute the total tci. \n------------------- \n")
        res.basic.tci.score = basic_tci_score( evidence.trans.matrix = ttci.score,
                                               basic.weights = lweigths.layer, # lweigths.layer = rep(0.25, 4)
                                               main = "Total TCI", 
                                               replace.NA.0 = TRUE,  
                                               p.value.permutation = FALSE, # The probabilities are computed below using the whole evidence matrix
                                               n.permutations = NULL )
        
        
        
        # Add the total tci to the table of tci
        ttci.score = cbind(rownames(evidence.trans.matrix), 
                           ttci.score, 
                           res.basic.tci.score$gene.tci.score)
        
        
        # Add colnames to the table and transform to data frame
        colnames(ttci.score) <- c("ensembl_gene_id", paste("tci_", levidence.layers, sep = ""), "tci_total")
        ttci.score = data.frame(ttci.score, stringsAsFactors = F)
        print( head(ttci.score) )  
        rownames(ttci.score) <- rownames(evidence.trans.matrix)
        
        # Prepare table of p.values of tci
        ttci.probability = rownames(ttci.score)
        
        # Compute the P(S>s) for each layer
        for( i in 1:length(levidence.layers)){
          
          cat(" ++ Computing P(S>s) for layer", levidence.layers[i], ". \n")
          
          
          # Use the evidence matrix to permute the scores (only columns of the layer)
          pos.layer = colnames(evidence.trans.matrix) %in% tdatasets$curated_dataset_id[tdatasets$evidence_layer == levidence.layers[i]] 
          res.basic.tci.score = basic_tci_score( evidence.trans.matrix = evidence.trans.matrix[, pos.layer, drop = F],
                                                 basic.weights = tdatasets$weights [ !is.na(tdatasets$weights) & pos.layer ] * lweigths.layer[i],
                                                 main = levidence.layers[i],
                                                 replace.NA.0 = TRUE, 
                                                 p.value.permutation = TRUE, 
                                                 n.permutations = n.permutations,
                                                 permute.option = "rows per column",
                                                 # If permute option is "rows per layer", this argument indicates the columns in each layer.
                                                 permute.col.groups = levidence.layers[i] )
          
          
          
          # Add the probability vector to the table of tci
          ttci.probability = cbind(ttci.probability, 
                                   res.basic.tci.score$gene.tci.probability)
          
          
          
        }
        
        # Compute the P(S>s) to the total tci
        cat(" ++ Computing P(S>s) for total tci. \n")
        
        # Use the whole evidence matrix to permute the scores
        res.basic.tci.score = basic_tci_score( evidence.trans.matrix = evidence.trans.matrix,
                                               basic.weights = tdatasets$weights [ !is.na(tdatasets$weights) ] * lweigths.layer,
                                               main = "Total TCI",
                                               replace.NA.0 = TRUE, 
                                               p.value.permutation = TRUE, 
                                               n.permutations = n.permutations,
                                               permute.option = "rows per column",
                                               # If permute option is "rows per layer", this argument indicates the columns in each layer.
                                               permute.col.groups = tdatasets$evidence_layer [ !is.na(tdatasets$weights_total) ] )
        
        
        # Add the probability vector to the table of tci
        ttci.probability = cbind(ttci.probability, 
                                 res.basic.tci.score$gene.tci.probability)
        
        # Add colnames to the table and transform to data frame
        colnames(ttci.probability) <- c("ensembl_gene_id",  paste("probability_tci_", levidence.layers, sep = ""), "probability_tci_total") # paste("probability_tci_", levidence.layers, sep = "")
        ttci.probability = data.frame(ttci.probability, stringsAsFactors = F)
        print( head(ttci.probability) )
        
        
        
        # Merge the tables of quantitative metrics
        cat("-------------------\n Mergin tables of quantitative metrics. \n------------------- \n")
        
        # Add together the tables of scores and the gene symbols
        all.ttci.score = full_join(ttci.score, tconfidence.score,  by = "ensembl_gene_id")
        all.ttci.score = full_join(all.ttci.score, ttci.probability,  by = "ensembl_gene_id")
        
        # Add column of total tci ranking
        all.ttci.score$ranking_tci_total <- rank( -as.numeric(all.ttci.score$tci_total), na.last = "keep", ties.method = "first")
        all.ttci.score = all.ttci.score[ order(all.ttci.score$ranking_tci_total, decreasing = F),]
        
        # Read the annotations of: Gene symbol
        tgenes.symbol = read.table("id_maps/ensembl_gene_symbol_map.txt", header = T)
        colnames(tgenes.symbol) <-  c("ensembl_gene_id", "selected_symbol")
        tgenes.symbol = data.frame(tgenes.symbol, stringsAsFactors = F)
        
        # Filter duplicated ensembl_gene_id per Gene symbol
        tgenes.symbol = tgenes.symbol[ !duplicated(tgenes.symbol$ensembl_gene_id), ]
        all.ttci.score = right_join(tgenes.symbol, all.ttci.score , by = "ensembl_gene_id")
        all.ttci.score = data.frame(all.ttci.score)
        # all.ttci.score <<- all.ttci.score
        
        
        cat( " \n ++ Computation of scores completed! \n ------------------- \n")
        
        return( all.ttci.score  )
        
        
      }), # End reactive valueExpr
            ignoreNULL =  TRUE # The scores are computed with default weights and shown when the app initialises
      
    ) # End eventReactive
    
    
    
    # Filtering scores by gene list ----
    filtered_scores.function <- eventReactive( 
      req(input$buttom_compute_scores), 
      valueExpr = ({
        
        # Call the function of query genes
        genes.query.0 = genes_query.function()
        all.ttci.score = compute_scores.function()
        
        # If the object of genes.query.0 has been created (by uploading a file)
        if( !is.null(genes.query.0) ){
          print(genes.query.0)
          all.ttci.score.filt = all.ttci.score[ all.ttci.score$ensembl_gene_id %in% genes.query.0, ] 
          
          # Otherwise show all genes in all datasets
        }else{
          all.ttci.score.filt = all.ttci.score # [1:50,]  
        }
        
        # Add column of gene list ranking
        all.ttci.score.filt$ranking <- rank( -as.numeric(all.ttci.score.filt$tci_total), na.last = "keep", ties.method = "first")
        
        cat( " \n ++ Filtering of scores completed! \n ------------------- \n")
        
        return( all.ttci.score.filt  )
        
      }),
      ignoreNULL =  TRUE 
    ) # End eventReactive 
    
    
    
    
    
    
    # Select columns to print ----
    select_filtered_scores.function <- reactive({
      
      req(input$select_quant_metrics)
      
      levidence.layers.total = c(levidence.layers, "total")
      selected.columns = c("ranking","ensembl_gene_id", "selected_symbol", paste(rep(input$select_quant_metrics, each = length(levidence.layers.total)), levidence.layers.total, sep = "_") )
      
      # Select columns to show
      df_sel <- filtered_scores.function() %>% select(selected.columns)
      
      # Sort columns by ranking
      
      # Rename columns
      colnames(df_sel) <- gsub(x = colnames(df_sel), pattern = "confidence_score", replacement = "compl.")
      colnames(df_sel) <- gsub(x = colnames(df_sel), pattern = "selected_symbol", replacement = "gene_symbol")
      colnames(df_sel) <- gsub(x = colnames(df_sel), pattern = "probability", replacement = "prob.")
      
      # Remove spaces and capitalise words
      colnames(df_sel) <- str_to_title(string = gsub(x = colnames(df_sel), pattern = "_", replacement = " "))
      colnames(df_sel) <- gsub(x = colnames(df_sel), pattern = "Tci", replacement = "Gene score")
      colnames(df_sel) <- gsub(x = colnames(df_sel), pattern = "In Vitro", replacement = "In-Vitro")
      
      return(df_sel)
      
    })
    
    
    
    # Send the table to the ui through a render table function.
    output$query_computed_scores <- renderDT({
      

      formated.scores.table <- DT::datatable( data =  as.data.frame( select_filtered_scores.function() ),
                                              rownames= FALSE,  
                                              options = list( order = list( list(0, 'asc') ),
                                                              pageLength = 5
                                                            )
                                              ) 
      #browser()
      # https://towardsdatascience.com/top-7-packages-for-making-beautiful-tables-in-r-7683d054e541
      # Give format to columns of TCI
      cols.tci = grep( pattern = "Gene score", x = colnames(formated.scores.table$x$data) ) 
      if( length(cols.tci) > 0 ){
        
        formated.scores.table <- formated.scores.table %>% 
                                  formatStyle(cols.tci,
                                              background = styleColorBar( seq(0,1,0.01), '#80ed99') ,
                                              backgroundSize = '98% 88%',
                                              backgroundRepeat = 'no-repeat',
                                              backgroundPosition = 'center'
                                  )
        
      }
      
      # Give format to columns of completeness
      cols.tci = grep( pattern = "Compl", x = colnames(formated.scores.table$x$data) )
      if( length(cols.tci) > 0 ){
        
        formated.scores.table <- formated.scores.table %>% 
          formatStyle(cols.tci,
                      background = styleColorBar( seq(0,1,0.01), '#e9c46a'),
                      backgroundSize = '98% 88%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          )
        
      }
      
      # Give format to columns of completeness
      cols.tci = grep( pattern = "Prob", x = colnames(formated.scores.table$x$data) ) 
      if( length(cols.tci) > 0 ){
        
        formated.scores.table <- formated.scores.table %>% 
          formatStyle(cols.tci,
                      background = styleColorBar( seq(0,1,0.01), '#f28482'),
                      backgroundSize = '98% 88%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'right'
          )
        
      }
      
      
      # Give format to columns of TCI total
      cols.tci = grepl( pattern = "TCI Total", x = colnames(formated.scores.table$x$data) ) & 
        !grepl( pattern = "Prob. TCI Total", x = colnames(formated.scores.table$x$data) )
      if( length(cols.tci) > 0 ){
        
        formated.scores.table <- formated.scores.table %>% 
          formatStyle(cols.tci,
                      background = styleColorBar( seq(0,1,0.01), '#48cae4'),
                      backgroundSize = '98% 88%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'left'
          )
        
      }
      
      return(formated.scores.table)
      
    }) # End renderDT

     
    # Downloadable tsv of filtered scores table ----
    output$buttom_download_scores <- downloadHandler(
      
      filename = function() {
        "tci_ranking_gene_list.tsv"
      },  
      content = function(file) {
        write.table(x = select_filtered_scores.function(), file, row.names = FALSE, sep = "\t", quote = FALSE )
      }
      
    )
    
    
    
    # Query the scores for this gene list when the user clicks the button ----
    genes_query.function <- eventReactive( 
      
          eventExpr = input$file_gene_list, 
          valueExpr =  ({
             
              # Read the gene list queried by the user
              if( !is.null( input$file_gene_list$datapath ) ){
                
                genes.query.0 <- read.table(file = input$file_gene_list$datapath,
                                            header = FALSE,
                                            sep = ",")
                
                genes.query.0 <- c(genes.query.0$V1)
                
              }else{
                
                genes.query.0 = NULL
                
              }
              
        
              # Filter table to only genes in the query of the user
              return( genes.query.0 )
      
      }), ignoreNULL = FALSE # End of valueExpr
    ) # End of eventReactive


     # Adapt the data for violin plots visualisation 
    violin_plot_table.function <- eventReactive( 
      
      eventExpr = c(input$buttom_compute_scores, input$layer_violin_plot), 
      valueExpr =  ({
        
        # Call the function of query genes and compute scores
        genes.query.0 = genes_query.function()
        ttci.score.background = compute_scores.function()
        ttci.score.background = ttci.score.background[, !( colnames(ttci.score.background) %in% "ranking_tci_total" ) ]
        
        # Adapt the input table to label genes either as gene list of background
        ttci.score.background$my_list = "Background"
        
        # If the object of genes.query.0 has been created (by uploading a file)
        if( !is.null(genes.query.0) ){
          ttci.score.background$my_list [ ttci.score.background$ensembl_gene_id %in% genes.query.0 ] <- "Gene list"
        }
        
        # Find the columns to plot
        sel.cols = grep(pattern = c( "total", levidence.layers) [ c(name.total.layer, levidence.layers.names)  == input$layer_violin_plot ], 
                        x = colnames(ttci.score.background))
        
        # Change to numeric
        for( j in sel.cols){
          ttci.score.background[ , j] <- as.numeric( ttci.score.background[ , j] )
        }
        
        # Select the column from the table
        ttci.score.background = ttci.score.background[, c( sel.cols, which( colnames(ttci.score.background) %in% "my_list" ) )]
       
         
        # Return the table for violinplots
        return(ttci.score.background)
        
      }), ignoreNULL = FALSE # End of valueExpr
        
    )  # End of eventReactive

    
    # Violin plots of quantitative metrics ----
    
    output$violin_plots_metrics <- renderPlot({
      
      
      # Number of permutations to compute the p-value
      bootstrap.size = 100
      
      # Call the function that adapts the data for plots
      ttci.score.background = violin_plot_table.function()
      
      # Select colours for boxplots
      colours = c( RColorBrewer::brewer.pal(n = 6, name = "Pastel1")  [c(2,1,3)], 
                   RColorBrewer::brewer.pal(n = 12, name = "Set3")  [c(12)]
      )
      colour.labels = rep("black", nrow(ttci.score.background)) 
      
      
      # Plot of overall tci vs background
      p <- ggplot( ttci.score.background, 
                   aes(x = my_list, y = ttci.score.background[,1])) + ylim(0,1.05) +
        theme_classic(base_size = 18) + labs(y="Gene score", x = "")
      
      if( input$layer_violin_plot == name.total.layer){
        p1 <- p + geom_violin(trim=T, show.legend = FALSE, color = colours[1], fill = colours[1]) + geom_boxplot(width=0.1)  
      }else{
        p1 <- p + geom_violin(trim=T, show.legend = FALSE, color = colours[3], fill = colours[3]) + geom_boxplot(width=0.1)
      }
      
      
      
      # Compute p-value for median differences between the gene list and the background
      if(  length(unique(ttci.score.background$my_list)) > 1 ){
        
        mygenelist.A = ttci.score.background[,1][ ttci.score.background$my_list != "Background"]
        mygenelist.B = ttci.score.background[,1][ ttci.score.background$my_list == "Background"]
        
        
        p.value = bootstrapping.median.differences( sample.A = mygenelist.A, 
                                                    sample.B = mygenelist.B, 
                                                    bootstrap.size = bootstrap.size, 
                                                    alternative = "greater")[[1]]
        
        # Add format to the empirical p-value
        if( p.value == 1/bootstrap.size ){ # If reached minimum (given the number of boostrap samplings)
          p.value = paste0("p<", p.value, ", Median differences")
        }else{
          p.value = paste0("p=", p.value, ", Median differences")
        }
        
        # Create tibble for the ggplot function that adds p-values to your plot
        stat.test = tibble( group1 = levels( as.factor(ttci.score.background$my_list) )[1] ,
                            group2 = levels( as.factor(ttci.score.background$my_list) ) [2],
                            p.value = p.value)
        
        p1 <- p1 + stat_pvalue_manual(
          stat.test, y.position = 1.05, label = "p.value", tip.length = 0.01
        )
      
      } # End adding median differences p-values
      
      # Plot of completeness
      p <- ggplot( ttci.score.background, 
                   aes(x = my_list, y = ttci.score.background[,2])) + ylim(0,1.05) +
        theme_classic(base_size = 18) + labs(y="Data completeness", x = "")
      
      
      p2 <- p + geom_violin(trim=T, show.legend = FALSE, color = colours[4], fill = colours[4]) + geom_boxplot(width=0.1)
      
      
      # Compute p-value for median differences between the gene list and the background
      if(  length(unique(ttci.score.background$my_list)) > 1 ){   
        mygenelist.A = ttci.score.background[,2][ ttci.score.background$my_list != "Background"]
        mygenelist.B = ttci.score.background[,2][ ttci.score.background$my_list == "Background"]
        
        
        p.value = bootstrapping.median.differences( sample.A = mygenelist.A, 
                                                    sample.B = mygenelist.B, 
                                                    bootstrap.size = bootstrap.size, 
                                                    alternative = "greater")[[1]]
        
        # Add format to the empirical p-value
        if( p.value == 1/bootstrap.size ){ # If reached minimum (given the number of boostrap samplings)
          p.value = paste0("p<", p.value, ", Median differences")
        }else{
          p.value = paste0("p=", p.value, ", Median differences")
        }
        
        # Create tibble for the ggplot function that adds p-values to your plot
        stat.test = tibble( group1 = levels( as.factor(ttci.score.background$my_list) )[1] ,
                            group2 = levels( as.factor(ttci.score.background$my_list) ) [2],
                            p.value = p.value)
        
        p2 <- p2 + stat_pvalue_manual(
          stat.test, y.position = 1.05, label = "p.value", tip.length = 0.01
        )
      
      } # End adding median differences p-values
         
      
      # Plot of p-value
      p <- ggplot( ttci.score.background, 
                   aes(x = my_list, y = ttci.score.background[,3])) + ylim(0,1.05) +
        theme_classic(base_size = 18) + labs(y="Probability of higher gene score", x = "")
      
      
      p3 <- p + geom_violin(trim=T, show.legend = FALSE, color = colours[2], fill = colours[2]) + geom_boxplot(width=0.1)
      
      # Compute p-value for median differences between the gene list and the background
      if(  length(unique(ttci.score.background$my_list)) > 1 ){    
        mygenelist.A = ttci.score.background[,3][ ttci.score.background$my_list != "Background"]
        mygenelist.B = ttci.score.background[,3][ ttci.score.background$my_list == "Background"]
        
        p.value = bootstrapping.median.differences( sample.A = mygenelist.A, 
                                                    sample.B = mygenelist.B, 
                                                    bootstrap.size = bootstrap.size, 
                                                    alternative = "less")[[1]]
        
        # Add format to the empirical p-value
        if( p.value == 1/bootstrap.size ){ # If reached minimum (given the number of boostrap samplings)
          p.value = paste0("p<", p.value, ", Median differences")
        }else{
          p.value = paste0("p=", p.value, ", Median differences")
        }
        
        # Create tibble for the ggplot function that adds p-values to your plot
        stat.test = tibble( group1 = levels( as.factor(ttci.score.background$my_list) )[1] ,
                            group2 = levels( as.factor(ttci.score.background$my_list) ) [2],
                            p.value = p.value)
        
        p3 <- p3 + stat_pvalue_manual(
          stat.test, y.position = 1.05, label = "p.value", tip.length = 0.01
        )
     
      } # End adding median differences p-values
      
      # Plot the grid
      p.all = gridExtra::grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)
      return( print(p.all) )
      
    }) # End of violin plots
    
    
    # This event updates the select input of for density plots using the background of genes
    # Update gene list for density plots ----
    observe({ 
      
      req(input$buttom_compute_scores, input$layer_violin_plot) 
      lgene.choices = compute_scores.function()$selected_symbol
      updateSelectInput(session = session, inputId = "gene_density_plot",
                        label = "Select the gene:",
                        choices = lgene.choices, 
                        selected = "APOB")
      
    })
    
    
    # Density plots ----
    output$density_plots <- renderPlot({
    
      if( input$gene_density_plot == "" ){
        
        return()
        
      }else{
      
        
        # Assign the table of scores
        ttci.score <- compute_scores.function()
        
        # Load object with the transformed data.
        load(file = "dummy_data/transformed_datasets.RData"  ) # load(file = "scripts/external_app/gallery/dummy_data/transformed_datasets.RData"  )
        
        #browser()
        
        # Find the position of the gene in the table
        pos.genes = match( input$gene_density_plot, ttci.score$selected_symbol )  
        ( lselected.gene = ttci.score$ensembl_gene_id [ pos.genes ])
        
        # Find the layer selected by the user
        select.layer.density = levidence.layers [ levidence.layers.names  == input$layer_density_plot ]
        
        # Filter the input tables
        ltrans.datasets = ltrans.datasets[ which( names(ltrans.datasets) %in% tdatasets$curated_dataset_id[ tdatasets$evidence_layer == select.layer.density ] ) ]
        
        # Onwards: Only one gene is accepted and its position in lselected.gene is given by
        g = 1
        
        
        # Select the data mode (layers or transformed data)
        ldata.mode = input$select_data_mode_density_plot
        
        for( m in 1:length(ldata.mode)){ # m = 1
          
          # Iniatialise list of datasets 
          ldatasets = NULL
          
          # Select data mode 
          data.mode = ldata.mode[m]
          cat(" +++ Plotting data mode: ", data.mode, "\n")
          
          
          # Construct list of datasets depending on data.mode:
          if( data.mode == "transformed dataset" ){
            
            
            ldatasets <- ltrans.datasets
            xlim.plot = c(0,1)
            ldatasets[ is.na(ldatasets)] <- 0
            
          }else if( data.mode == "layer" ){
            
            
            # Move the columns to lists
            ltci.layers = paste("tci_", select.layer.density, sep = "")
            k = 1
            for( j in which( colnames(ttci.score) %in% ltci.layers) ){
              
              ldatasets[[k]] = data.frame(ttci.score[,j, drop = F])
              rownames(ldatasets[[k]]) <- as.character(ttci.score$ensembl_gene_id) # To make unique names: make.names(as.character(ttci.score$ensembl_gene_id), unique = TRUE)
              names(ldatasets)[k] <- colnames(ttci.score)[j]
              k = k + 1
              
            }
            xlim.plot = c(0,1)
            
          }
          
          
          # Explore the raw data:
          length( ldatasets )
          names(ldatasets) 
          
          
          # Loop across layers plotting patterns in raw data
          par(mfrow = c(1,1))
          plot_list = list()
          
          # Loop across dataset
          
          for( i in 1:length(ldatasets) ){
            
            # Extract the dataset
            dataset = ldatasets[[i]] # names(ldatasets)
            
            
            # Find name of the layer:
            ( layer.name = tdatasets$evidence_layer[ tdatasets$curated_dataset_id == names(ldatasets)[i] ] )
            
            # Change to data.frame
            dataset = as.data.frame(dataset, stringsAsFactors = F)
            if( !( "ensembl_gene_id" %in% colnames(dataset) ) ){
              dataset$ensembl_gene_id = rownames(dataset)
            }else{
              dataset <- dataset %>% relocate(ensembl_gene_id, .after = last_col())
            }
            
            
            # Change from wide table to long table
            dataset.long <- pivot_longer(dataset, cols = 1:(ncol(dataset)-1), names_to = "data_column", values_to = "value")
            dataset.long$value = as.numeric(dataset.long$value)
            
            
            # Impute missing data for those modes after raw data
            if( data.mode != "raw dataset" ){
              dataset.long$value[ is.na( dataset.long$value ) ] <- 0
              
              # Add the ranking of the gene within each data_column
              if( sum( dataset.long$value == 0, na.rm = T) <= length(dataset.long$value) ){
                
                dataset.long <- dataset.long %>% 
                  group_by(data_column) %>% 
                  # mutate(percentile = rank( -value, ties.method = "max" ) )
                  mutate(percentile = paste0("p", ceiling(rank( -value, ties.method = "max" )/length(value)*100) ))
                # mutate(percentile = ceiling( ( 1-ecdf(value)(value) )*100) )
                
              }
              
            }else{
              
              # dataset.long$percentile = 1
              
            }
            
            
            # Filter the dataset to contain only genes in the gene query
            dataset = dataset[ dataset$ensembl_gene_id %in% ttci.score$ensembl_gene_id, ]
            
            
            # Trim names of columns
            if( sum( nchar(dataset.long$data_column) > 150 ) > 1 ){
              dataset.long$data_column <- paste( gsub( x = substr( x = dataset.long$data_column, start = 1, stop = 17), pattern = "[.]", replacement = "_"),
                                                 "_column", as.numeric(as.factor(dataset.long$data_column)), sep ="" )
              dataset.long$data_column = as.factor(dataset.long$data_column)
              head(dataset.long)
              
            }
            
            # Font size in the plots
            if( ncol(dataset) > 10 ){
              axis.text.y_size = 6
            }else{
              axis.text.y_size = 12
            }
            
            
            # Density plots are actually half Violin plots
            p <- ggplot( dataset.long, 
                         aes(x = data_column, y = value)) +  
              geom_half_violin(trim = T, side = "right", 
                               show.legend = FALSE, color = "deepskyblue1", 
                               fill = "lightblue")  +
              theme_classic(base_size = 14) +  coord_flip() +
              theme( axis.text.y = element_text(size = axis.text.y_size)) + labs(y="Gene score", x = "", title = NULL) 
            
            
            
            # Adapt limits
            if( data.mode != "raw dataset" ){
              p <- p + ylim( xlim.plot) 
              size.point = 8
            }else{
              size.point = 4
            }
            
            
            
            # Create label of title
            if( data.mode != "layer" ){
              
              # Find subtitle with names of datsets
              title.label = stringr::str_wrap(tdatasets$dataset_name [ match( names(ldatasets)[i], tdatasets$curated_dataset_id ) ], width=40)
              
            }else{
              
              title.label = levidence.layers.names [ paste0("tci_", levidence.layers) == layer.name ]
              
            }
            
            # Add label of title
            p <- p + labs(title = title.label) + 
              theme(plot.subtitle = element_text(size = 8))
            
            # Add points for selected genes
            dataset.long.select = dataset.long[ dataset.long$ensembl_gene_id %in% lselected.gene[g], ]
            dataset.long.select$point_colour = "black"
            dataset.long.select$text_colour = "white"
            
            # Highlight in green if the gene value is within the 10th percentile 
            if( data.mode != "raw dataset" ){ 
              dataset.long.select$point_colour [ as.numeric( gsub(dataset.long.select$percentile, pattern = "p", replacement = "") ) <= 10 ] <- "green"
              dataset.long.select$text_colour [ as.numeric( gsub(dataset.long.select$percentile, pattern = "p", replacement = "") ) <= 10 ] <- "black" 
            }
            
            dataset.long.select = arrange(dataset.long.select, data_column)
            
            p <- p + geom_point(data = dataset.long.select,
                                aes(x = as.factor(data_column), y = value, group = ensembl_gene_id), 
                                size = size.point, 
                                color = dataset.long.select$point_colour
            ) 
            
            # For data mode different than raw data, add the percentile to the plot
            # In raw data, the percentile interpretation could be misleading due to the different trends for good/bad evidence
            if( data.mode != "raw dataset" ){
              
              p <- p + geom_text(data = dataset.long.select,
                                 aes(x = as.factor(data_column), y = value, label = percentile, group = ensembl_gene_id), 
                                 size = 2.5,
                                 color = dataset.long.select$text_colour,
                                 alpha = 1 )
              
              
            }
            
            # Add plot to the list
            plot_list[[i]] <- p
            
            
          } # End loop across dataset
          
          # Plot in a grid
          n.ncol = 1
          
        } # End loop across data modes
        
        return( print( do.call( grid.arrange, c( plot_list, ncol = n.ncol) )) )
        
      }# End if gene selected != ""
      
    }) # End of density plots
    
    
  }
)