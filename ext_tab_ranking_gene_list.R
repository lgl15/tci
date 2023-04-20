ranking_gene_list_tab <- argonTabItem(tabName ="tab_ranking_gene_list",
             
                                      
                                    h1("Prioritise targets using evidence layers for atherosclerosis disease"),
                                    h4("This release includes 4 layers and +16,000 genes." ),
                                      
                                    argonRow( 
                                      argonCard(
                                        width = 12,
                                        title = div(strong("Weighting of evidence layers"), style = "font-size: 160%; text-transform: none"),
                                        src = NULL,
                                        hover_lift = F,
                                        shadow = TRUE,
                                        shadow_size = NULL,
                                        hover_shadow = FALSE,
                                        border_level = 0,
                                        icon = argonIcon("ui-04"),
                                        status = "primary",
                                        background_color = NULL,
                                        gradient = FALSE, 
                                        floating = FALSE,
                                        
                                        h4("Set the weights (%) of layers in sequence and the last layer will be fixed to make the total 100%"),
                                        br(),
                                        argonRow(
                                          argonColumn(
                                              sliderInput("layer1", div( paste0("1. ", tdatasets$evidence_layer_name[1]), style = "font-size: 100%; text-transform: none") , min=0, max=100, value=0,
                                                          step=1)
                                          ),
                                          argonColumn(
                                              sliderInput("layer2", paste0("2. ", tdatasets$evidence_layer_name[2]), min=0, max=100, value=0,
                                                          step=1)
                                          ), 
                                          argonColumn(
                                              sliderInput("layer3", paste0("3. ", tdatasets$evidence_layer_name[3]), min=0, max=100, value=0,
                                                          step=1)
                                          ), 
                                          argonColumn(
                                              sliderInput("layer4", paste0("4. ", tdatasets$evidence_layer_name[4]), min=0, max=100, value=0,
                                                          step=1)
                                          ) 
                                        ) # End of argonRow
                                      ), # End of argonCard
                                        
                                      argonCard(
                                        width = 7,
                                        title = div(strong("Upload your list of genes"), style = "font-size: 160%; text-transform: none"),
                                        src = NULL,
                                        hover_lift = F,
                                        shadow = TRUE,
                                        shadow_size = NULL,
                                        hover_shadow = FALSE,
                                        border_level = 0,
                                        icon = argonIcon("bag-17"),
                                        status = "primary",
                                        background_color = NULL,
                                        gradient = FALSE, 
                                        floating = FALSE,
                                        
                                        argonRow(center = FALSE,
                                          argonColumn(center = FALSE,
                                            
                                            h4("Currently accepting a list of ENSEMBL ids, one gene per line. Optionally, just click on compute scores to rank +16,000 genes found in all input datasets."),
                                            fileInput(inputId = "file_gene_list",
                                                      label = "",
                                                      multiple = FALSE,
                                                      accept = c("text/csv",
                                                                 "text/comma-separated-values,text/plain",
                                                                 ".csv", ".txt", ".tsv")),
                                            
                                          ),
                                          
                                          argonColumn(center = TRUE,
                                            
                                            # Compute scores and download files
                                            actionButton(inputId = "buttom_compute_scores", 
                                                       label = "Compute scores", 
                                                       style = "simple", size = "sm",
                                                       color = "royal"),
                                            br(), br(), 
                                            
                                            downloadButton(outputId = "buttom_download_scores", 
                                                           label = "Download scores",
                                                           style = "simple", size = "sm",
                                                           color = "royal")
                                            
                                            
                                          ) # End of argonColumn
                                        ) # End of argonRow
                                      ), # End of argonCard
                                    
                                      
                                      argonCard(
                                        width = 4,
                                        title = div(strong("Change the table view"), style = "font-size: 160%; text-transform: none"),
                                        src = NULL,
                                        hover_lift = F,
                                        shadow = TRUE,
                                        shadow_size = NULL,
                                        hover_shadow = FALSE,
                                        border_level = 0,
                                        icon = argonIcon("settings"),
                                        status = "primary",
                                        background_color = NULL,
                                        gradient = FALSE, 
                                        floating = FALSE,
                                                    
                                        h4("Select quantiative metrics:"),
                                        checkboxGroupInput(inputId = "select_quant_metrics", 
                                                           label = "", 
                                                           choices = c("Gene score" = "tci", "Data completeness" = "confidence_score", "Probability of higher gene score" = "probability_tci" ),
                                                           selected = "tci"
                                        )
                                              
                                      ) # End of argonCard
                                      
                                    ), # End of argonRow
                                      
                                      
             # # Panel for selecting weights
             # sidebarLayout(
             #   
             #   sidebarPanel(
             #   
             #          h2("Weighting of evidence layers:"),
             #          sliderInput("layer1", tdatasets$evidence_layer_name[1], min=0, max=100, value=25,
             #                      step=1),
             #          sliderInput("layer2", tdatasets$evidence_layer_name[2], min=0, max=100, value=25,
             #                      step=1),
             #          sliderInput("layer3", tdatasets$evidence_layer_name[3], min=0, max=100, value=25,
             #                      step=1),
             #          sliderInput("layer4", tdatasets$evidence_layer_name[4], min=0, max=100, value=25,
             #                      step=1),
             #          
             #          width = 6
             #          
             #   ),
             #   
             #   # Main panel for uploading gene lists, comput scores and download files
             #   mainPanel(
             #  
             #      # Input gene list: Select a file ----
             #      h2("Upload your list of genes:"),
             #      h4("Currently accepting a list of ENSEMBL ids, one gene per line"),
             #      fileInput(inputId = "file_gene_list",
             #                label = "",
             #                multiple = FALSE,
             #                accept = c("text/csv",
             #                           "text/comma-separated-values,text/plain",
             #                           ".csv", ".txt", ".tsv")),
             #      
             #      # radioButtons(inputId = "radio_button_gene_list",
             #      #              label = "Select view",
             #      #              choices = c("Show gene list" = "genelist0",
             #      #                          "Show all genes" = "all"
             #      #              ),
             #      #              selected = "all"
             #      # ),
             #      #     
             # 
             #      # Horizontal line to separate panels
             #      hr(),
             #      
             #     # Compute scores and download files
             #     actionBttn(inputId = "buttom_compute_scores", 
             #                label = "Compute scores",
             #                style = "simple", size = "sm",
             #                color = "royal"),
             #     actionBttn(inputId = "buttom_dowload_scores", 
             #                label = "Download scores",
             #                style = "simple", size = "sm",
             #                color = "royal"),
             #     
             #     # Horizontal line to separate panels
             #     hr(),
             #     
             #     h2("Select quantiative metrics:"),
             #     checkboxGroupInput(inputId = "select_quant_metrics", 
             #                        label = "", 
             #                        choices = c("Gene score" = "tci", "Data completeness" = "confidence_score", "Probability of higher gene score" = "probability_tci" ),
             #                        selected = "tci"
             #     ),
             #     
             #     width = 6
             #     
             #   )
             # ), # End side bar layout
             # 
             # 
             #    
             # # Horizontal line to separate panels
             # hr(),

             # Output: Table of scores ----
             DT::dataTableOutput(outputId = "query_computed_scores")
            
            
) # End of argonTabItem

