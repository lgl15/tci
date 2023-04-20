visualisations_tab <- argonTabItem(tabName = "tab_visualisations",
           
           argonRow( 
             argonCard(
               width = 12,
               title = div(strong("Comparing quantitative metrics"), style = "font-size: 190%; text-transform: none"),
               src = NULL,
               hover_lift = F,
               shadow = TRUE,
               shadow_size = NULL,
               hover_shadow = FALSE,
               border_level = 0,
               icon = argonIcon("album-2"),
               status = "primary",
               background_color = NULL,
               gradient = FALSE, 
               floating = FALSE,
               
               h4("After ranking your gene list, here you can compare the quantiative metrics versus the background of +16,000 genes."),
               br(),
               # List to select the layer to plot
               selectInput(inputId = "layer_violin_plot", multiple = FALSE, width = '20%',
                           label = "Select the layer:",
                           selected = name.total.layer,
                           choices =  c(name.total.layer, tdatasets$evidence_layer_name )
               ),
               
               
               # Output: Violin plots ----
               plotOutput(outputId = "violin_plots_metrics")
               
             ) # End of argonCard
           ), # End of argonRow 
             
           argonRow( 
             argonCard(
               width = 12,
               title = div(strong("How a single gene compares with the background?"), style = "font-size: 190%; text-transform: none"),
               src = NULL,
               hover_lift = F,
               shadow = TRUE,
               shadow_size = NULL,
               hover_shadow = FALSE,
               border_level = 0,
               icon = argonIcon("square-pin"),
               status = "primary",
               background_color = NULL,
               gradient = FALSE, 
               floating = FALSE,
               
               h4("Observe the percentile of the distribution for a gene across layers and transformed data sets."),
               br(),
               sidebarLayout(
                 
                 sidebarPanel( width = '75%',
                               # List to select the layer to plot
                               selectInput(inputId = "gene_density_plot", multiple = FALSE, 
                                           label = "Select the gene:", selectize = TRUE, size = NULL,
                                           selected = "APOB",
                                           choices =  "" 
                               ),
                               
                               # List to select the layer to filter the plots
                               selectInput(inputId = "layer_density_plot", multiple = FALSE,
                                           label = "Select the layer:",
                                           selected = tdatasets$evidence_layer_name[1],
                                           choices = tdatasets$evidence_layer_name 
                               ),
                               
                               # Select data mode (layers or transformed data)
                               radioButtons(inputId = "select_data_mode_density_plot", 
                                            label = "Select the level of score aggregation:", 
                                            choices = c("Gene scores in datasets within this layer" = "transformed dataset", "Gene scores in this layer" = "layer"),
                                            selected = "transformed dataset"
                               ),
                               
                               img(
                                 src = "legend_density_plots.png",height = 240, width = 300
                                 
                               )
                               
                               
                 ), # End of sidebarPanel
                 
                 
                 # Show a plot of the generated distribution
                 mainPanel(
                   # Output: Density plots ----
                   plotOutput(outputId = "density_plots")
                 )
                 
               ) # End of sidebarLayout
               
             ) # End of argonCard
           ) # End of argonRow                          
                                   
              
          
             
            
             
) # End of Tab for downloading all data
