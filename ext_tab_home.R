# Load metadata for datasets
tdatasets = fread("dummy_data/ext_mapping_dataset_layer_atherosclerosis.txt", header = TRUE)
tdatasets = data.frame(tdatasets)

# Sort metadata for Home tabe
tabText1 <- paste0(tdatasets$dataset_name[1], "\n  Source: ", tdatasets$relevant_PMIDs_links[1] )

tabText2 <- paste0(tdatasets$dataset_name[2], "\n Source: ", tdatasets$relevant_PMIDs_links[2] )

tabText3 <- paste0(tdatasets$dataset_name[3], "\n Source: ", tdatasets$relevant_PMIDs_links[3] )

tabText4 <- paste0(tdatasets$dataset_name[4] )

# Extract and format the names of layers
tdatasets$evidence_layer_name = str_to_title(string = gsub(x = tdatasets$evidence_layer, pattern = "_", replacement = " "))
levidence.layers = unique(tdatasets$evidence_layer)
levidence.layers.names = unique(tdatasets$evidence_layer_name)
name.total.layer = "Overall"

# Declare the homa tab
home_tab <- argonTabItem(tabName ="tab_home",
             
               argonDashHeader(
                 gradient = TRUE,
                 color = "purple",
                 separator = FALSE,
                 separator_color = "info",
                 top_padding = 4,
                 bottom_padding = 4,          
                           
                 argonCard(
                   title = div(strong("Target Credibility Indicator"), style = "font-size: 300%; text-transform: none; color: black"),
                   hover_lift = FALSE,
                   shadow = TRUE,
                   width = 12,
                   #src = "",
                   #btn_text = "More",
                   shadow_size = NULL,
                   hover_shadow = FALSE,
                   border_level = 4,
                   icon = NULL,
                   status = "primary",
                   background_color = NULL,
                   gradient = TRUE, 
                   floating = FALSE,
                   div("Rank your gene list by combining layers of evidence in atherosclerosis. This tool allows you to aggregate pre-curated biological data and identify the most credible targets in your list. ")
                     
                 ),
               
               ),
               
               argonDashHeader(
                 gradient = TRUE,
                 color = "info",
                 separator = TRUE,
                 separator_color = "secondary",
                 top_padding = 5,
                 bottom_padding = 5,
               
                 argonRow(
                   
                   # Horizontal Tabset
                   argonColumn(
                     width = 12,
                     argonH1("Biological layers of gene-disease association", display = 4),
                     h3("We have downloaded datasets of relevance for atherosclerosis, curated the data and assigned biological layers:"),
                     argonTabSet(
                       id = "tab-1",
                       card_wrapper = TRUE,
                       horizontal = TRUE,
                       circle = FALSE,
                       size = "sm",
                       width = 15,
                       iconList = lapply(X = 1:4, FUN = argonIcon, name = "atom"),
                       argonTab(
                         tabName = tdatasets$evidence_layer_name[1],
                         active = TRUE,
                         tabText1
                       ),
                       argonTab(
                         tabName = tdatasets$evidence_layer_name[2] ,
                         active = FALSE,
                         tabText2
                       ),
                       argonTab(
                         tabName = tdatasets$evidence_layer_name[3],
                         active = FALSE,
                         tabText3
                       ),
                       argonTab(
                         tabName = tdatasets$evidence_layer_name[4],
                         active = FALSE,
                         tabText4
                       )
                     ),
                     br(), 
                     
                     h3("New layers coming soon! including differential expression data from patient samples."),
                     br(),
                     br(),
                     #argonH1("Weigthing the data layers", display = 4),
                     # "Set the sliders to weight each layer and compute scores:",
                     
                   )
                 ),
               ),
               
               argonDashHeader(
                 gradient = TRUE,
                 color = "secondary",
                 separator = TRUE,
                 separator_color = "warning",
                 top_padding = 6,
                 bottom_padding = 2,
                 width = 12,
                 argonRow(
                   argonColumn(
                     width = 12,
                     h1("The algorithm for data integration"),
                     h3("Our Basic algorithm has the purpose of scoring each biological layer with user-defined weights.  Internally, the algorithm uses data transformation functions specific to each data type: the raw data is transformed into the same interval of values [0-1] and later aggregated to score every gene independently."),  
                     br(),br(),
                     div(img(tags$figure(
                       tags$img(
                         src = "basic_algorithm_flowchart.png",
                         width = 600
                       )
                     )), style="text-align: center;")
                   ) 
                  )
               ),
               argonDashHeader(
                 gradient = TRUE,
                 separator = TRUE,
                 color = "warning",
                 top_padding = 5,
                 bottom_padding = 4,
                 separator_color = "success",
                 width = 12,
                 argonRow(
                   argonColumn(
                     width = 12,
                     h1("Quantitative metrics for target prioritisation"),
                     h3("The ranking is interpreted together with the amount of data available for each gene (i.e., data completeness). Another metric for interpreting the ranking is the probability of observing a higher overall score by random chance."),
                     br(),br(),
                     tags$figure(
                       align = "center",
                       tags$img(
                         src = "quantitative_metrics.png",
                         width = 800
                       )
                     )
                   )
                 )
               ),
               argonDashHeader(
                 gradient = TRUE,
                 separator = TRUE,
                 color = "success",
                 top_padding = 5,
                 bottom_padding = 5,
                 width = 12,
                 separator_color = "secondary",
                 argonRow(
                   argonColumn(
                     width = 12,
                     h1("Visualisations bringing interpretatibility to prioritisations results"),
                     h3("Compare your gene scores with the background of all genes in our datasets. Also, pick a gene and check the datasets where it shows stronger evidence. "),
                     br(),br(),
                     div( tags$img(
                         src = "visualisations_home_example.png",
                         width = 800
                       ), style="text-align: center;"),
                     br(),br()
                     
                   )
                 )
               ),
               
               argonDashHeader(
                 gradient = TRUE,
                 separator = TRUE,
                 color = "gray200",
                 top_padding = 5,
                 bottom_padding = 5,
                 width = 12,
                 separator_color = "primary",
                 argonRow(
                   argonColumn(
                     width = 12,
                     h3("We hope you would find this tool useful for reducing long lists of genes, taking promising candidates to your validation assays or complementing the package of information around your favourite targets!")
                     
                   )
                 )
               )
               
               
               
               
                         
                  
               
               
                      
             # argonRow(
             #   argonCard(
             #     title = div(actionLink("tab2link", "Prioritise targets using pre-computed scores for atherosclerosis disease"), style = "font-size: 200%; text-transform: none"), 
             #     #src = "",
             #     #btn_text = "Load",
             #     description = "", 
             #     icon = argonIcon("planet"), 
             #     icon_background = "danger",
             #     hover_lift = TRUE,
             #     background_color = "white",
             #     width = 5,
             #     "This release includes 4 biological layers of evidence in atherosclerosis and +16,000 genes." 
             #   )
               # argonCard(
               #   title = div("Prioritise targets using my own data", style = "font-size: 200%; text-transform: none"), 
               #   #src = "",
               #   #btn_text = "Ready!",
               #   icon = icon("users"), 
               #   icon_background = "warning",
               #   background_color = "Info",
               #   shadow = TRUE,
               #   width = 3,
               #   "Coming soon! Load your own datasets to compute rankings for your gene lists."
               # )
             #)
             
            
) # End of argonTabItem

