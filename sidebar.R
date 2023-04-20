argonSidebar <- argonDashSidebar(

    vertical = TRUE,
    skin = "light",
    background = "white",
    size = "md",
    side = "left",
    id = "my_sidebar",
    #brand_url = "http://www.google.com",
    brand_logo = "https://i.ibb.co/dP44c72/tci-logo3.png",
    argonSidebarHeader(title = "Main Menu"),
    argonSidebarMenu(
      
      
      argonSidebarItem(
        tabName = "tab_home",
        icon = argonIcon(name = "compass-04", color = "black"),
        "Home"
      ),
      argonSidebarItem(
        tabName = "tab_ranking_gene_list",
        icon = argonIcon(name = "bullet-list-67", color = "purple"),
        "Ranking my gene list"
      ),
      argonSidebarItem(
        tabName = "tab_visualisations",
        icon = argonIcon(name = "spaceship", color = "black"),
        "Visualisations"
      ),
      argonSidebarItem(
        tabName = "tab_about",
        icon = argonIcon(name = "compass-04", color = "black"),
        "Aknowledgements"
      )
     
     
    ),
    argonSidebarDivider(),
    br(),br(),br(),br(),br(), 
    br(),br(),br(),br(),br(), br(), 

    
    div( tags$img(
      src = "novo_nordisk_logo.png",
      height="65%", width="65%", align="center"
      
    ),  style="width: 20px", align="center"),
    
    br()
    
    
    
)
