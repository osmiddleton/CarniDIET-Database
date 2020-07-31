
# Code for producing a shiny app

# Aim of CarniDIET Shiny app:
# - Interactive map of all study (unique lat/long) locations
#   - Users can subset by family, genus or individual species
#   - families are coloured differently 
#   - When a species is selected, a silhouette of the species is downloaded from Phylopic and placed in the top
#   - Markers will be stock silhouettes of the species in question or the genus (use the Phylopic extraction function)
#
# - Studies can be selected and information extracted from each of the studies
#   - Diet composition
#   - Prey species list
#   - Sample size
#   - Duration of study
#   
# - It will be have an interactive query from that allows the user to download what they want from the dataset

# Packages ----
library(pacman)
pacman::p_load(shiny, leaflet, sf, sp, tidyverse, shinythemes)

# Loading data -----
carnidiet <- read.csv("./Data/data.csv") # CarniDIET
current.ranges <- st_read("./Data", "current.ranges") # Current geographic ranges
pn.ranges <- st_read("./Data", "pn.ranges") # Potentual geographic ranges

# Tidy data
temporal.diet <- carnidiet %>% 
  dplyr::group_by(Family,CarniBinom, Age, Sex,
                  PrimaryRef, Country, Region,
                  SiteName, AltitudeMaximum,
                  Season,Start.Year, End.Year,
                  Longitude.Centroid.Mean.x, Latitude.Centroid.Mean.x,Scat.Stomach.Tissue) %>% 
  dplyr::summarise(x = length(CarniBinom))
# Quick tidy/checks
temporal.diet <- temporal.diet[!is.na(temporal.diet$Latitude.Centroid.Mean.x),]
temporal.diet <- temporal.diet[!is.na(temporal.diet$Longitude.Centroid.Mean.x),]
temporal.diet <- temporal.diet[!temporal.diet$Latitude.Centroid.Mean.x == 0,]
temporal.diet <- temporal.diet[!temporal.diet$Longitude.Centroid.Mean.x == 0,]
str(temporal.diet)

#temporal.diet %>% filter(Family == "Mephitidae") %>% pull(CarniBinom)
tib <- as_tibble(temporal.diet)
colnames(temporal.diet)[13:14] <- c("long","lat")
colnames(tib)[13:14] <- c("long","lat")

# ---- Basic shiny app ----
# # ui.R ----
# ui <- bootstrapPage(
#   
#   titlePanel("CarniDIET: a quantitative database of the diet of the world's carnivorous mammal"),
#     
#   leafletOutput("mymap"),
#     
#   p(),
#     
#   headerPanel("Species selection"),
#     
#   sidebarPanel(
#       
#     uiOutput("select_var1"), 
#     uiOutput("select_var2")
#     
#     ),
#     
#     mainPanel(
#       tableOutput("table")
#     )
#   )
#   
# # server.R ----
# server <- function(input, output, session) {
#     
#     # current.ranges <- reactive({
#     #   ranges %>% filter(Species == input$var2)
#     # })
#     # 
#     tab <- reactive({ 
#       
#       tib %>% 
#         filter(Family == input$var1) %>% 
#         filter(CarniBinom == input$var2)
#       
#     })
#     
#     output$select_var1 <- renderUI({
#       
#       selectizeInput('var1', 'Select family', choices = c("select" = "", levels(tib$Family)))
#       
#     })
#     
#     output$select_var2 <- renderUI({
#       
#       
#       choice_var2 <- reactive({
#         tib %>% 
#           filter(Family == input$var1) %>% 
#           pull(CarniBinom) %>% 
#           as.character()
#         
#       })
#       
#       selectizeInput('var2', 'Select species', choices = c("select" = "", choice_var2())) # <- put the reactive element here
#       
#     })
#     
#     output$table <- renderTable({ 
#       
#       tab()
#       
#     })
# 
#     # Need to have a map in the first place which show ALL data
#     # Also filter through IUCN range maps and add the species' geographic range to the map
#     output$mymap <- renderLeaflet ({
#       leaflet() %>%
#         addTiles() %>%
#         addCircles(data = tib,
#                          lng = ~long,
#                          lat = ~lat,
#                          popup = ~CarniBinom,
#                          color = "lightgrey") %>% 
#         addCircles(data = tab(),
#                          lng = ~long,
#                          lat = ~lat,
#                    color = "black",
#                    weight = 5) 
#         # addPolygons(data = range(),
#         #             color = NULL,
#         #             smoothFactor = 1,
#         #             opacity = 0.5)
#         
#     })
#       
#     }
#     
#     # This apparently can update the map
#         
#     # observe({
#     #   leafletProxy("map", data = tab()) %>%
#     #     clearShapes() %>%
#     #     addCircleMarkers(lat = ~lat,
#     #                      lng = ~long,
#     #                      popup = ~CarniBinom)
#     # })
# 
# # Run the app ----
# shinyApp(ui = ui, server = server)
# 

# ---- 2. Developing shiny app ----
# ui.R ----
ui <- bootstrapPage(
  
  theme = shinytheme("united"),
  
  tabsetPanel(
  
  tabPanel("Species select", "Overview of the distribution of studies for a selected species."),
  tabPanel("Diet", "Diet overview of the selected species.")),
  
  column(3,offset = 5, titlePanel("CarniDIET 1.0")),

  leafletOutput("mymap"),
  
  p(),
  
  headerPanel("Species selection"),
  
  sidebarPanel(
    
    uiOutput("select_var1"), 
    uiOutput("select_var2")
    
  ),
  
  mainPanel(
    tableOutput("table")
  )
)

# server.R ----
server <- function(input, output, session) {
  
  # current.ranges <- reactive({
  #   ranges %>% filter(Species == input$var2)
  # })
  # 
  tab <- reactive({ 
    
    tib %>% 
      filter(Family == input$var1) %>% 
      filter(CarniBinom == input$var2)
    
  })
  
  output$select_var1 <- renderUI({
    
    selectizeInput('var1', 'Select family', choices = c("select" = "", levels(tib$Family)))
    
  })
  
  output$select_var2 <- renderUI({
    
    
    choice_var2 <- reactive({
      tib %>% 
        filter(Family == input$var1) %>% 
        pull(CarniBinom) %>% 
        as.character()
      
    })
    
    selectizeInput('var2', 'Select species', choices = c("select" = "", choice_var2()))
    
  })
  
  output$table <- renderTable({ 
    
    tab()
    
  })
  
  # Create map to show all data...
  # and updated as you select species
  output$mymap <- renderLeaflet ({
    leaflet() %>%
      addTiles() %>%
      addCircles(data = tib,
                 lng = ~long,
                 lat = ~lat,
                 popup = ~CarniBinom,
                 color = "lightgrey") %>% 
      addCircles(data = tab(),
                 lng = ~long,
                 lat = ~lat,
                 color = "black",
                 weight = 5) 
    # addPolygons(data = range(),
    #             color = NULL,
    #             smoothFactor = 1,
    #             opacity = 0.5)
    
  })
  
}

# This apparently can update the map

# observe({
#   leafletProxy("map", data = tab()) %>%
#     clearShapes() %>%
#     addCircleMarkers(lat = ~lat,
#                      lng = ~long,
#                      popup = ~CarniBinom)
# })

# Run the app ----
shinyApp(ui = ui, server = server)
