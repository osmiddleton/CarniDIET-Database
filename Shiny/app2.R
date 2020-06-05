
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

# Libraries
library(shiny)
library(leaflet)
library(sf)
library(sp)
library(tidyverse)

# Load data
ecos <- st_read("./Data/GIS/ecoregions_biomes/biomes", "global_biomes")

#IUCN ranges
ranges <- st_read("./Data/GIS/TERRESTRIAL_MAMMALS", "TERRESTRIAL_MAMMALS")

# # CarniDIET
carnidiet <- read.csv("./Data/CarniDIET Database/CarniDIET_Full_Current.csv")
carnidiet <- carnidiet[c(2:75)]

# Tidy binomial names to match IUCN
carnidiet$CarniBinom <- as.character(carnidiet$CarniBinom)
carnidiet$CarniBinom <- gsub("_", " ", carnidiet$CarniBinom)
carnidiet$CarniBinom <- as.factor(carnidiet$CarniBinom)

temporal.diet <- carnidiet %>% 
  dplyr::group_by(Family,CarniBinom, Age, Sex,
                  PrimaryRef, Country, Region,
                  SiteName, AltitudeMaximum,
                  Season,Start.Year, End.Year,BIOME,
                  Longitude.Centroid.Mean.x, Latitude.Centroid.Mean.x,Scat.Stomach.Tissue) %>% 
  dplyr::summarise(x = length(CarniBinom))

temporal.diet <- temporal.diet[!is.na(temporal.diet$Latitude.Centroid.Mean.x),]
temporal.diet <- temporal.diet[!is.na(temporal.diet$Longitude.Centroid.Mean.x),]

temporal.diet$Latitude.Centroid.Mean.x <- as.numeric(as.character(temporal.diet$Latitude.Centroid.Mean.x))
temporal.diet$Longitude.Centroid.Mean.x <- as.numeric(as.character(temporal.diet$Longitude.Centroid.Mean.x))
str(temporal.diet)

temporal.diet %>% filter(Family == "Mephitidae") %>% pull(CarniBinom)
tib <- as_tibble(temporal.diet)
show(tib)
show(temporal.diet)

colnames(temporal.diet)[14:15] <- c("long","lat")
colnames(tib)[14:15] <- c("long","lat")

# The temporary CarniDIET table has been made, now let's try to visualise:
# Shiny includes making a UI side of the app followed by a server side

# App 1: just returns  the data frame for the species that was selected...
shinyApp(
  ui = bootstrapPage(
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
  ),
  
  server = function(input, output, session) {
    
    range <- reactive({
      ranges %>% filter(binomial == input$var2)
    })
    
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
      
      selectizeInput('var2', 'Select species', choices = c("select" = "", choice_var2())) # <- put the reactive element here
      
    })
    
    output$table <- renderTable({ 
      
      tab()
      
    })

    # Need to have a map in the first place which show ALL data
    # Also filter through IUCN range maps and add the species' geographic range to the map
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
                   weight = 5) %>% 
        addPolygons(data = range(),
                    color = NULL,
                    smoothFactor = 1,
                    opacity = 0.5)
        
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
)


