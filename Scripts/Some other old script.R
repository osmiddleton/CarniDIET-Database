#————————————————————————————————————————————————————————————————————————————————————————————####
#
#  ---- CarniDIET Version 1.0 ----
#   
#   A database of the diets of the world's terrestrial, carnivorous mammals
#   
#   Script purpose: 
#   
#   All code used to summarise data in CarniDIET and produce figures 
#   for the the manuscript:
#   Middleton, O.S, Svensson, H, Scharlemann, J.P.W,Faurby, S, Sandom, C.J.
#   CarniDIET 1.0: A database of terrestrial carnivorous#   mammal diets.
#   Global Ecology and Biogeography. In review.
#   DOI: 
#   
#   Sections:
#   
#   - 1. Database metadata values
#   
#   - 2. Carnivorous mammals taxonomic representation:
#        - Percentage of species in a Family with at least one dietary study completed
#        - Geographic range of missing species (which biomes have the most species unstudied?)
#        - Number of studies per species, highlighting those species which make up the majority of the total studies (50%)
#   
#   - 3. Spatio-temporal information of studies:
#        - Spatio-temporal coverage (latitude~Years studied) (i.e. PREDICTS)
#        - Global maps of study distributions
#   
#   - 4. Resolution of dietary composition for carnivorous mammals:
#        - Number of species
#        - Percentage recorded to different taxonomic resolutions
#        
#   - 5. Methodology used to quantify diets:
#        - Percentage of records made using different methods
#       
#   A note on taxonomy:
#   Taxonomy follows that of Phylacine v1.2 (Faurby et al., 2018) due to the
#   purpose of it being built to be a standardised database for macroecolgoical
#   research on mammals during the Late Quaternary.
#————————————————————————————————————————————————————————————————————————————————————————————####

#————————————————————————————————————————————————————————————————————————————————————————————####
# Load packages and data  ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# Packages ----
library(pacman)
pacman::p_load(tidyverse,plyr,ggrepel,sp,sf,viridis,rgdal,
               rgeos,maptools,reshape2,RColorBrewer,gridExtra,
               grid,doParallel,effects,raster,cluster,factoextra,
               vegan,FD,factoextra,pscl,boot,mapproj,ggplot2, dplyr,
               stringr, gridExtra,viridisLite,raster,rgdal,maptools,
               ape, ggtree,ggpubr, update = F)


# Data ----

# In the CarniDIET database

# CarniDIET 1.0:
carnidiet <- read.csv("./Version 1.0/CarniDIET 1.0.csv")

# List of potential mammal-consumers:
# Includes species selected as primary consumers of mammals from MammalDIET (Kissling et al., 2014)
cd.wos.hits <- read.csv("./Version 1.0/Supplementary data/Potential species list.csv", header = TRUE)

# External downloads
# Note: you will need to have downloaded these datasets externally to CarniDIET and have
# them stored in a directory immediately before thiscurrent directory.

# Phylacine (https://github.com/MegaPast2Future/PHYLACINE_1.2):
phylacine <- read.csv("../Phylacine/Trait_data.csv")

# IUCN species ranges (https://www.iucnredlist.org/resources/spatial-data-download):
ranges <- st_read("../GIS/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp") 

# Ecoregions (https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world):
ecoregions <- st_read("../GIS/ecoregions_biomes/ecoregions/wwf_terr_ecos.shp")

#————————————————————————————————————————————————————————————————————————————————————————————####
# 1. CarniDIET Metadata ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# 1. Mammal-consumer species ----

# How many ... in CarniDIET
nlevels(carnidiet$scientificNameCarni) # 103 species names
nlevels(carnidiet$commonNameCarni) # 103 common names
nlevels(carnidiet$orderCarni) # 5 orders
nlevels(carnidiet$familyCarni)  # 15 families

# How many levels of:
nlevels(carnidiet$lifeStageCarni) # 4 ages/life-stages
nlevels(carnidiet$sex) # Male and female and both



# 2. Diet resolution -----

# How many levels for food-type variables:
# the '-1' accounts for "" entries which could be treated as NA's...
nlevels(carnidiet$foodType) # 15 Food types
nlevels(carnidiet$orderPrey) # 60 prey orders
nlevels(carnidiet$familyPrey) # 121 prey families
nlevels(carnidiet$genusPrey) # 472 prey genera
nlevels(carnidiet$speciesPrey) # 712 prey species
nlevels(carnidiet$scientificNamePrey) # 824 prey binomials
nlevels(carnidiet$commonNameprey) # 1635 common names

# How many domestic animals
dom <- carnidiet %>% filter(domesticOrAgricultural == 1) %>% droplevels()
levels(dom$scientificNamePrey)

# Categories of taxonomic precision
nlevels(carnidiet$taxonRankPrey) # 10 categories of taxonomic precision

# Range of quantitative importance of prey species or food type
range(carnidiet$percentage) # 0.005% to 100%

# Range of sample sizes used for quantifying diet compositions
range(carnidiet[!is.na(carnidiet$sampleSizeScatStomachTissue),]$sampleSizeScatStomachTissue) # 2 to 11,478
range(carnidiet[!is.na(carnidiet$sampleSizeKillsPreyItems),]$sampleSizeKillsPreyItems) # 3 to 23,824



# 3. Temporal ----

# Range in study years
range(carnidiet[!is.na(carnidiet$startYear),]$startYear) # 1933 to 2017
range(carnidiet[!is.na(carnidiet$endYear),]$endYear) # 1946-2017

# Months
nlevels(carnidiet[!is.na(carnidiet$startMonth),]$startMonth)
nlevels(carnidiet[!is.na(carnidiet$endMonth),]$endMonth)
# All months represented in start and month of study

# Day of the month the study started
nlevels(factor(carnidiet[!is.na(carnidiet$startDayOfYear),]$startDayOfYear))
nlevels(factor(carnidiet[!is.na(carnidiet$endDayOfYear),]$endDayOfYear))

# Seasons - these are pretty crude
levels(carnidiet$season) # 45 seasons


# 4. Methods ----

# Quantification methods and sampling protocols
levels(carnidiet$methodQuantification) # 9 Quantification methods
levels(carnidiet$samplingProtocol) # 20 diet sampling protocols



# 5. Spatial information -----

# Altitude
range(carnidiet[!is.na(carnidiet$maximumElevationInMeters),]$maximumElevationInMeters) # 0 - 8156m
range(carnidiet[!is.na(carnidiet$minimumElevationInMeters),]$minimumElevationInMeters) # 0 - 4543m


# Latitude and longitude
range(carnidiet[!is.na(carnidiet$decimalLatitude),]$decimalLatitude) # -55S to 81N
range(carnidiet[!is.na(carnidiet$decimalLongitude),]$decimalLongitude) # -169W to 174E
levels(carnidiet$georeferenceSources) # 4 Types of geo-referencing


# Descriptions of study area
levels(carnidiet$geographicRegion) # 139 broad geographic regions
levels(carnidiet$country) # 98 countries
levels(carnidiet$stateProvince) # 156 higher level administrative units 
levels(carnidiet$county) # 81 lower level administrative units
levels(carnidiet$municipality) # 118 lower level administrative units
levels(carnidiet$verbatimLocality) # 344 general descriptions or non-official name places
levels(carnidiet$protectedAreaHigher) # 256 areas of some level of protection
levels(carnidiet$protectedAreaLower) # 3 lower level descriptions of site in protected area
levels(carnidiet$islandGroup) # 3 island groups
levels(carnidiet$island) # 28 islands


# Range of study area size
carnidiet$studyAreaSize <- as.numeric(as.character(carnidiet$studyAreaSize))
range(carnidiet[!is.na(carnidiet$studyAreaSize),]$studyAreaSize) # 0.03 to 100,000km2
summary(carnidiet$studyAreaSize)

# ---- 6. Bibliographic information ----
levels(carnidiet$sourcePrimaryReference) # 719
range(carnidiet$sourceYear) # 1952 to 2019
levels(carnidiet$sourceJournal) # 196 Journals
levels(carnidiet$sourceCollectionReference) # 690 Collection Refs

# 7. Number of studies ----

# A study is a subset of a source and can contain multiple breakdowns of a species' diet composition between years, or seasons and with multiple
# quantification methods.
colnames(carnidiet)

# Find all 
carni.studies <- carnidiet %>% 
  dplyr::group_by(familyCarni, scientificNameCarni, lifeStageCarni, sexCarni, # Unique species and demography
                  
                  # Unique methods
                  samplingProtocol, methodQuantification,
                  
                  # Unique spatial location
                  geographicRegion,country,stateProvince,county,municipality,verbatimLocality,protectedAreaHigher,
                  protectedAreaLower,islandGroup,island,maximumElevationInMeters,
                  decimalLatitude, decimalLongitude,
                  
                  # Unique temporal information
                  startYear, endYear, startMonth, endMonth, season, 
                  
                  # Unique source
                  sourcePrimaryReference) %>%  
  
  dplyr::summarise(x = length(scientificNameCarni))


# Need to identify how many of these studies have multiple seasons, multiple quantification methods, and multiple years

carni.studies <- carni.studies %>% 
  dplyr::group_by(familyCarni, scientificNameCarni, lifeStageCarni, sexCarni, # Unique species and demography
                           
                           # Unique spatial location
                           geographicRegion,country,stateProvince,county,municipality,verbatimLocality,protectedAreaHigher,
                           protectedAreaLower,islandGroup,island,maximumElevationInMeters,
                           decimalLatitude, decimalLongitude,
                           
                           # Unique source
                  sourcePrimaryReference) %>% 
  dplyr::summarise(seasons = length(unique(season)),
                   years = length(unique(startYear)),
                   methods = length(unique(methodQuantification)),
                   sampling = length(unique(samplingProtocol)))


# How many are singles:
nrow(carni.studies %>% filter(seasons == 1 & years == 1)) # 196
994/1310



# How many are seasonal comparisons:
nrow(carni.studies %>% filter(seasons > 1 & years == 1)) # 196
196/1310

# How many are time-series?
nrow(carni.studies %>% filter(seasons == 1 & years > 1)) # 51
51/1310

# How many are seasonal-comparison time-series?
nrow(carni.studies %>% filter(seasons > 1 & years > 1)) # 71
71/1310


#————————————————————————————————————————————————————————————————————————————————————————————####
# Potential mammal-consumer species covered by diet studies  ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# How many potential ...?
# Species
length(unique(cd.wos.hits$Bin.)) # 210 species
nlevels(cd.wos.hits$Order) # 9 orders
nlevels(cd.wos.hits$Family) # 23 families

# How many species in each family?
table(cd.wos.hits$Family)

# ---- Family representation in diet sources -----
# Number of POTENTIAL species in CarniDIET (by family)
family.species <- cd.wos.hits %>% dplyr::group_by(Family) %>% dplyr::summarize(No.Species = length(Family))
family.species
colnames(family.species)[1] <- "familyCarni"

# Number of ACTUAL species in CarniDIET (by family)
family.CD <- carnidiet[c(3:4)]
family.CD <- unique(family.CD)
family.CD <- family.CD %>% dplyr::group_by(familyCarni) %>% dplyr::summarise(No.Species = length(familyCarni))
family.CD

# Number of WoS hits per Family
family.wos.hits <- cd.wos.hits %>% dplyr::group_by(Family) %>% dplyr::summarize(No.Hits = sum(WoS.Hits))
family.wos.hits
colnames(family.wos.hits)[1] <- "familyCarni"

# Number of sources included in CarniDIET (by Family)
family.papers <- carnidiet[c(3,55)]
family.papers <- unique(family.papers)
family.papers <- family.papers %>% dplyr::group_by(familyCarni) %>% dplyr::summarise(No.papers = length(familyCarni))
family.papers

# Median and max of number of sources per species (per Family)
family.species.papers <- carnidiet[c(3:4,55)]
family.species.papers <- unique(family.species.papers)
family.species.papers <- family.species.papers %>% dplyr::group_by(familyCarni,scientificNameCarni) %>% dplyr::summarise(studiesperspecies = length(scientificNameCarni))
(family.species.papers <- family.species.papers %>% dplyr::group_by(familyCarni) %>% dplyr::summarise(median.papers = median(studiesperspecies),
                                                                                                max.papers = max(studiesperspecies)))

# Building Family-level table
family.info <- merge(family.species,family.CD, by = "familyCarni", all = TRUE)
family.info$Perc.Families <- (family.info$No.Species.y/family.info$No.Species.x)*100
family.info <- merge(family.info, family.wos.hits, by = "familyCarni", all = TRUE)
family.info <- merge(family.info, family.papers, by = "familyCarni", all = TRUE)
family.info$Hits_Papers_Ratio <- (family.info$No.papers/family.info$No.Hits)*100
family.info <- merge(family.info, family.species.papers, by = "familyCarni", all = TRUE)
family.info[is.na(family.info)] <- 0

# How many sources were used compared to actual hits?
sum(family.info$No.Hits)
sum(family.info$No.papers)
median(family.info$Hits_Papers_Ratio)

# Merge in Order information
colnames(phylacine)
order.family <- unique(phylacine[c(2,3)])
colnames(order.family) <- c("orderCarni", "familyCarni")

family.info <- merge(family.info,order.family, by = "familyCarni", all.x = TRUE)
family.info <- family.info[c(10,1:9)]

family.info <- family.info %>% arrange(orderCarni, familyCarni)


# ---- Genus representation in diet sources -----

# Number of potential species in each genus
genus.species <- cd.wos.hits %>% dplyr::group_by(Family, Genus_V_1.2) %>%
              dplyr::summarize(No.Species = length(Genus_V_1.2))

# Number of actual species included in carnidiet
diet.2 <- carnidiet %>% separate(scientificNameCarni, into = paste(c("Genus_V_1.2","Species"), sep = "_"))

genus.CD <- diet.2[c(3:4,5)]
genus.CD <- unique(genus.CD)
genus.CD <- genus.CD %>% dplyr::group_by(familyCarni,Genus_V_1.2) %>% dplyr::summarise(No.Species = length(Genus_V_1.2))

# Number of hits per genus
genus.wos.hits <- cd.wos.hits %>% dplyr::group_by(Family,Genus_V_1.2) %>% dplyr::summarize(No.Hits = sum(WoS.Hits))

# Number of actual sources per genus
genus.papers <- diet.2[c(3:4,56)]
genus.papers <- unique(genus.papers)
genus.papers <- genus.papers %>% dplyr::group_by(familyCarni,Genus_V_1.2) %>% dplyr::summarise(No.papers = length(Genus_V_1.2))

# Median and max of number of sources per species per Genus
genus.species.papers <- diet.2[c(3:5,56)]
genus.species.papers <- unique(genus.species.papers)
genus.species.papers <- genus.species.papers %>% dplyr::group_by(familyCarni,Genus_V_1.2,Species) %>% dplyr::summarise(studiesperspecies = length(Species))
genus.species.papers <- genus.species.papers %>% dplyr::group_by(familyCarni,Genus_V_1.2) %>% dplyr::summarise(median.papers = median(studiesperspecies),
                                                                                                          max.papers = max(studiesperspecies))
# Building Genus-level table
genus.info <- merge(genus.species,genus.CD, by = "Genus_V_1.2", all = TRUE)
genus.info$Perc.Genus <- (genus.info$No.Species.y/genus.info$No.Species.x)*100
genus.info <- merge(genus.info, genus.wos.hits, by = "Genus_V_1.2", all = TRUE)
genus.info <- merge(genus.info, genus.papers, by = "Genus_V_1.2", all = TRUE)
genus.info$Hits_Papers_Ratio <- (genus.info$No.papers/genus.info$No.Hits)*100
genus.info <- merge(genus.info, genus.species.papers, by = "Genus_V_1.2", all = TRUE)
genus.info <- genus.info[c(1:3,5:6,8,10,11,13:14)]
genus.info[is.na(genus.info)] <- 0
colnames(genus.info)[2] <- "familyCarni"
genus.info <- merge(genus.info, order.family, by = "familyCarni", all.x = TRUE)
genus.info <- genus.info[c(11,1:10)]

genus.info <- genus.info %>% arrange(orderCarni, familyCarni, Genus_V_1.2)


# ---- Tables: Family & genus level representation in dietary studies ----
#write.csv(family.info, "../Tables/Family.Information.Summary_UPDATE.csv")
#write.csv(genus.info, "../Tables/Genus.Information.Summary_UPDATE.csv")


# ---- Supplementary figure: Scatterplot family-level diet representation ----

# Get diagonal representation lines for plot
# 75%
y = log10(7.5) # Species included 0.8750613
x = log10(10) # Total species 1
c75 = y - (1*x)

# 50%
y = log10(5) # Species included 0.8750613
x = log10(10) # Total species 1
c50 = y - (1*x)

# 25%
y = log10(2.5) # Species included 0.8750613
x = log10(10) # Total species 1
c25 = y - (1*x)

# 10%
y = log10(1) # Species included 0.8750613
x = log10(10) # Total species 1
c10 = y - (1*x)

# 5%
y = log10(0.5) # Species included 0.8750613
x = log10(10) # Total species 1
c5 = y - (1*x)

# 1%
y = log10(0.1) # Species included 0.8750613
x = log10(10) # Total species 1
c1 = y - (1*x)

(plot <- ggplot(data = family.info[!family.info$Perc.Families == 0,], aes(x = log10(No.Species.x), y = log10(No.Species.y), # Removed species with no studies
                                                                fill = familyCarni, size = log10(No.papers)), label = familyCarni) +
  geom_abline(intercept = 0, slope = 1, colour = "grey") +
  geom_abline(intercept = c75, slope = 1, colour = "grey", linetype = 2) +
  geom_abline(intercept = c50, slope = 1, colour = "grey", linetype = 3) +
  geom_abline(intercept = c25, slope = 1, colour = "grey", linetype = 4) +
  geom_abline(intercept = c10, slope = 1, colour = "grey", linetype = 5) +
  annotate("text", x = 1.7, y=1.75, label = "100%", size = 3) +
  annotate("text", x = 1.7, y=1.6, label = "75%", size = 3) +
  annotate("text", x = 1.7, y=1.45, label = "50%", size = 3) +
  annotate("text", x = 1.7, y=1.15, label = "25%", size = 3) +
  annotate("text", x = 1.7, y=0.8, label = "10%", size = 3) +
  geom_point(aes(fill = familyCarni), alpha = 0.5, pch=21, colour = "black") +
  scale_fill_manual(values = colours, guide = "none") +
  geom_text_repel(aes(label = familyCarni),
                  size = 2, colour = "black", force = 4) +
  xlim(0,1.8) + ylim(0,1.75)+
  theme_bw() + xlab("log10(number of mammal-consumers in Family)") +
    ylab("log10(number of mammal-consumers with diet study)") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9)))

# Save figure
ggsave(filename = paste0("../Final figures/Figure S1/PercentFamiliesStudiedUpdated.tiff"), plot = plot, device = NULL, path = NULL,
       scale = 1, width = 10.5, height = 10.5, units = c("cm"),
       dpi = 300)


# Supplementary figure: Barplot for top-20 studied species in database ----

# There are 719 sources used. How many feature specific species?

# Number of sources feature each species
carnivore.summary <- carnidiet %>%
                     dplyr::group_by(familyCarni, scientificNameCarni) %>%
                     dplyr::summarise(papers = length(unique(sourcePrimaryReference)))
 
# There are 1312 studies. How many for each species?

# Number of studies for each species
carni.studies <- carnidiet %>% 
  dplyr::group_by(familyCarni, scientificNameCarni, lifeStageCarni, sexCarni, # Unique species and demography
                  
                  # Unique spatial location
                  geographicRegion,country,stateProvince,county,municipality,verbatimLocality,protectedAreaHigher,
                  protectedAreaLower,islandGroup,island,maximumElevationInMeters,
                  decimalLatitude, decimalLongitude,
                  
                  # Unique source
                  sourcePrimaryReference) %>% 
  dplyr::summarise(x = length(scientificNameCarni))

carnivore.summary.3 <- carni.studies %>%
  dplyr::group_by(familyCarni, scientificNameCarni) %>%
  dplyr::summarise(studies = length(scientificNameCarni)) # Number of studies

carnivore.summary.master <- merge(carnivore.summary, carnivore.summary.3, by = "scientificNameCarni")

# How many studies are there
sum(carnivore.summary.master$studies) # 1312 - This is all individual studies
1310/719 # On average, there are 1.8 studies per source

# What makes 50% studies?
# 656 = 50% studies
# 984 = 75% studies

# Tidy table
colnames(carnivore.summary.master)[1] <- "CarniBinom"

# Get top-20 species
carnivore.summary.master$CarniBinom <- reorder(carnivore.summary.master$CarniBinom, carnivore.summary.master$studies)
x <- carnivore.summary.master[order(carnivore.summary.master$studies, decreasing = TRUE),]
y <- x[1:20,]
y$CarniBinom <- factor(y$CarniBinom)
levels(y$CarniBinom)

# How many species need to be included to get 50% of the studies?
(sum(x[1:5,]$studies)/1310)*100 #  9 species includes just over 50% studies
(sum(x[1:12,]$studies)/1310)*100 #  9 species includes just over 50% studies
(sum(x[1:20,]$studies)/1310)*100 # 20 species includes 70% studies

# Rename species factor levels for figure
y$CarniBinom <- revalue(y$CarniBinom, c("Vulpes_vulpes" = "V. vulpes", "Lynx_lynx" = "L. lynx",
                                        "Canis_lupus" = "C. lupus", "Felis_silvestris" = "F. silvestris",
                                        "Panthera_pardus" = "P. pardus","Neovison_vison" = "N. vison",
                                        "Puma_concolor" = "P. concolor","Martes_foina" = "M. foina",
                                        "Canis_latrans" = "C. latrans","Acinonyx_jubatus" = "A. jubatus",
                                        "Ursus_arctos" = "U. arctos","Canis_aureus" = "C. aureus",
                                        "Martes_martes" = "M. martes","Mustela_erminea" = "M. erminea",
                                        "Panthera_tigris" = "P. tigris","Panthera_uncia" = "P. uncia",
                                        "Martes_americana" = "M. americana","Lynx_rufus" = "L. rufus",
                                        "Crocuta_crocuta" = "C. crocuta","Panthera_leo" = "P. leo"))


# Main plot
(p10 <- ggplot(data = y[1:20,], aes(x = reorder(CarniBinom, studies), y = studies)) + # 20 species refers to 66% all papers
  geom_bar(stat = "identity", fill = "darkgrey") +
  #annotate("rect", xmin = 19.5, xmax = 20.5, ymin = -Inf, ymax =Inf, fill = "red") +
  #annotate("rect", xmin = 15.5, xmax = 19.5, ymin = -Inf, ymax =Inf, fill = "orange") +
  #annotate("rect", xmin = 7.5, xmax = 15.5, ymin = -Inf, ymax =Inf, fill = "yellow") +
  coord_flip() +
  xlab("") + ylab("Number of studies") +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Sources") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 17.5, colour = "black", linetype = 2, alpha = 0.8) + # 15% studies
  geom_vline(xintercept = 8.5, colour = "black", linetype = 2, alpha = 0.8) + # 50% studies
  annotate("text", x = 17, y = 190, label = "25%", size = 3) +
  annotate("text", x = 8, y = 190, label = "50%", size = 3) +
  
    # Add text for each species
  annotate("text", x = 20, y = 175, label = "(92)", size = 2, colour = "black") +
    annotate("text", x = 19, y = 92, label = "(61)", size = 2, colour = "black") +
    annotate("text", x = 18, y = 70, label = "(55)", size = 2, colour = "black") +
    annotate("text", x = 17, y = 65, label = "(41)", size = 2, colour = "black") +
    annotate("text", x = 16, y = 57, label = "(35)", size = 2, colour = "black") +
    annotate("text", x = 15, y = 45, label = "(31)", size = 2, colour = "black") +
    annotate("text", x = 14, y = 43, label = "(19)", size = 2, colour = "black") +
    annotate("text", x = 13, y = 43, label = "(21)", size = 2, colour = "black") +
    annotate("text", x = 12, y = 40, label = "(20)", size = 2, colour = "black") +
    annotate("text", x = 11, y = 36, label = "(21)", size = 2, colour = "black") +
    annotate("text", x = 10, y = 36, label = "(22)", size = 2, colour = "black") +
    annotate("text", x = 9, y = 36, label = "(25)", size = 2, colour = "black") +
    annotate("text", x = 8, y = 34, label = "(25)", size = 2, colour = "black") +
    annotate("text", x = 7, y = 32, label = "(15)", size = 2, colour = "black") +
    annotate("text", x = 6, y = 32, label = "(20)", size = 2, colour = "black") +
    annotate("text", x = 5, y = 32, label = "(14)", size = 2, colour = "black") +
    annotate("text", x = 4, y = 31, label = "(19)", size = 2, colour = "black") +
    annotate("text", x = 3, y = 29, label = "(16)", size = 2, colour = "black") +
    annotate("text", x = 2, y = 27, label = "(16)", size = 2, colour = "black") +
    annotate("text", x = 1, y = 26, label = "(11)", size = 2, colour = "black") +
  
    # Add changes to axes  
  ylim(0,200) + theme(axis.text.y = element_text(face = "italic")))

# Save plot
ggsave(filename = paste0("../Final figures/Figure S2/Top-20 studied species.tiff"), plot = p10, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c("cm"),
       dpi = 200)



# Figure: Phylogenetic tree and number of studies per species ----

# Download phylogeny from PHYLACINE
# Save  phylogeny in  repository behind the one with this project
forest <- read.nexus("../Phylacine/Phylogenies/Complete_phylogeny.nex")
names(forest) <- NULL
set.seed(42)
forest <- forest[sample(1:1000, 30)]

# Species list in CarniDIET
(sp <- levels(carnivore.summary.master$Binomial.1.2))

# Trim tree down to CarniDIET species list
pruned.forest <- lapply(forest, 
                        function(x) {
                          drop.tip(x, x$tip.label[!x$tip.label %in% sp])
                        }
)
class(pruned.forest) <- "multiPhylo"

# Pick a single tree for overplotting and labelling
tree <- pruned.forest[1]
plot(tree)
summary(tree)

# Turn tree to dataframe
tree <- fortify(tree)

# Change columns name in carnivore summary to match Phylacine
colnames(carnivore.summary.master)[1] <- "Binomial.1.2"

# Add number of studies to the tree
tree <- left_join(tree, carnivore.summary.master, by = c("label" = "Binomial.1.2"))

# Reverse age scale
tree$x <- tree$x - max(tree$x)

# Convert multiPhylo to data.frame
#pruned.forest <- fortify(pruned.forest)

# Reverse ages to show time BP
#pruned.forest$x <- pruned.forest$x - max(pruned.forest$x)

# Create new name binomial names
foo <- data.frame(do.call('rbind', strsplit(as.character(tree$label),'_',fixed=TRUE)))
tree$Genus <- substring(foo$X1,1,1)
tree$Species <- foo$X2
tree$label2 <- paste(tree$Genus, tree$Species, sep = ". ")

# What happens if you change the branch length?
tree$branch.length <- tree$branch.length/2

# Update colour scheme
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

# Plot phylogeny
(p <- ggtree(tree, layout='circular', ladderize = FALSE) +
    geom_tiplab(size=1.5, aes(label = label2, parse = T, size= 3.5, angle=angle), offset = 10) +
    geom_tippoint(aes(x, y, fill = familyCarni.x, size = studies), pch = 21, colour = "black", alpha = 0.75) +
    scale_fill_manual(values = colours) +
    #xlim(-160, 30) +
    theme(panel.background = element_blank(),
          legend.position = "none"))

# Save phylogenetic tree
ggsave("../Final figures/Figure 2/Phylo tree_updated.pdf", p,
       width = 4, height = 4, units = "in")


# Barplot showing decline in studies and coloured by family

(barplot <- ggplot(data = carnivore.summary.master, aes(x = reorder(Binomial.1.2, -studies),
                                                        fill = familyCarni.x,
                                                     y = studies)) +
    geom_bar(stat = "identity", width = 1) +
    # scale_y_continuous(labels = c("0" = "1",
    #                                "1" = "10",
    #                                "2" = "100",
    #                                "3" = "1000")) +
    scale_fill_manual(values = colours, guide = "none") +
    labs(x = "Species", y = "Number of studies") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 5)))

ggsave("../Final figures/Figure 2/barplot updated.pdf", barplot,
       width = 1.75, height = 1.75, units = "in")



# Add figures together later in Inkscape


# Supplementary figure: Repeat with all mammal-consumers ----

# Reload phylogeny
forest <- read.nexus("../Phylacine/Phylogenies/Complete_phylogeny.nex")
names(forest) <- NULL
set.seed(42)
forest <- forest[sample(1:1000, 30)]

# Get all potential mammal-consumers
(sp <- levels(cd.wos.hits$Bin.))

# Trim tree down to the species list of MEMs created earlier
pruned.forest <- lapply(forest, 
                        function(x) {
                          drop.tip(x, x$tip.label[!x$tip.label %in% sp])
                        }
)
class(pruned.forest) <- "multiPhylo"

# Pick a single tree for overplotting and labelling
tree <- pruned.forest[1]
plot(tree)
summary(tree)

# Turn tree to dataframe
tree <- fortify(tree)

# Change columns name in carnivore summary
colnames(carnivore.summary.master)[1] <- "Binomial.1.2"

# Add number of studies to tree
tree <- left_join(tree, carnivore.summary.master, by = c("label" = "Binomial.1.2"), all.x = TRUE)

# Reverse age scale
tree$x <- tree$x - max(tree$x)

# Convert multiPhylo to data.frame
pruned.forest <- fortify(pruned.forest)

# Reverse ages to show time BP
pruned.forest$x <- pruned.forest$x - max(pruned.forest$x)

# Update binomial names for aesthetics
foo <- data.frame(do.call('rbind', strsplit(as.character(tree$label),'_',fixed=TRUE)))
tree$Genus <- substring(foo$X1,1,1)
tree$Species <- foo$X2
tree$label2 <- paste(tree$Genus, tree$Species, sep = ". ")

# What happens if you change the branch length?
tree$branch.length <- tree$branch.length/2

# Merge family info to new tree with all mammal-consumers
tree <- left_join(tree, cd.wos.hits, by = c("label" = "Bin."))
levels(tree$Family)

# Update colour scheme to match previous, although requires some hackery
levels(carnivore.summary.master$Family.x)
levels(tree$Family.x)

colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")
# Add random colours into the plot but keep colours for families in CarniDIET constant
colours.updated <- c("#CCCCCC", # Canidae IN
             "#CCCCCC", # Cricetidae OUT
             "#99CC00", # Dasyuridae IN
             "#999999", # Didelphidae IN
             "#CCCCCC", # Erinaceidae OUT
             "#339900", # Eupleridae IN
             "#CC9900", # Felidae IN
             "#FFCC00", # Gliridae IN
             "#FF6600", # Herpestidae IN
             "#993300", # Hyaenidae IN
             "#FFCCCC", # Lorisidae IN
             "#CCCCCC", # Megadermatidae OUT
             "#CC6699", # Mephitidae IN
             "#FF99FF", # Mustelidae IN
             "#663399", # Phyllostomidae IN
             "#CCCCCC", # Pitheciidae OUT
             "#9999CC", # Procyonidae IN
             "#CCCCCC", # Soricidae OUT
             "#CCCCCC", # Tenrecidae OUT
             "#CCCCCC", # Tupaiidae OUT
             "#3399FF", # Ursidae IN
             "#CCFFFF", # Viverridae IN
             "#003333") # Spare colour...
# Could have made a data frame of colour assignments to families and filtered, but meh!

# Plot phylogenetic tree
(p <- ggtree(tree, layout='circular', ladderize = FALSE) +
    geom_tiplab(size=1.5, aes(label = label2, parse = T, size= 3.5, angle=angle), offset = 10) +
    geom_tippoint(aes(x, y, fill = familyCarni.x, size = studies), pch = 21, colour = "black", alpha = 0.75) +
    scale_fill_manual(values = colours.updated) +
    #xlim(-160, 30) +
    theme(panel.background = element_blank(),
          legend.position = "none"))

# Save plot
ggsave("../Final figures/Figure S3/Phylo tree_ALL Species.pdf", p,
       width = 5, height = 5, units = "in")


#————————————————————————————————————————————————————————————————————————————————————————————####
# Spatial and temporal distribution of diet studies ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# Need to talk about this in the context of where species are distributed:
ranges$binomial <- gsub(" ", "_", ranges$binomial)
ranges$binomial <- as.factor(ranges$binomial)
ranges.1 <- ranges[ranges$binomial %in% c(levels(cd.wos.hits$Bin.)),]

# Convert ecoregion polygons into a shapefile per ecoregion
ecoregions.2 <- ecoregions %>% dplyr::group_by(ECO_NAME) %>% dplyr::summarize()

# Convert carnivore range polygons into a shapefile per carnivore
ranges.2 <- ranges %>% dplyr::group_by(binomial) %>% dplyr::summarize()

# Identify which ecoregions a species' range overlaps                 
new.df <- matrix(ncol = nlevels(cd.wos.hits$Bin.),
                 nrow = nrow(ecoregions.2))
colnames(new.df) <- levels(cd.wos.hits$Bin.)

# Run for-loop
for(i in 1:length(levels(cd.wos.hits$Bin.))) {
  
  if(length(unique(i == c(9,122,170,194))) > 1) {
    
    print(i)
    
  } else {
    
    df <- data.frame(eco = 1:nrow(ecoregions.2))
    test <- as.data.frame(st_intersects(ecoregions.2, ranges.2[ranges.2$binomial == levels(cd.wos.hits$Bin.)[i],]))
    test$col.id <- 1
    test <- unique(test)
    
    colnames(test)[1] <- "eco"
    df.1 <- merge(df, 
                  test, by = "eco", all.x = TRUE)
    
    new.df[,i] <- df.1$col.id
    
    print(i)
     
  }
  
}

# Tidy dataframe
new.df[is.na(new.df)] <- 0
ecoregions.2$SpecRich <- NA
ecoregions.2$SpecRich <- rowSums(new.df)

# Species richness of missing species
miss <- which(!colnames(new.df) %in% levels(carnidiet$scientificNameCarni))
missing.species <- new.df[, miss]

ecoregions.2$MissingSpecRich <- NA
ecoregions.2$MissingSpecRich <- rowSums(missing.species)

# Remove 'Rock & Ice' and 'Lake'
levels(ecoregions.2$ECO_NAME)
ecor.2.tidy <- ecoregions.2[!ecoregions.2$SpecRich > 30,]
str(ecor.2.tidy)
ecor.2.tidy$MissingSpecRich <- as.integer(ecor.2.tidy$MissingSpecRich)

# Plot ecoregion missing species richness
eco.missing.species <- ggplot() +
  geom_sf(data = ecor.2.tidy, aes(fill = MissingSpecRich), colour = NA) +
  geom_sf(data = world, fill = NA, colour = "black",lwd = 0.1) +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Species richness", limits = c(0,13), breaks = c(0,2,4,6,8,10,12)) +
  coord_sf(crs = "+proj=moll") +
  theme_classic() +
  theme(panel.grid = element_blank())

ggsave("../Final figures/Figure SX/Missing species richness.tiff", p1)

# Plot all ecoregion species richnes
eco.specie.richness <- ggplot() +
  geom_sf(data = ecor.2.tidy, aes(fill = SpecRich), colour = NA) +
  geom_sf(data = world, fill = NA, colour = "black",lwd = 0.1) +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Species richness", limits = c(0,30), breaks = c(0,5,10,15,20,25,30)) +
  coord_sf(crs = "+proj=moll") +
  theme_classic() +
  theme(panel.grid = element_blank())

ggsave("../Final figures/Figure SX/All species richness.tiff", p2)


# Get unique studies ----

# Find the number of studies but keep year information
temporal.diet <- carnidiet %>% 
  dplyr::group_by(familyCarni, scientificNameCarni, lifeStageCarni, sexCarni, # Unique species and demography
                  
                  # Unique spatial location
                  geographicRegion,country,stateProvince,county,municipality,verbatimLocality,protectedAreaHigher,
                  protectedAreaLower,islandGroup,island,maximumElevationInMeters,
                  decimalLatitude, decimalLongitude,
                  
                  # Unique temporal information
                  startYear, endYear, 
                  
                  # Unique source
                  sourcePrimaryReference) %>%  
  
  dplyr::summarise(x = length(scientificNameCarni))

nrow(temporal.diet) # 1574 studies

# Get midpoints for each study:
temporal.diet$midYear <- (temporal.diet$startYear + temporal.diet$endYear)/2

# Figure: Spatio-temporal distribution ----

# Find median and IQR of study mid-year
(midYear <- temporal.diet$midYear[!is.na(temporal.diet$midYear)])

(median.year <- quantile(midYear)[3])
(LQ.year <- quantile(midYear)[2])
(UQ.year <- quantile(midYear)[4])

# Find median and IQR of study latitudes
latitudes <- temporal.diet$decimalLatitude[!is.na(temporal.diet$decimalLatitude)]

(median.lat <- quantile(abs(latitudes))[3])
(LQ.lat <- quantile(abs(latitudes))[2])
(UQ.lat <- quantile(abs(latitudes))[4])

# Get colour scheme consistent again:
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

# Plot figure showing temporal distribution and |latitude| of studies
p1 <- ggplot() +
  # annotate("rect",xmin =-Inf, xmax = Inf, ymin = LQ.year, ymax =UQ.year, alpha = .2) +
  # annotate("rect",xmin = LQ.lat, xmax = UQ.lat, ymin = -Inf, ymax =Inf, alpha = .2) +
  # annotate("rect",xmin = LQ.lat, xmax = UQ.lat,ymin = LQ.year, ymax = UQ.year, alpha = .2) +
  ylab("Year") + xlab("|Latitude|") +
  coord_flip() +
  scale_y_continuous(breaks = c(1940,1960,1980,2000,2020)) + 
  geom_errorbar(data = temporal.diet, 
                    aes(ymin = startYear, ymax = endYear, 
                        x = abs(decimalLatitude), colour = familyCarni), alpha = .5, width = 0) +
  geom_point(data = temporal.diet, 
               aes(y = midYear, x = abs(decimalLatitude),fill = familyCarni),
               pch = 21, alpha = .5, colour = "black", size = 1.5) +
  # geom_vline(xintercept = median.lat) + # Median latitude
  # geom_vline(xintercept = 23.5, linetype = 2) + # Tropic lines
  # annotate("text", x = 22, y = 1936.5, label = "Tropic lines", colour = "black", size = 3) +
  # geom_vline(xintercept = 0, linetype = 2) + # Equator
  # annotate("text", x = 2.5, y = 1935, label = "Equator", colour = "black", size = 3) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  #geom_hline(yintercept = median.year, colour = "black") + # Median year
  theme_bw() + theme(legend.position = "top",
                     legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank())

# Save figure
ggsave("../Final figures/Figure 1/Spatio-temporal information_UPDATED Again.tiff", p1,
        width = 10, height = 8, dpi = 100)



# Figure: Map showing study distribution ----

# Get biomes and dissolve into a world map
biomes <- st_read("../GIS/ecoregions_biomes/biomes/global_biomes.shp")
crs(biomes)
class(biomes)
world <- st_union(biomes)
class(world)

# Practice plot with different projection system
ggplot() + 
  geom_sf(data = world) +
  coord_sf(crs = "+proj=moll") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank())

# Convert study information into a spatial points dataframe
colnames(temporal.diet)

# Remove NAs in the dataframe
temporal.diet.no.nas <- temporal.diet[!is.na(temporal.diet$decimalLatitude),]
temporal.diet.no.nas <- temporal.diet[!is.na(temporal.diet$decimalLongitude),]

# change to sf object and will change to Mollweide projection...difficult
coords <- as.data.frame(temporal.diet.no.nas[c(17,16)])

# Get Spatial points
temporal.diet.sp <- SpatialPointsDataFrame(coords = coords,
                                  data = temporal.diet.no.nas,
                                  proj4string = CRS(crs(biomes)))

plot(temporal.diet.sp)

print(temporal.diet.sp)

temporal.diet.sp <- st_as_sf(temporal.diet.sp)

# Need to get colours right again...
levels(factor(temporal.diet$familyCarni))
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

levels(factor(temporal.diet.sp$familyCarni))
colours.updated.again <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
                           "#FFCCCC","#CC6699","#FF99FF","#9999CC","#3399FF","#CCFFFF","#003333")

# Plot map
p2 <- ggplot() +
  geom_sf(data = ecor.2.tidy, aes(fill = SpecRich), colour = NA, alpha =0.7) +
  scale_fill_gradient(low = "white", high = "black", guide = "Ecoregion SR") +
  geom_sf(data = world, fill = NA, colour = "black", lwd = 0.1, alpha = 0.4) +
  new_scale_fill() +
  geom_sf(data = temporal.diet.sp, aes(fill = familyCarni), 
          alpha = 0.5, pch = 21, colour = "black", size = 1.5) +
  coord_sf(crs = "+proj=moll") +
  labs(x = "", y = "") +
  scale_fill_manual(values = colours.updated.again) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank())

ggsave("../Final figures/Figure 1/Map but fpor the additional legend.pdf", p2,
       width = 8, height = 8, units = "in")

p3 <- ggarrange(p1,
                p2,
                ncol=1, nrow=2, labels = "auto",
                common.legend = TRUE, legend="right")
   
ggsave("../Final figures/Figure 1/Combined spatio-temporal_mollweide UPDATED AGAIN.pdf", p3, 
       width = 8, height = 8, units = "in")


# ---- Data summary: More information on temporal information ----

# What is the average length of time a study is reported for
x <- temporal.diet[!is.na(temporal.diet$endYear),]
y <- x[!is.na(x$startYear),]

median(y$endYear - y$startYear)
quantile(y$endYear - y$startYear)


# How many studies are there per ecoregion and is this associated with mammal-consumer species richness?
diet.sf <- st_as_sf(temporal.diet.sp)

# Get the intersection of points with ecoregions
ecor.2.tidy

test <- head(diet.sf)
test <- test[1,]

# get intersection between ecoregions and study points:
t <- as.data.frame(st_intersects(ecor.2.tidy, diet.sf))
t$count <- 1

# Summarise
eco.r.studies <- t %>% dplyr::group_by(row.id) %>% dplyr::summarise(study.count = sum(count))
colnames(eco.r.studies)[1] <- "ECO.ROW"

# Merge with ecoregions file
ecor.2.tidy$ECO.ROW <- 1:nrow(ecor.2.tidy)
ecoregions.3 <- merge(ecor.2.tidy, eco.r.studies, by = "ECO.ROW", all.x = TRUE)

ecoregions.3[is.na(ecoregions.3$study.count),]$study.count <- 0

cor.test(ecoregions.3$SpecRich, ecoregions.3$study.count)

(cor <- ggplot(data = ecoregions.3[ecoregions.3$SpecRich < 30,], aes(x = SpecRich, y = study.count)) +
  geom_jitter(pch = 21, colour = "black", fill = "lightgrey", size = 1.75, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") + 
  labs(x = "Mammal-consumer species richness per ecoregion", y = "Number of studies per ecoregion") +
  theme_bw() +
  theme(panel.grid = element_blank()))

ggsave("../Final figures/Figure SX/Ecoregion species richness plot.tiff", cor,
       width = 4, height = 4, units = "in")



# Plot map of number of studies in ecoregions
eco.studies <- ggplot() +
  geom_sf(data = ecoregions.3[ecoregions.3$study.count > 0,], aes(fill = study.count), colour = NA) +
  geom_sf(data = world, fill = NA, colour = "black",lwd = 0.1) +
  scale_fill_gradient(low = "darkgrey", high = "black", name = "No. of studies") +
  coord_sf(crs = "+proj=moll") +
  theme_classic() +
  theme(panel.grid = element_blank())

ggsave("../Final figures/Figure SX/test again.tiff", eco.studies)

# Plot all in a single plot:
all.ecoregion.info <- cowplot::plot_grid(eco.specie.richness, eco.missing.species,eco.studies,
                                         ncol = 1, nrow = 3, labels = "auto")

ggsave("../Final figures/Figure SX/Ecoregion information.tiff", all.ecoregion.info,
       width = 6, height = 6, units = "in")

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- Diet resolution ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Summarising data ----
# For each species get:
       
#       - Percentage of records at certain taxonomic resolutions
#       - Number of records using different sources
#       - Number of records using different quantification methods

# Number of records per carnivore and for each of the different methods
diet.records_perSource <- carnidiet %>%
                          dplyr::group_by(scientificNameCarni, samplingProtocol) %>%
                          dplyr::summarise(records = length(samplingProtocol))

# Get percentage of records for each sample source
perc <- diet.records_perSource %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarise(num = sum(records))
(perc$perc <- (perc$num/sum(perc$num))*100)

diet.records_perSource <- spread(diet.records_perSource, samplingProtocol, records)

# Methods
diet.records_perMethod <- carnidiet %>%
                          dplyr::group_by(scientificNameCarni, methodQuantification) %>%
                          dplyr::summarise(records = length(methodQuantification))

# Get percentage of records for each sample source
perc <- diet.records_perMethod %>% dplyr::group_by(methodQuantification) %>% dplyr::summarise(num = sum(records))
perc
perc$perc <- (perc$num/sum(perc$num))*100
perc

diet.records_perMethod <- spread(diet.records_perMethod, methodQuantification, records)

# Number of unique prey species identified in species diets
carnidiet$scientificNamePrey <- as.factor(carnidiet$scientificNamePrey)
nlevels(carnidiet$scientificNamePrey)
levels(carnidiet$scientificNamePrey)

colnames(carnidiet)
prey <- carnidiet[c(4,13)]
prey <- unique(prey)
prey <- prey[!is.na(prey$scientificNamePrey),]
prey$value <- 1

prey.consumed <- prey %>%
                 dplyr::group_by(scientificNameCarni) %>%
                 dplyr::summarise(prey.SR = length(scientificNameCarni))

# Percentage of records from prey orders
colnames(carnidiet)
prey <- carnidiet %>% filter(taxonRankPrey == "Species")

prey <- prey %>% dplyr::group_by(orderPrey) %>% dplyr::summarize(count = length(orderPrey)) 

prey$total <- sum(prey$count)

prey$taxa.perc <- (prey$count/prey$total)*100


# Percentage of records which are at different taxonomic resolutions
colnames(carnidiet)
prey <- carnidiet[c(4,8:13,16)]
prey$value <- 1

prey <- prey %>%
        group_by(scientificNameCarni, taxonRankPrey) %>%
        dplyr::summarise(value = length(taxonRankPrey))

# Get percentage of records for each taxonomic precision
perc <- prey %>% dplyr::group_by(taxonRankPrey) %>% dplyr::summarise(num = sum(value))
perc
perc$perc <- (perc$num/sum(perc$num))*100


prey <- as.data.frame(prey[c(1:3)])
prey

prey <- spread(prey, taxonRankPrey, value)
colnames(prey)
prey <- prey[c(1,9,4,3,6,2,8,5)]
prey[is.na(prey)] <- 0
total.obvs <- as.data.frame(rowSums(prey[c(2:8)]))

prey.percentage <- as.matrix(prey)
row.names(prey.percentage) <- prey.percentage[,1]
prey.percentage <- prey.percentage[,-c(1)]
prey.percentage[is.na(prey.percentage)] <- 0

prey.percentage <- as.data.frame(prey.percentage)

prey.percentage.2 <- prey.percentage
prey.percentage.2[] <- NA

for(i in 1:nrow(prey.percentage.2)) {
  x <- as.numeric(prey.percentage[i,])
  y <- sum(x)
  x <- x/y*100
  prey.percentage.2[i,] <- x
  print(i)
}

prey.percentage.2 <- as.data.frame(prey.percentage.2)
prey.percentage.2$CarniBinom <- row.names(prey.percentage.2)
row.names(prey.percentage.2) <- NULL

diet.summary <- cbind(as.data.frame(diet.records_perSource),
                      as.data.frame(diet.records_perMethod[c(2:10)]))
diet.summary <- cbind(prey.percentage.2[c(8,1:7)], diet.summary[c(2:30)])
colnames(diet.summary)[1] <- "scientificNameCarni"

# Merge the number of prey species
diet.summary <- merge(diet.summary, prey.consumed, by.y = "scientificNameCarni", all = TRUE)

diet.summary <- cbind(diet.summary[1:103,], total.obvs)
diet.summary[is.na(diet.summary)] <- 0

# Save this table
#write.csv(diet.summary, "./Results/Tables/Diet Summary.csv")

# Get some descriptive statistics for these
colSums(diet.summary[9:39])
colnames(diet.summary)

# Source of data collection
data.sources <- diet.summary[c(9:28)]
data.sources <- melt(data.sources)
data.sources <- data.sources %>% dplyr::group_by(variable) %>% dplyr::summarise(sum = sum(value))
data.sources <- data.sources[data.sources$sum > 100,]

data.sources$variable <- factor(data.sources$variable)
data.sources$variable <- reorder(data.sources$variable, -data.sources$sum)

# Methods of data quantification
data.sources <- diet.summary[c(29:39)]
data.sources <- melt(data.sources)
data.sources <- data.sources %>% dplyr::group_by(variable) %>% dplyr::summarise(sum = sum(value))
data.sources <- data.sources[data.sources$sum > 100,]

data.sources$variable <- factor(data.sources$variable)
data.sources$variable <- reorder(data.sources$variable, -data.sources$sum)



# ---- Figure: Number of records used from each source and methods used in each source ----
diet.method.summary <- carnidiet %>% dplyr::group_by(samplingProtocol, methodQuantification) %>% dplyr::summarise(count = length(methodQuantification))
diet.method.summary <- diet.method.summary[diet.method.summary$count > 100,]

diet.method.summary$samplingProtocol <- factor(diet.method.summary$samplingProtocol)
diet.method.summary$samplingProtocol <- reorder(diet.method.summary$samplingProtocol, -diet.method.summary$count)
diet.method.summary$samplingProtocol <- factor(diet.method.summary$samplingProtocol,
                                            levels(diet.method.summary$samplingProtocol)[c(2,1,4,5,3,6:7)])


# Need to show as percentages
# Get counts for eaxh sample source
sample.count <- diet.method.summary %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarize(total.count = sum(count)) 

# Merge
diet.method.summary <- merge(diet.method.summary, sample.count, by = "samplingProtocol", all.x = TRUE)

# Get percentages
diet.method.summary$Perc <- (diet.method.summary$count/diet.method.summary$total.count)*100

(p1 <- ggplot(diet.method.summary, aes(x = samplingProtocol, y = Perc, fill = methodQuantification)) +
  geom_bar(stat = "identity", colour = "black", lwd = 0.1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  #scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sampling protocol") + ylab("Records (%)") +
  # Add number of 
  scale_fill_brewer(palette = "Greys") +
  annotate("text", x = 1, y = 105, label = "(17266)", size = 2.6) +
    annotate("text", x = 2, y = 105, label = "(6238)", size = 2.6) +
    annotate("text", x = 3, y = 105, label = "(2378)", size = 2.6) +
    annotate("text", x = 4, y = 105, label = "(1266)", size = 2.6) +
    annotate("text", x = 5, y = 105, label = "(936)", size = 2.6) +
    annotate("text", x = 6, y = 105, label = "(123)", size = 2.6) +
    annotate("text", x = 7, y = 105, label = "(108)", size = 2.6) +
    
    theme(axis.text = element_text(size = 8),
          legend.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.title = element_blank(),
          #axis.text.x = element_blank(),
          #axis.title.x = element_blank(),
          legend.position="top",
          legend.key.size =  unit(0.1, "in")))

ggsave("../Final figures/test.tiff", p1, width = 3, height = 3.5)


# ---- Supplementary Figure: Number of records used from each source and methods used in each source ----
diet.method.summary <- carnidiet %>% dplyr::group_by(samplingProtocol, methodQuantification) %>% dplyr::summarise(count = length(methodQuantification))


# Need to show as percentages
# Get counts for eaxh sample source
sample.count <- diet.method.summary %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarize(total.count = sum(count)) 


# Merge
diet.method.summary <- merge(diet.method.summary, sample.count, by = "samplingProtocol", all.x = TRUE)
# Get percentages
diet.method.summary$Perc <- (diet.method.summary$count/diet.method.summary$total.count)*100

diet.method.summary.sum <- diet.method.summary %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarise(sum = sum(count))
diet.method.summary.sum$samplingProtocol <- reorder(diet.method.summary.sum$samplingProtocol, -diet.method.summary.sum$sum)
diet.method.summary.sum <- arrange(diet.method.summary.sum, -diet.method.summary.sum$sum)

levels(sampling.protocol.sum$samplingProtocol)

diet.method.summary$samplingProtocol <- factor(diet.method.summary$samplingProtocol, levels(diet.method.summary.sum$samplingProtocol))


# Get percentages
diet.method.summary$Perc <- (diet.method.summary$count/diet.method.summary$total.count)*100

(p1 <- ggplot(diet.method.summary, aes(x = samplingProtocol, y = Perc, fill = methodQuantification)) +
    geom_bar(stat = "identity", colour = "black", lwd = 0.1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
    #scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Sampling protocol") + ylab("Records (%)") +
    # Add number of 
    scale_fill_brewer(palette = "Paired") +
    annotate("text", x = 1, y = 105, label = "(17490)", size = 2.6) +
    annotate("text", x = 2, y = 105, label = "(6249)", size = 2.6) +
    annotate("text", x = 3, y = 105, label = "(2444)", size = 2.6) +
    annotate("text", x = 4, y = 105, label = "(1301)", size = 2.6) +
    annotate("text", x = 5, y = 105, label = "(936)", size = 2.6) +
    annotate("text", x = 6, y = 105, label = "(188)", size = 2.6) +
    annotate("text", x = 7, y = 105, label = "(128)", size = 2.6) +
    annotate("text", x = 8, y = 105, label = "(119)", size = 2.6) +
    annotate("text", x = 9, y = 105, label = "(77)", size = 2.6) +
    annotate("text", x = 10, y = 105, label = "(48)", size = 2.6) +
    annotate("text", x = 11, y = 105, label = "(40)", size = 2.6) +
    annotate("text", x = 12, y = 105, label = "(26)", size = 2.6) +
    annotate("text", x = 13, y = 105, label = "(17)", size = 2.6) +
    annotate("text", x = 14, y = 105, label = "(15)", size = 2.6) +
    annotate("text", x = 15, y = 105, label = "(11)", size = 2.6) +
    annotate("text", x = 16, y = 105, label = "(10)", size = 2.6) +
    annotate("text", x = 17, y = 105, label = "(8)", size = 2.6) +
    annotate("text", x = 18, y = 105, label = "(6)", size = 2.6) +
    annotate("text", x = 19, y = 105, label = "(6)", size = 2.6) +
    annotate("text", x = 20, y = 105, label = "(2)", size = 2.6) +
    theme(axis.text = element_text(size = 8),
          legend.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.title = element_blank(),
          #axis.text.x = element_blank(),
          #axis.title.x = element_blank(),
          legend.position="top",
          legend.key.size =  unit(0.1, "in")))

ggsave("../Final figures/Figure S9.tiff", p1,
       width = 8, height = 4)


# ---- Figure: Number of records used for each source and the breakdown of taxonomic resolution ----
diet.resolution  <- carnidiet %>% dplyr::group_by(samplingProtocol, taxonRankPrey) %>% dplyr::summarise(count = length(taxonRankPrey))
diet.resolution <- diet.resolution[diet.resolution$count > 100 |
                                     diet.resolution$samplingProtocol == "Scat/Stomach",]
diet.resolution <- diet.resolution[!is.na(diet.resolution$samplingProtocol),]
diet.resolution[is.na(diet.resolution$taxonRankPrey),]$taxonRankPrey <- "Other"


levels(factor(diet.resolution$samplingProtocol))
diet.resolution$samplingProtocol <- factor(diet.resolution$samplingProtocol,
                                        levels(factor(diet.resolution$samplingProtocol))[c(4,2,6,1,3,5,7)])

levels(factor(diet.resolution$taxonRankPrey))
diet.resolution$taxonRank <- factor(diet.resolution$taxonRankPrey,
                                             levels(factor(diet.resolution$taxonRankPrey))[c(4,1,5,2,3,7,6)])

# Need to show as percentages
# Get counts for eaxh sample source
sample.count <- diet.resolution %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarize(total.count = sum(count)) 

# Merge
diet.resolution <- merge(diet.resolution, sample.count, by = "samplingProtocol", all.x = TRUE)

# Get percentages
diet.resolution$Perc <- (diet.resolution$count/diet.resolution$total.count)*100

# Percentage of records in database for each taxon rank out of 29121 records

sum(diet.resolution[diet.resolution$taxonRankPrey == "Species",]$count)/29121
sum(diet.resolution[diet.resolution$taxonRankPrey == "Class",]$count)/29121
sum(diet.resolution[diet.resolution$taxonRankPrey == "Genus",]$count)/29121


# Plot figure
(p <- ggplot(diet.resolution, aes(x = samplingProtocol, y = Perc, fill = taxonRankPrey)) +
  geom_bar(stat = "identity", colour = "black", lwd = 0.1) +
    theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  scale_fill_brewer(palette = "Greys") +
  #scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sampling protocol") + ylab("Records (%)") +
    annotate("text", x = 1, y = 105, label = "(17408)", size = 2.6) +
    annotate("text", x = 2, y = 105, label = "(6242)", size = 2.6) +
    annotate("text", x = 3, y = 105, label = "(2364)", size = 2.6) +
    annotate("text", x = 4, y = 105, label = "(1156)", size = 2.6) +
    annotate("text", x = 5, y = 105, label = "(682)", size = 2.6) +
    annotate("text", x = 6, y = 105, label = "(188)", size = 2.6) +
    annotate("text", x = 7, y = 105, label = "(119)", size = 2.6) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.position="top",
        legend.key.size =  unit(0.1, "in")))

ggsave("../Final figures/test2.tiff", p,
       width = 3, height = 3.5)

# Plot both figures together
(p4 <- ggarrange(p + theme(axis.ticks.x = element_blank()),
                p1 + theme(axis.ticks.x = element_blank()),
                      ncol=2, nrow=1, labels = "auto"))

ggsave("../Final figures/test3.tiff", p4,
       width = 6, height = 3.5)

# Save plot
ggsave("../Final figures/Figure 2/DietResolutionSummary_long.pdf", p4,
       width = 6, height = 3.5, units = "in")


# ---- S figure: ALL RECORDS for sampling protocol and the breakdown of taxonomic resolution ----
diet.resolution  <- carnidiet %>% dplyr::group_by(samplingProtocol, taxonRankPrey) %>% dplyr::summarise(count = length(taxonRankPrey))

# Need to show as percentages
# Get counts for eaxh sample source
sample.count <- diet.resolution %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarize(total.count = sum(count)) 

# Merge
diet.resolution <- merge(diet.resolution, sample.count, by = "samplingProtocol", all.x = TRUE)

# Get percentages
diet.resolution$Perc <- (diet.resolution$count/diet.resolution$total.count)*100

# Percentage of records in database for each taxon rank out of 29121 records

sum(diet.resolution[diet.resolution$taxonRankPrey == "Species",]$count)/29121
sum(diet.resolution[diet.resolution$taxonRankPrey == "Class",]$count)/29121
sum(diet.resolution[diet.resolution$taxonRankPrey == "Genus",]$count)/29121

sampling.protocol.sum <- diet.resolution %>% dplyr::group_by(samplingProtocol) %>% dplyr::summarise(sum = sum(count))
sampling.protocol.sum$samplingProtocol <- reorder(sampling.protocol.sum$samplingProtocol, -sampling.protocol.sum$sum)
sampling.protocol.sum <- arrange(sampling.protocol.sum, -sampling.protocol.sum$sum)

levels(sampling.protocol.sum$samplingProtocol)


diet.resolution$samplingProtocol <- factor(diet.resolution$samplingProtocol, levels(sampling.protocol.sum$samplingProtocol))

# Plot figure
(p <- ggplot(diet.resolution, aes(x = samplingProtocol, y = Perc, fill = taxonRankPrey)) +
    geom_bar(stat = "identity", colour = "black", lwd = 0.1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
    scale_fill_brewer(palette = "Paired") +
    #scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Sampling protocol") + ylab("Records (%)") +
    annotate("text", x = 1, y = 105, label = "(17490)", size = 2.6) +
    annotate("text", x = 2, y = 105, label = "(6249)", size = 2.6) +
    annotate("text", x = 3, y = 105, label = "(2444)", size = 2.6) +
    annotate("text", x = 4, y = 105, label = "(1301)", size = 2.6) +
    annotate("text", x = 5, y = 105, label = "(936)", size = 2.6) +
    annotate("text", x = 6, y = 105, label = "(188)", size = 2.6) +
    annotate("text", x = 7, y = 105, label = "(128)", size = 2.6) +
    annotate("text", x = 8, y = 105, label = "(119)", size = 2.6) +
    annotate("text", x = 9, y = 105, label = "(77)", size = 2.6) +
    annotate("text", x = 10, y = 105, label = "(48)", size = 2.6) +
    annotate("text", x = 11, y = 105, label = "(40)", size = 2.6) +
    annotate("text", x = 12, y = 105, label = "(26)", size = 2.6) +
    annotate("text", x = 13, y = 105, label = "(17)", size = 2.6) +
    annotate("text", x = 14, y = 105, label = "(15)", size = 2.6) +
    annotate("text", x = 15, y = 105, label = "(11)", size = 2.6) +
    annotate("text", x = 16, y = 105, label = "(10)", size = 2.6) +
    annotate("text", x = 17, y = 105, label = "(8)", size = 2.6) +
    annotate("text", x = 18, y = 105, label = "(6)", size = 2.6) +
    annotate("text", x = 19, y = 105, label = "(6)", size = 2.6) +
    annotate("text", x = 20, y = 105, label = "(2)", size = 2.6) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size = 8),
          legend.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.title = element_blank(),
          #axis.text.x = element_blank(),
          #axis.title.x = element_blank(),
          legend.position="top",
          legend.key.size =  unit(0.1, "in")))

ggsave("../Final figures/Figure S8.tiff", p,
       width = 8, height = 3.5)

# Plot both figures together
(p4 <- ggarrange(p + theme(axis.ticks.x = element_blank()),
                 p1 + theme(axis.ticks.x = element_blank()),
                 ncol=2, nrow=1, labels = "auto"))

ggsave("../Final figures/test3.tiff", p4,
       width = 6, height = 3.5)

# Save plot
ggsave("../Final figures/Figure 2/DietResolutionSummary_long.pdf", p4,
       width = 6, height = 3.5, units = "in")



# Summary statistics
table(carnidiet$Sample.Source)
table(carnidiet$Method)
table(carnidiet$TaxonomicPrecision)



# ---- Simulation of prey species richness with random samples ----
# First, descriptive statistics of the number of prey species recorded in diet analyses
x <- carnidiet[carnidiet$TaxonomicPrecision == "Species",]
colnames(x)
prey.info <- x[c(7:13)]

# Proportion of predator-orey records - which orders... 
t <- as.data.frame(table(prey.info$PreyOrder))
t$perc <- (t$Freq/sum(t$Freq))*100

table(prey.info$PreyFamily)
table(prey.info$PreyGenus)

colnames(carnidiet)
x <- carnidiet[c(3:4, 12)]
x$count <- 1
x <- unique(x)

specinteractions <- x %>%
                    dplyr::group_by(Family, CarniBinom) %>%
                    dplyr::summarise(species = sum(count))

carni.diet.subset <- carnidiet[carnidiet$CarniBinom == "Panthera_pardus" | # 54 studies
                                 carnidiet$CarniBinom == "Vulpes_vulpes" | # 66 studies
                                 carnidiet$CarniBinom == "Puma_concolor" | # 40 studies
                                 carnidiet$CarniBinom == "Canis_lupus" | # 59 studies
                                 carnidiet$CarniBinom == "Panthera_leo",] # 25 studies

carni.diet.subset$CarniBinom <- factor(carni.diet.subset$CarniBinom)

species.only <- carni.diet.subset[carni.diet.subset$TaxonomicPrecision == "Species",]
species.only <- species.only[!species.only$PreyBinom == "",]
species.only <- species.only[!is.na(species.only$CarniBinom),]

nlevels(factor(species.only[species.only$CarniBinom == "Panthera_leo",]$PrimaryRef))

accumulation.matrix <- matrix(ncol = 66, nrow = length(levels(carni.diet.subset$CarniBinom)))
row.names(accumulation.matrix) <- levels(carni.diet.subset$CarniBinom)
colnames(accumulation.matrix)[1:66] <- 1:66
data <- data.frame()

for(a in 1:100) {
for(i in 1:length(levels(carni.diet.subset$CarniBinom))) {
  
  carnivore <- levels(carni.diet.subset$CarniBinom)[i]
  
  x <- species.only[species.only$CarniBinom == carnivore,]
  x <- x[!is.na(x$CarniBinom),]
  
  studies <- levels(factor(x$PrimaryRef))
  species.list <- data.frame(ncol = 1)
  colnames(species.list) <- "Prey"
  print("NEW SPECIES")
  
  if(length(studies) > 0){
    for(j in 1:length(studies)) {
      
      stud <- sample(studies,1)
      y <- x[x$PrimaryRef == stud,]
      prey <- as.data.frame(y$PreyBinom)
      colnames(prey) <- "Prey"
      
      species.list <- rbind(species.list, prey)
      species.list <- unique(species.list)
      
      z <- nrow(species.list) - 1
      accumulation.matrix[i,j] <- z
      remove <- as.numeric(which(studies %in% stud))
      studies <- studies[-c(remove)]
      print(j)
    }
  }
  else {}
}

resultsss <- melt(accumulation.matrix)
resultsss <- resultsss[!is.na(resultsss$Var2),] 
resultsss$Repetition <- a
data <- rbind(resultsss,data)
print(paste(a,"YAAAAAAAAAY"))

}

# Get average and standard error around following the 100 samples
data.updated <- data %>% 
                dplyr::group_by(Var1, Var2) %>% 
                dplyr::summarise(mean.sr = mean(value),
                          sd = sd(value))
data.updated$Var1 <- revalue(data.updated$Var1, c("Vulpes_vulpes" = "V. vulpes",
                                                  "Canis_lupus" = "C. lupus",
                                                  "Panthera_leo" = "P. leo",
                                                  "Panthera_pardus" = "P. pardus",
                                                  "Puma_concolor" = "P. concolor"))

data$Var1 <- revalue(data$Var1, c("Vulpes_vulpes" = "V. vulpes",
                                                  "Canis_lupus" = "C. lupus",
                                                  "Panthera_leo" = "P. leo",
                                                  "Panthera_pardus" = "P. pardus",
                                                  "Puma_concolor" = "P. concolor"))

# Add a common grouping variable for species and the repetition
data$rep.spec <- paste0(data$Var1,data$Repetition)

# ---- Figure: Simulation of top 5 species with the highest prey species richness recorded ----
p.plot <- ggplot() +
  geom_line(data = data, aes(x = Var2, y = value, group = rep.spec, colour = Var1), alpha = .1) +
  geom_line(data = data.updated, aes(x = Var2, y = mean.sr, group = Var1, colour = Var1)) +
  scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(title = "Species")) +
  theme_bw() + labs(x = "Number of studies", y = "Prey species richness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(face = "italic"))

ggsave(filename = paste0("./Results/Figures/ppinteractions_simulatios.tiff"), plot = p.plot, device = NULL, path = NULL,
       scale = 1, width = 16, height = 8, units = c("cm"),
       dpi = 300)    

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- GRAVEYARD SHIFT Current not in use ----
# ---- Extra predator-prey code ----

bodymass <- phylacine[c(1,12)]

predprey <- carnidiet[c(4,12,19,51:52)]
colnames(predprey)[1] <- "Binomial.1.2"
predprey <- merge(predprey, bodymass, by.y = "Binomial.1.2")

colnames(predprey)[1:6] <- c("Predator","Binomial.1.2","Interaction","Source","Method","Pred.Mass.g")

predprey <- merge(predprey, bodymass, by.y = "Binomial.1.2")

colnames(predprey)[1] <- "Prey"
colnames(predprey)[7] <- "Prey.Mass.g"

colnames(predprey)

predprey <- predprey[c(2,6,1,7,3,4,5)]

ggplot(data = predprey, aes(x = log10(Pred.Mass.g), y = log10(Prey.Mass.g), colour = log10(Interaction))) +
  geom_point() +
  geom_vline(xintercept = 4.332438, linetype = 2) +
  facet_wrap(.~Source) +
  geom_smooth(method = "loess") + 
  theme_bw()

predprey2 <- predprey %>% dplyr::group_by(Predator,Prey,Pred.Mass.g,Prey.Mass.g) %>% dplyr::summarise(freq = length(Prey))  

ggplot(data = predprey2, aes(x = log10(Pred.Mass.g), y = log10(Prey.Mass.g), colour = log10(freq + 1), size = log10(freq + 1))) +
  geom_point() +
  geom_vline(xintercept = 4.332438, linetype = 2) +
  geom_smooth(method = "loess") + 
  theme_bw()

# Tidying up prey species taxonomy 

# Which prey species do not merge with Phylacine
split <- carnidiet[carnidiet$TaxonomicPrecision == "Species",]
which(!carnidiet$PreyBinom %in% levels(phylacine$Binomial.1.2))

bad.taxonomy <- split[which(!split$PreyBinom %in% levels(phylacine$Binomial.1.2)),]

write.csv(bad.taxonomy, "TidyingBadTaxonomy.csv")