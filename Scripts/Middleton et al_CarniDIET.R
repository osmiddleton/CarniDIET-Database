#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- SCRIPT OVERVIEW ----
# 
#   Chapter 3: CarniDIET: establishing a macroecologcal dataset on the diets of Carnivora
#
#   The purpose of this script is to summarise the data collected within CarniDIET. This will involve evaluating several
#   key aspects:-
#         
#   - Taxonomic coverage:
#        > Percentage of species in a Family with at least one dietary study completed
#        > Geographic range of missing species (which biomes have the most species unstudied?)
#        > Number of studies per species, highlighting those species which make up the majority of the total studies (50%)
#   
#   - Spatio-temporal coverage:
#        > Spatio-temporal coverage (latitude~Years studied) (i.e. PREDICTS)
#        > Global maps of study distributions
#   
#   - Dietary coverage:
#        > Number of species
#        > Percentage recorded to different taxonomic resolutions
#        
#   - Methdological information
#        > Percentage of records made using different methods
# 


#  Terminology:
# 
#  WHAT MAKES A PAPER AND A STUDY?
# 
#     - Paper: Primary reference from where dietary information was collected.
#     - Study: Finest coverage of dietary information from a paper (No.Studies from paper => 1 paper)
#               > Multiple years
#               > Annual
#               > Seasonal
#               > Finer
#       Note: This is also influenced by the sex/age of the individuals studied in 
#       the papers (i.e. male/female or adult/juvenile will be split into multiple studies within a paper).
#       
#   Taxonomy:
#   This study will follow the taxonomy of the Phylacine v1.2 database due to the purpose of being built to be a standardised
#   database for macroecolgoical research on mammals in the Late Quaternary
#————————————————————————————————————————————————————————————————————————————————————————————####

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- BOOK-KEEPING ----
#————————————————————————————————————————————————————————————————————————————————————————————####
#---- Functions ----
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#---- Colour schemes ----
colours = c("#99CCFF","#3399FF", "#003366", # Blues
            "#CCFFCC","#66CC00", "#336600", # Greens
            "#FFCCCC","#FF99CC", "#660033", # Pink/Red
            "#FFCC99","#FF9933", "#996666", # Oranges
            "#9999FF","#9933FF", "#333399") # Purple

#---- Load packages ----
library(tidyverse)
library(plyr)
library(ggrepel)
library(sp)
library(sf)
library(viridis)
library(rgdal)
library(rgeos)
library(maptools)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(doParallel)
library(effects)
library(raster)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library("vegan")
library(FD)
library(factoextra)
library(pscl)
library(boot)
library(mapproj)

#---- Load data ----
# Load Phylacine - for name standardisation:
#phylacine <- read.csv("./Data/TraitData/phylacine/Trait_data.csv", header = TRUE)

# Load CarniDIET
carnidiet <- read.csv("./Version 1.0/CarniDIET 1.0.csv")
#carnidiet <- carnidiet[c(2:75)]

# Load potential species list
# Includes species selected as primary consumers from MammalDIET (Kissling et al., 2014)
cd.wos.hits <- read.csv("./Version 1.0/Supplementary data/CarniDiet_PotentialSpecies.csv", header = TRUE)

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- CarniDIET Metadata ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- 1. Carnivore species information ----

# How many potential ...?
nrow(cd.wos.hits) # 213 species
levels(cd.wos.hits$Order) # 9 orders
levels(cd.wos.hits$Family) # 23 families

# How many ... in CarniDIET
levels(carnidiet$CarniBinom) # 104 species
levels(carnidiet$Order) # 5 orders
levels(carnidiet$Family)  # 15 families

# How many levels of:
levels(carnidiet$Age) # 4 ages/life-stages
levels(carnidiet$Sex) # Both sexes and both/unknown

# ---- 2. Diet resolution information -----

# How many levels for food-type variables:
# the '-1' accounts for "" entries which could be treated as NA's...
nlevels(carnidiet$Food.type) - 1 # 16 Food types
nlevels(carnidiet$PreyOrder) - 1 # 62 prey orders
nlevels(carnidiet$PreyFamily) - 1 # 130 prey families
nlevels(carnidiet$PreyGenus) - 1 # 490 prey genera
nlevels(carnidiet$PreySpecies) - 1 # 735 prey species
nlevels(carnidiet$PreyBinom) - 1 # 866 prey binomials
nlevels(carnidiet$Common.name) - 1 # 1649 common names

# Categories of taxonomic precision
nlevels(carnidiet$TaxonomicPrecision) - 1 # 11 categories of taxonomic precision

# Range of taxonomic importance
range(carnidiet$QuantificationImportance) # 0.005% to 100%

# Range of sample sizes used for quantifying diet compositions
range(carnidiet[!is.na(carnidiet$Scat.Stomach.Tissue),]$Scat.Stomach.Tissue) # 1 to 11,478
range(carnidiet[!is.na(carnidiet$Prey.Items.Kill.Sample.Size),]$Prey.Items.Kill.Sample.Size) # 1 to 23,824

# ---- 3. Temporal information ----

# Range in study years
range(carnidiet[!is.na(carnidiet$Start.Year),]$Start.Year) # 1933 to 2017
range(carnidiet[!is.na(carnidiet$End.Year),]$End.Year) # 1946-2017

# Months
levels(carnidiet[!is.na(carnidiet$Start.Month),]$Start.Month)
levels(carnidiet[!is.na(carnidiet$End.Month),]$End.Month)
# All months represented in start and month of study

# Day of the month the study started
levels(factor(carnidiet[!is.na(carnidiet$Start.Day),]$Start.Day))
levels(factor(carnidiet[!is.na(carnidiet$End.Day),]$End.Day))

# Seasons - these are pretty crude
levels(carnidiet$Season)

# ----- 3. Methods information ----
levels(carnidiet$Method)
levels(carnidiet$Sample.Source)

# ----- 4. Spatial information -----
# Altitude
range(carnidiet[!is.na(carnidiet$AltitudeMinimum),]$AltitudeMinimum) # 0 - 4543m
range(carnidiet[!is.na(carnidiet$AltitudeMaximum),]$AltitudeMaximum) # 0 - 8156m

# Latitude and longitude
range(carnidiet[!is.na(carnidiet$Latitude.Centroid.Mean.x),]$Latitude.Centroid.Mean.x) # -55S to 81N
range(carnidiet[!is.na(carnidiet$Longitude.Centroid.Mean.x),]$Longitude.Centroid.Mean.x) # -169W to 178E
levels(carnidiet$Coordinates.Source) # Some extra ones

# Descriptions of study area
levels(carnidiet$SiteName) # 791 sites
levels(carnidiet$Region) # 410 Regions
levels(carnidiet$Country) # 107 "countries"

# Range of study area size
levels(carnidiet$Study.Area.Size.km2)

# Remove descriptive study area sizes
study.area <- carnidiet %>%
  filter(!Study.Area.Size.km2 %in% c("50km radius",
                                   "50 - 100",
                                   "3km radius")) %>% 
  droplevels()
# Convert to numeric
vec <- unique(as.numeric(as.character(study.area[!is.na(study.area$Study.Area.Size.km2),]$Study.Area.Size.km2)))
# Distribution on log scale
hist(log10(vec[!is.na(vec)]))

# Spread
summary(log10(vec[!is.na(vec)]))
range(vec[!is.na(vec)]) # 0.03 to 100000km2


# ---- 4. Bibliographic information ----
levels(carnidiet$PR.Author)
range(carnidiet$PR.Year)
levels(carnidiet$PR.Title)
levels(carnidiet$PR.Journal)
levels(carnidiet$PrimaryRef)
levels(carnidiet$CollectionRef)
levels(carnidiet$Collector)


#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- CARNIVOROUS MAMMAL DIET COVERAGE  ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Family representation in dietary studies -----
# Number of POTENTIAL species in CarniDIET (by family)
family.species <- cd.wos.hits %>% dplyr::group_by(Family) %>% dplyr::summarize(No.Species = length(Family))
family.species

# Number of ACTUAL species in CarniDIET (by family)
family.CD <- carnidiet[c(3:4)]
family.CD <- unique(family.CD)
family.CD <- family.CD %>% dplyr::group_by(Family) %>% dplyr::summarise(No.Species = length(Family))
family.CD

# Number of WoS hits per Family
family.wos.hits <- cd.wos.hits %>% dplyr::group_by(Family) %>% dplyr::summarize(No.Hits = sum(WoS.Hits))
family.wos.hits

# Number of papers included in CarniDIET (by Family)
family.papers <- carnidiet[c(3,44)]
family.papers <- unique(family.papers)
family.papers <- family.papers %>% dplyr::group_by(Family) %>% dplyr::summarise(No.papers = length(Family))
family.papers

# Median and max of number of studies per species (per Family)
family.species.papers <- carnidiet[c(3:4,44)]
family.species.papers <- unique(family.species.papers)
family.species.papers <- family.species.papers %>% dplyr::group_by(Family,CarniBinom) %>% dplyr::summarise(studiesperspecies = length(CarniBinom))
family.species.papers <- family.species.papers %>% dplyr::group_by(Family) %>% dplyr::summarise(median.papers = median(studiesperspecies),
                                                                                                max.papers = max(studiesperspecies))
family.species.papers

# Building Family-level table
family.info <- merge(family.species,family.CD, by = "Family", all = TRUE)
family.info$Perc.Families <- (family.info$No.Species.y/family.info$No.Species.x)*100
family.info <- merge(family.info, family.wos.hits, by = "Family", all = TRUE)
family.info <- merge(family.info, family.papers, by = "Family", all = TRUE)
family.info$Hits_Papers_Ratio <- (family.info$No.papers/family.info$No.Hits)*100
family.info <- merge(family.info, family.species.papers, by = "Family", all = TRUE)
family.info[is.na(family.info)] <- 0
family.info

# How many papers were used compared to actual hits?
sum(family.info$No.Hits)
sum(family.info$No.papers)
median(family.info$Hits_Papers_Ratio)

# ---- Genus representation in dietary studies -----
# Number of potential species in each genus
genus.species <- cd.wos.hits %>% dplyr::group_by(Family, Genus_V_1.2) %>%
              dplyr::summarize(No.Species = length(Genus_V_1.2))

# Number of actual species included in carnidiet
diet.2 <- carnidiet %>% separate(CarniBinom, into = paste(c("Genus_V_1.2","Species"), sep = "_"))

genus.CD <- diet.2[c(3:4,5)]
genus.CD <- unique(genus.CD)
genus.CD <- genus.CD %>% dplyr::group_by(Family,Genus_V_1.2) %>% dplyr::summarise(No.Species = length(Genus_V_1.2))

# Number of hits per genus
genus.wos.hits <- cd.wos.hits %>% dplyr::group_by(Family,Genus_V_1.2) %>% dplyr::summarize(No.Hits = sum(WoS.Hits))

# Number of actual papers per genus
genus.papers <- diet.2[c(3:4,45)]
genus.papers <- unique(genus.papers)
genus.papers <- genus.papers %>% dplyr::group_by(Family,Genus_V_1.2) %>% dplyr::summarise(No.papers = length(Genus_V_1.2))

# Median and max of number of studies per species per Genus
genus.species.papers <- diet.2[c(3:5,45)]
genus.species.papers <- unique(genus.species.papers)
genus.species.papers <- genus.species.papers %>% dplyr::group_by(Family,Genus_V_1.2,Species) %>% dplyr::summarise(studiesperspecies = length(Species))
genus.species.papers <- genus.species.papers %>% dplyr::group_by(Family,Genus_V_1.2) %>% dplyr::summarise(median.papers = median(studiesperspecies),
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
genus.info

# ---- Table: Family & genus level representation in dietary studies ----
write.csv(family.info, "./Results/Tables/Family.Information.Summary_UPDATE.csv")
write.csv(genus.info, "./Results/Tables/Genus.Information.Summary_UPDATE.csv")

# ---- Figure: Family-level dietary representation ----
# Family info
# Find 75%,50%,25,10% slopes
# line: y = mx + c
# m = (y2 - y1/x2 - x1)

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
                                                                fill = Family, size = log10(No.papers)), label = Family) +
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
  geom_point(aes(fill = Family), alpha = 0.5, pch=21, colour = "black") +
  scale_fill_manual(values = colours, guide = "none") +
  geom_text_repel(aes(label = Family),
                  size = 2, colour = "black", force = 4) +
  xlim(0,1.8) + ylim(0,1.75)+
  theme_bw() + xlab("log10(number of mammal-consumers in Family)") +
    ylab("log10(number of mammal-consumers with diet study)") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()))

ggsave(filename = paste0("./Results/Figures/PercentFamiliesStudied.tiff"), plot = plot, device = NULL, path = NULL,
       scale = 1, width = 10.5, height = 10.5, units = c("cm"),
       dpi = 300)

# ---- Figure: Where are potential, but missing, species located? ----
# This section produces a map of the potential species richness of a cell which have not been studied.

# Load Phylacine species range rasters 
current.range.full <- read.csv("./Data/GIS/Phylacine ranges/Current_Species_Matrix.csv", header = TRUE)

# Potential species richness
# First, which species do not match those in phylacine range data
which(!levels(cd.wos.hits$Bin.) %in% colnames(current.range.full[c(9:5839)]))

carnivore.potential <- which(colnames(current.range.full[c(9:5839)]) %in% cd.wos.hits$Bin.)
potential.range.full <- cbind(current.range.full[c(1:8)],current.range.full[c(9:5839)][c(carnivore.potential)])

pixel.id <- potential.range.full[c(2:8)]
species.ranges <- as.matrix(potential.range.full[c(9:212)])
species.richness <- rowSums(species.ranges)

row.names(species.ranges) <- paste0("comm",1:51120)

test <- melt(species.ranges)
colnames(test)[1] <- "Community"

# Merge to get pixel information for all communities in the melted dataframe
potential.range.full <- merge(test,pixel.id, by.y = "Community") # This may take some time

# Group by community to get species richness per cell
potential.SR <- potential.range.full %>% dplyr::group_by(Community,x,y) %>% dplyr::summarise(SR = sum(value))
potential.SR <- merge(potential.SR, pixel.id, by = "Community")

# Actual species richness
carnivore.actual <- which(colnames(current.range.full[c(9:5839)]) %in% levels(carnidiet$CarniBinom))
actual.range.full <- cbind(current.range.full[c(1:8)],current.range.full[c(9:5839)][c(carnivore.actual)])

pixel.id <- actual.range.full[c(2:8)]
species.ranges <- as.matrix(actual.range.full[c(9:112)])

row.names(species.ranges) <- paste0("comm",1:51120)
test <- melt(species.ranges)
colnames(test)[1] <- "Community"

# Merge to get pixel information for all communities in the melted dataframe
actual.range.full <- merge(test,pixel.id, by.y = "Community") # This will also likely take some time

# Group by community to get species richness per cell
actual.SR <- actual.range.full %>% dplyr::group_by(Community,x,y) %>% dplyr::summarise(SR = sum(value))
actual.SR <- merge(actual.SR, pixel.id, by = "Community")

# Calculate species richness of species not studied
deficit.SR <- actual.SR
deficit.SR$SR <- potential.SR$SR - actual.SR$SR

potential.SR$data <- "Potential"
actual.SR$data <- "Actual"
deficit.SR$data <- "Deficit"

# Merge all for facet'ing
all.Sr <- rbind(potential.SR,actual.SR, deficit.SR)

# Load in biomes using the sf package
biomes <- st_read("./Data/GIS/ecoregions_biomes/biomes/global_biomes_behrmann.shp")

# Plot this
(deficit <- ggplot() +
  geom_raster(data = deficit.SR[!is.na(deficit.SR$Continent),], aes(x = x.x, y = y.y, fill = SR)) +
  geom_sf(data = biomes, fill = NA, colour = "black") +
  scale_fill_viridis(option = "inferno", name = "SR",
                     limits = c(0,10), breaks = c(0, 2, 4, 6, 8, 10)) + 
  #facet_wrap(.~data, ncol = 1) +
  theme_bw() + xlab("") + ylab("") + 
  theme(legend.position = "right") +
  guides(fill = guide_colourbar(barwidth = .5)) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank()))

# Save this plot 
ggsave(filename = paste0("./Results/Figures/Species richnes not studied_UPDATED.tiff"), plot = deficit, device = NULL, path = NULL,
       scale = 1, width = 20, height = 10, units = c("cm"),
       dpi = 200)

# ---- Figure: Bias by conservation status? ----
colnames(carnidiet)
x <- carnidiet[c(3,4,57)]
x <- unique(x)

species.included.paper <- x %>% dplyr::group_by(Family,CarniBinom) %>% dplyr::summarise(papers = length(CarniBinom))
colnames(species.included.paper)[2] <- "Bin."
species.included.paper$inclusion <- "Included"
y <- merge(species.included.paper, cd.wos.hits[c(7,10)], by = "Bin.", all = TRUE)
y <- y[c(5,1,3,4)]
y$papers[is.na(y$papers)] <- 0
y$inclusion[is.na(y$inclusion)] <- "NotIncluded"

colnames(y)[2] <- "CarniBinom"

# Merge ICUN statuses to new datasets
colnames(phylacine)[1] <- "CarniBinom"
y <- merge(y, phylacine[c(1,18)])

# find length of species included/not included
nrow(y[y$inclusion == "Included",])
nrow(y[y$inclusion == "NotIncluded",])

# Sum by IUCN status
pred.iucn.summary <- y %>% dplyr::group_by(inclusion, IUCN.Status.1.2) %>% dplyr::summarise(species = length(CarniBinom))

# Make species missing minus
x <- pred.iucn.summary[pred.iucn.summary$inclusion == "NotIncluded",]
x$perc <- (x$species/sum(x$species))*100
x$species = -(x$species)
x$perc = -(x$perc)
y <- pred.iucn.summary[pred.iucn.summary$inclusion == "Included",]
y$perc <- (y$species/sum(y$species))*100
pred.iucn.summary <- rbind(x,y)

pred.iucn.summary$IUCN.Status.1.2 <- ordered(pred.iucn.summary$IUCN.Status.1.2, levels = c("DD","CR","EN",
                                                                                           "VU","NT","LC"))
(iucn.plot <- ggplot(pred.iucn.summary[!is.na(pred.iucn.summary$IUCN.Status.1.2),], aes(x = IUCN.Status.1.2, y = perc, fill = inclusion)) +
  geom_bar(stat = "identity") + theme_bw() + coord_flip() +
  scale_fill_manual(values = c("#000000","#666666"), guide = "none") + ylab("Number of species (%)") +
  xlab("IUCN Red List Category") + 
  geom_hline(yintercept = 0, linetype = 2, colour = "#CCCCCC") +
  annotate("text",x = "DD", y = -40, label = "Absent", size = 5, colour = "#666666") +
  annotate("text",x = "DD", y = 40, label = "Present", size = 5, colour = "#000000") +
  scale_y_continuous(limits = c(-80,80),
                     breaks = c(-80, -60, -40, -20, 0 , 20, 40, 60, 80),
                     label = c("80","60", "40", "20", "0", "20", "40", "60", "80")))

ggsave(filename = paste0("./Results/Figures/IUCN Categories missed and present.tiff"), plot = iucn.plot, device = NULL, path = NULL,
       scale = 1, width = 10, height = 10, units = c("cm"),
       dpi = 200)

# ---- Figure: Best represented species in database ----
# How many papers are there?
nlevels(carnidiet$PrimaryRef) # 749 primary references

# Get the number of papers focussing on a species
carnivore.summary <- carnidiet %>%
                     dplyr::group_by(Family, CarniBinom) %>%
                     dplyr::summarise(papers = length(unique(PrimaryRef))) # Number of papers
 
# Number of studies
carnivore.summary.2 <- carnidiet %>% 
                       dplyr::group_by(Family, CarniBinom, Age, Sex, Scat.Stomach.Tissue,Prey.Items.Kill.Sample.Size,
                                       PR.Title, Country, Region, SiteName, AltitudeMaximum,
                                       Start.Year,End.Year, Start.Month,End.Month, Season, Sample.Source, Method, ECO_NAME,
                                       BIOME, Latitude.Centroid.Mean.x, Longitude.Centroid.Mean.x) %>% 
                       dplyr::summarise(x = length(CarniBinom))

carnivore.summary.3 <- carnivore.summary.2 %>%
  dplyr::group_by(Family, CarniBinom) %>%
  dplyr::summarise(studies = length(CarniBinom)) # Number of papers

carnivore.summary.master <- merge(carnivore.summary, carnivore.summary.3, by = "CarniBinom")

# How many studies are there
sum(carnivore.summary.master$papers) # 950 - although this includes papers with multiple carnivores
length(unique(carnidiet$PrimaryRef)) # 749 - This is individual papers

# How many studies are there
sum(carnivore.summary.master$studies) # 3161 - This is all individual studies

3161/749 # On average, there are 4 studies per paper

# 475 papers = 50% papers
# 713 = 75% papers

#---- Figure: Phylogenetic tree and number of studies per species ----

# 1. Histogram of studies per species
(hist <- ggplot(data = carnivore.summary.master, aes(x = log10(studies))) +
  geom_histogram(binwidth = 0.2, colour = "black", fill = "lightgrey") +
  scale_x_continuous(labels = c("0" = "1",
                                "1" = "10",
                                "2" = "100",
                                "3" = "1000")) +
  labs(x = "Number of studies", y = "Number of species") +
  theme_bw() +
  theme(panel.grid = element_blank()))

ggsave("./Results/Figures/Final figures/Figure 1/Histogram.pdf", hist,
       width = 2, height = 2, units = "in")

# Quickly make master figure for the top 20 species
colnames(carnivore.summary.master)[1] <- "CarniBinom"
carnivore.summary.master$CarniBinom <- reorder(carnivore.summary.master$CarniBinom, carnivore.summary.master$studies)

x <- carnivore.summary.master[order(carnivore.summary.master$studies, decreasing = TRUE),]
y <- x[1:20,]
y$CarniBinom <- factor(y$CarniBinom)
levels(y$CarniBinom)


# How many species need to be included to get 50% of the studies?
(sum(x[1:9,]$studies)/3161)*100 #  9 species includes just over 50% studies
(sum(x[1:24,]$studies)/3161)*100 # 24 species includes 75% studies

# Find percentage lines
(sum(x[1:1,]$studies)/3161)*100 # 14.8% studies
(sum(x[1:5,]$studies)/3161)*100 # 40.5% studies
(sum(x[1:10,]$studies)/3161)*100 # 62.9% studies

# Rename factor levels to make a nicer figure.
y$CarniBinom <- revalue(y$CarniBinom, c("Vulpes_vulpes" = "V. vulpes","Lynx_lynx" = "L. lynx",
                                        "Canis_lupus" = "C. lupus","Felis_silvestris" = "F. silvestris",
                                        "Panthera_pardus" = "P. pardus","Neovison_vison" = "N. vison",
                                        "Puma_concolor" = "P. concolor","Martes_foina" = "M. foina",
                                        "Canis_latrans" = "C. latrans","Genetta_genetta" = "G. genetta",
                                        "Ursus_arctos" = "U. arctos","Canis_aureus" = "C. aureus",
                                        "Martes_martes" = "M. martes","Mustela_erminea" = "M. erminea",
                                        "Panthera_tigris" = "P. tigris","Panthera_uncia" = "P. uncia",
                                        "Martes_americana" = "M. americana","Lynx_rufus" = "L. rufus",
                                        "Martes_pennanti" = "M. pennanti","Crocuta_crocuta" = "C. crocuta"))


(p10 <- ggplot(data = y[1:20,], aes(x = reorder(CarniBinom, studies), y = studies, fill = papers)) + # 20 species refers to 66% all papers
  geom_bar(stat = "identity") +
  #annotate("rect", xmin = 19.5, xmax = 20.5, ymin = -Inf, ymax =Inf, fill = "red") +
  #annotate("rect", xmin = 15.5, xmax = 19.5, ymin = -Inf, ymax =Inf, fill = "orange") +
  #annotate("rect", xmin = 7.5, xmax = 15.5, ymin = -Inf, ymax =Inf, fill = "yellow") +
  coord_flip() +
  xlab("") + ylab("Number of studies") +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Papers") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 19.5, colour = "black", linetype = 2, alpha = 0.8) + # 10% studies
  geom_vline(xintercept = 14.5, colour = "black", linetype = 2, alpha = 0.8) + # 30% studies
  geom_vline(xintercept = 5.5, colour = "black", linetype = 2, alpha = 0.8) + # 50% studies 
  annotate("text", x = 19, y = 380, label = "10%", size = 3) +
  annotate("text", x = 14, y = 380, label = "30%", size = 3) +
  annotate("text", x = 5, y = 380, label = "50%", size = 3) +
  ylim(0,400) + theme(axis.text.y = element_text(face = "italic")))

ggsave(filename = paste0("./Results/Figures/Primary species in dataset.tiff"), plot = p10, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c("cm"),
       dpi = 200)

pacman::p_load(ggplot2,
               dplyr,
               stringr,
               gridExtra,
               viridisLite,
               raster,
               rgdal,
               maptools,
               ape,
               ggtree, update = F)


forest <- read.nexus("./Data/Phylogenies/Complete_phylogeny.nex")
names(forest) <- NULL
set.seed(42)
forest <- forest[sample(1:1000, 30)]

# I will be colouring by number of papers
# Current species list
sp <- carnivore.summary.master$Binomial.1.2

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

# Group tree by extant species for marking them later
tree <- groupOTU(tree,
                 tree$tip.label %in% sp)

# Turn tree to dataframe
tree <- fortify(tree)

# Change columns name in carnivore summary
colnames(carnivore.summary.master)[1] <- "Binomial.1.2"

# Add trait info
tree <- left_join(tree, carnivore.summary.master, by = c("label" = "Binomial.1.2"))

# Reverse age scale
tree$x <- tree$x - max(tree$x)

# Prepare studies based legend
#mass.breaks <- marsupials$Mass.g/1000
#mass.breaks <- ceiling(c(min(mass.breaks), median(mass.breaks), max(mass.breaks)))

# Prepare IUCN status color legend. We use a modified version of IUCN's official color palette
# status.colors <- c("EP" = "#87421F", "EX" = "#8F47B3", "EW" = "#8F47B3", 
#                    "CR" = "#D81E05", "EN" = "#FC7F3F", "VU" = "#F9E814", 
#                    "NT" = "#CCE226", "LC" = "#60C659", "DD" = "#D1D1C6")

# Convert multiPhylo to data.frame
pruned.forest <- fortify(pruned.forest)

# Reverse ages to show time BP
pruned.forest$x <- pruned.forest$x - max(pruned.forest$x)

# Plot forest of uncertainty (based on only 30 out of 1000 trees for speed)
# p.tree <- ggplot(pruned.forest) +
#   geom_tree(col = "lightblue", alpha = .3, multiPhylo = TRUE, layout = "circular") +
#   theme_tree2() +
#   scale_x_continuous("Time BP (Ma)",
#                      limits = c(min(tree$x), 23), breaks = seq(-50, 0, 10),
#                      label = abs(seq(-50, 0, 10)))

# Create new name binomial names
foo <- data.frame(do.call('rbind', strsplit(as.character(tree$label),'_',fixed=TRUE)))
tree$Genus <- substring(foo$X1,1,1)
tree$Species <- foo$X2
tree$label2 <- paste(tree$Genus, tree$Species, sep = ". ")

# What happens if you change the branch length?
tree$branch.length <- tree$branch.length/2

colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

(p <- ggtree(tree, layout='circular', ladderize = FALSE) +
    geom_tiplab(size=1.5, aes(label = label2, parse = T, size= 3.5, angle=angle), offset = 10) +
    geom_tippoint(aes(x, y, fill = Family.x, size = studies), pch = 21, colour = "black", alpha = 0.75) +
    scale_fill_manual(values = colours) +
    xlim(-160, 30) +
    theme(panel.background = element_blank(),
          legend.position = "none"))

ggsave("./Results/Figures/Final figures/Figure 1/Phylo tree.pdf", p,
       width = 4, height = 4, units = "in")


#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- || SPATIAL AND TEMPORAL DISTRIBUTION OF DIETARY STUDIES||  ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Summarising data ----
# Find the number of studies
levels(carnidiet$Latitude.Centroid.Mean.x)
str(carnidiet$Latitude.Centroid.Mean.x)
as.numeric(as.character(carnidiet$Latitude.Centroid.Mean.x))

temporal.diet <- carnidiet %>% 
  dplyr::group_by(Family, CarniBinom, Age, Sex, Scat.Stomach.Tissue,Prey.Items.Kill.Sample.Size,
                  PR.Title, Country, Region, SiteName, AltitudeMaximum,
                  Start.Year,End.Year, Start.Month,End.Month, Season, Sample.Source, Method, ECO_NAME,
                  BIOME, Latitude.Centroid.Mean.x, Longitude.Centroid.Mean.x) %>% 
  dplyr::summarise(x = length(CarniBinom))
nrow(temporal.diet)

# Find the midpoint of each study
temporal.diet$Mid.Year <- (temporal.diet$Start.Year + temporal.diet$End.Year)/2
temporal.diet <- temporal.diet[!temporal.diet$BIOME > 14,]
temporal.diet$BIOME <- as.factor(temporal.diet$BIOME)
temporal.diet$Latitude.Centroid.Mean.x <- as.numeric(as.character(temporal.diet$Latitude.Centroid.Mean.x))

temporal.diet <- temporal.diet[!is.na(temporal.diet$CarniBinom),]
temporal.diet <- temporal.diet[!temporal.diet$Latitude.Centroid.Mean.x == "#N/A",]
temporal.diet$Latitude.Centroid.Mean.x <- as.numeric(as.character(temporal.diet$Latitude.Centroid.Mean.x))
str(temporal.diet)

# Max year 
min(temporal.diet[!is.na(temporal.diet$Start.Year),]$Start.Year)
max(temporal.diet[!is.na(temporal.diet$End.Year),]$End.Year)
2017-1933 # 84 years


# ---- Figure: Distribution of data in time and space ----

# Find median and IQR for study latitude and year
mid.year <- temporal.diet$Mid.Year
median.year <- median(mid.year[!is.na(mid.year)])
quantile(mid.year[!is.na(mid.year)])

Latitude.Centroid.Mean.x <- temporal.diet$Latitude.Centroid.Mean.x
median.lat <- median(Latitude.Centroid.Mean.x[!is.na(Latitude.Centroid.Mean.x)])
quantile(abs(Latitude.Centroid.Mean.x[!is.na(Latitude.Centroid.Mean.x)]))

# Create shading dataframe
shading <- data.frame(xmin = c(0,12.04167,12.04167),
                      xmax = c(80,49.05833,49.05833),
                      ymin = c(1994,1930,1994),
                      ymax = c(2006,2010,2006),
                      rect = c("A","A","B"))

# Get colour scheme properly:

# Remove 9th and 12th colour from previous colour scheme to account for the fact that Lorisidae and 
# Phylloestomidae are not shown here as they do not have coords or year of study... 
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#CC6699","#FF99FF","#9999CC","#3399FF","#CCFFFF","#003333")

# Plot Predicts-style figure which shows the year and location of studies
p1 <- ggplot() +
  annotate("rect",xmin =-Inf, xmax = Inf, ymin = 1995, ymax =2008, alpha = .2) +
  annotate("rect",xmin = 30.857, xmax = 48,ymin = -Inf, ymax =Inf, alpha = .2) +
  annotate("rect",xmin = 30.857, xmax = 48,ymin = 1995, ymax =2008, alpha = .2) +
  ylab("Year") + xlab("|Latitude|") +
  coord_flip() +
  scale_y_continuous(breaks = c(1940,1960,1980,2000,2020)) + 
  # geom_pointrange(data = temporal.diet, 
  #                 aes(y = Mid.Year, ymin = Start.Year, ymax = End.Year, 
  #                     x = abs(Latitude.Centroid.Mean.x), colour = Family),
  #                 pch = 21, alpha = .5, shape = NULL) +
    
    geom_errorbar(data = temporal.diet, 
                    aes(ymin = Start.Year, ymax = End.Year, 
                        x = abs(Latitude.Centroid.Mean.x), colour = Family), alpha = .5, width = 0) +
    geom_point(data = temporal.diet, 
               aes(y = Mid.Year, x = abs(Latitude.Centroid.Mean.x),fill = Family),
               pch = 21, alpha = .5, colour = "black", size = 2) +
  #scale_fill_manual(values = c("#CCCCCC","#666666")) +
  geom_vline(xintercept = median.lat) + # Median latitude
  geom_vline(xintercept = 23.5, linetype = 2) + # Tropic lines
  annotate("text", x = 22, y = 1936.5, label = "Tropic lines", colour = "black", size = 3) +
  geom_vline(xintercept = 0, linetype = 2) + # Equator
  annotate("text", x = 2.5, y = 1935, label = "Equator", colour = "black", size = 3) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  #geom_vline(xintercept = 12.04167, linetype = 2) +
  #geom_vline(xintercept = 49.05833, linetype = 2) +
  geom_hline(yintercept = median.year, colour = "black") + # Median year
  #geom_hline(yintercept = 1994, linetype = 2, colour = "black") +
  #geom_hline(yintercept = 2006, linetype = 2, colour = "black") +
  #scale_size_continuous(values = log10(Scat.Stomach.Tissue)/2) +
  #geom_line() +
  theme_bw() + theme(legend.position = "top",
                     legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank())

# ggsave("./Results/Figures/Spatio-temporal information_UPDATED.tiff", p1,
#        width = 10, height = 8, dpi = 100)

ggsave("./Results/Figures/Final figures/Figure 2/Test.pdf", p1, 
       width = 6, height = 5, units = "in")

# Desciptive data
temporal.diet$Abs.Lat <- abs(as.numeric(as.character(temporal.diet$Latitude.Centroid.Mean.x)))

median(temporal.diet[!is.na(temporal.diet$Abs.Lat),]$Abs.Lat)
quantile(temporal.diet[!is.na(temporal.diet$Abs.Lat),]$Abs.Lat)[2]
quantile(temporal.diet[!is.na(temporal.diet$Abs.Lat),]$Abs.Lat)[4]

median(temporal.diet[!is.na(temporal.diet$Mid.Year),]$Mid.Year)
quantile(temporal.diet[!is.na(temporal.diet$Mid.Year),]$Mid.Year)[2]
quantile(temporal.diet[!is.na(temporal.diet$Mid.Year),]$Mid.Year)[4]

#temporal.diet <- as.data.frame(unique(temporal.diet[c(1:4)]))
temporal.diet <- temporal.diet[!is.na(temporal.diet$Start.Year),] 
temporal.diet$YearstoAdd = temporal.diet$End.Year - temporal.diet$Start.Year
temporal.diet$CarniBinom <- factor(temporal.diet$CarniBinom)

# Make world maps of the distribution of studies
biomes <- st_read("./Data/GIS/ecoregions_biomes/biomes/global_biomes.shp")
str(biomes)

# Transform biomes into a 
# mollweile.proj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
# biomes2 <- st_transform(biomes, mollweile.proj)

str(temporal.diet)
temporal.diet$Longitude.Centroid.Mean.x <- as.numeric(as.character(temporal.diet$Longitude.Centroid.Mean.x))
hist(temporal.diet$Latitude.Centroid.Mean.x)
hist(temporal.diet$Longitude.Centroid.Mean.x)


p2 <- ggplot() +
  geom_sf(data = biomes, fill = "lightgrey", alpha = 0.4) +
  geom_point(data = temporal.diet, aes(x = x,
                                       y = y, fill = Family), 
             alpha = .5, pch = 21, colour = "black", size = 2) +
  labs(x = "", y = "") +
  geom_hline(yintercept = 0, linetype = 1, alpha = .5) +
  geom_hline(yintercept = 23.5, linetype = 2, alpha = .5) +
  geom_hline(yintercept = -23.5, linetype = 2, alpha = .5) +
  scale_fill_manual(values = colours) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank())

ggsave("./Results/Figures/Final figures/Figure 2/Testfig2.pdf", p2, 
       width = 15, height = 8, units = "in")

# mylegend<-g_legend(p1)
# 
# grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
#                                p2 + theme(legend.position="none")),
#                    mylegend, nrow=3,heights=c(1, 1))

p3 <- ggarrange(p1,
                p2,
                ncol=1, nrow=2, labels = "auto",
                common.legend = TRUE, legend="right")
   
ggsave("./Results/Figures/Final figures/Figure 2/Combined spatio-temporal_molweile.pdf", p3, 
       width = 8, height = 8, units = "in")

# ---- Data summary: More information on temporal information ----
# What is the average length of time a study is reported for
x <- temporal.diet[!is.na(temporal.diet$End.Year),]
y <- x[!is.na(x$Start.Year),]

median(y$End.Year - y$Start.Year)
quantile(y$End.Year - y$Start.Year)

table(temporal.diet$Season)

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- || DIET RESOLUTION || ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Summarising data ----
# For each species get:
       
#       - Percentage of records at certain taxonomic resolutions
#       - Number of records using different sources
#       - Number of records using different quantification methods

# Number of records per carnivore and for each of the different methods
diet.records_perSource <- carnidiet %>%
                          dplyr::group_by(CarniBinom, Sample.Source) %>%
                          dplyr::summarise(records = length(Sample.Source))

# Get percentage of records for each sample source
perc <- diet.records_perSource %>% dplyr::group_by(Sample.Source) %>% dplyr::summarise(num = sum(records))
perc
perc$perc <- (perc$num/sum(perc$num))*100
perc

diet.records_perSource <- spread(diet.records_perSource, Sample.Source, records)

# Methods
diet.records_perMethod <- carnidiet %>%
                          dplyr::group_by(CarniBinom, Method) %>%
                          dplyr::summarise(records = length(Method))

# Get percentage of records for each sample source
perc <- diet.records_perMethod %>% dplyr::group_by(Method) %>% dplyr::summarise(num = sum(records))
perc
perc$perc <- (perc$num/sum(perc$num))*100
perc

diet.records_perMethod <- spread(diet.records_perMethod, Method, records)

# Number of unique prey species identified in species diets
carnidiet$PreyBinom <- as.factor(carnidiet$PreyBinom)
levels(carnidiet$PreyBinom)

colnames(carnidiet)
prey <- carnidiet[c(4,12)]
prey <- unique(prey)
prey <- prey[!prey$PreyBinom == "",]
prey$value <- 1

prey.consumed <- prey %>%
                 dplyr::group_by(CarniBinom) %>%
                 dplyr::summarise(prey.SR = length(CarniBinom))

# Percentage of records which are at different taxonomic resolutions
colnames(carnidiet)
prey <- carnidiet[c(4,7:12,17)]
prey$value <- 1

prey <- prey %>%
        group_by(CarniBinom, TaxonomicPrecision) %>%
        dplyr::summarise(value = length(TaxonomicPrecision))

# Get percentage of records for each taxonomic precision
perc <- prey %>% dplyr::group_by(TaxonomicPrecision) %>% dplyr::summarise(num = sum(value))
perc
perc$perc <- (perc$num/sum(perc$num))*100


prey <- prey[!prey$TaxonomicPrecision == "",]
prey <- prey[!is.na(prey$TaxonomicPrecision),]
prey <- as.data.frame(prey[c(1:3)])
prey

prey <- spread(prey, TaxonomicPrecision, value)
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
                      as.data.frame(diet.records_perMethod[c(2:11)]))
diet.summary <- cbind(prey.percentage.2[c(8,1:7)], diet.summary[c(2:32)])

# Merge the number of prey species
diet.summary <- merge(diet.summary, prey.consumed, by.y = "CarniBinom", all = TRUE)

diet.summary <- cbind(diet.summary[1:104,], total.obvs)
diet.summary[is.na(diet.summary)] <- 0
# Save this table
#write.csv(diet.summary, "./Results/Tables/Diet Summary.csv")

# Get some descriptive statistics for these
colSums(diet.summary[9:41])
colnames(diet.summary)

# Source of data collection
data.sources <- diet.summary[c(9:28)]
data.sources <- melt(data.sources)
data.sources <- data.sources %>% dplyr::group_by(variable) %>% dplyr::summarise(sum = sum(value))
data.sources <- data.sources[data.sources$sum > 100,]

data.sources$variable <- factor(data.sources$variable)
data.sources$variable <- reorder(data.sources$variable, -data.sources$sum)

# ---- Figure: Sum of records for each source of diet information ----
ggplot(data.sources, aes(x = variable, y = sum)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))

# Methods of data quantification
data.sources <- diet.summary[c(30:39)]
data.sources <- melt(data.sources)
data.sources <- data.sources %>% dplyr::group_by(variable) %>% dplyr::summarise(sum = sum(value))
data.sources <- data.sources[data.sources$sum > 100,]

data.sources$variable <- factor(data.sources$variable)
data.sources$variable <- reorder(data.sources$variable, -data.sources$sum)

# ---- Figure: Sum of records for each diet quantification method used ----
ggplot(data.sources, aes(x = variable, y = sum)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))

# ---- Figure: Number of records used from each source and methods used in each source ----
diet.method.summary <- carnidiet %>% dplyr::group_by(Sample.Source, Method) %>% dplyr::summarise(count = length(Method))
diet.method.summary <- diet.method.summary[diet.method.summary$count > 100,]

diet.method.summary$Sample.Source <- factor(diet.method.summary$Sample.Source)
diet.method.summary$Sample.Source <- reorder(diet.method.summary$Sample.Source, -diet.method.summary$count)
diet.method.summary$Sample.Source <- factor(diet.method.summary$Sample.Source,
                                            levels(diet.method.summary$Sample.Source)[c(2,1,4,5,3,6:7)])

p1 <- ggplot(diet.method.summary, aes(x = Sample.Source, y = count, fill = Method)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sample source") + ylab("Records") +
  scale_fill_manual(values = colours) +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 4),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.position=c(.6,.85),
        legend.key.size =  unit(0.1, "in")) # Change key size in the legend)
p1

ggsave("./Results/Figures/Final figures/Figure 1/DietMethodSummary.pdf", p1,
       width = 2.5, height = 3.5, units = "in")

# ---- Figure: Number of records used for each source and the breakdown of taxonomic resolution ----
diet.resolution  <- carnidiet %>% dplyr::group_by(Sample.Source, TaxonomicPrecision) %>% dplyr::summarise(count = length(TaxonomicPrecision))
diet.resolution <- diet.resolution[diet.resolution$count > 100 |
                                     diet.resolution$Sample.Source == "Scat/Stomach",]
diet.resolution <- diet.resolution[!is.na(diet.resolution$Sample.Source),]
diet.resolution[is.na(diet.resolution$TaxonomicPrecision),]$TaxonomicPrecision <- "Other"


levels(factor(diet.resolution$Sample.Source))
diet.resolution$Sample.Source <- factor(diet.resolution$Sample.Source,
                                        levels(factor(diet.resolution$Sample.Source))[c(4,2,6,1,3,5,7)])

levels(factor(diet.resolution$TaxonomicPrecision))
diet.resolution$TaxonomicPrecision <- factor(diet.resolution$TaxonomicPrecision,
                                             levels(factor(diet.resolution$TaxonomicPrecision))[c(4,1,5,2,3,7,6)])

## ---- Figure: Final summary of record/diet information ----
p <- ggplot(diet.resolution, aes(x = Sample.Source, y = count, fill = TaxonomicPrecision)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sample source") + ylab("Records") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 4),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.position=c(.78,.76),
        legend.key.size =  unit(0.1, "in"))
p

p3 <- ggarrange(p, p1, ncol=1, nrow=2, labels = "auto")
p3

ggsave("./Results/Figures/Final figures/Figure 1/DietResolutionSummary.pdf", p3,
       width = 3, height = 6, units = "in")

p4 <- ggarrange(p + theme(axis.ticks.x = element_blank()),
                p1 + theme(axis.text.x = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.ticks.x = element_blank()),
                      ncol=2, nrow=1, labels = "auto")

ggsave("./Results/Figures/Final figures/Figure 1/DietResolutionSummary_long.pdf", p4,
       width = 4.5, height = 2.5, units = "in")


# Alte
p4 <- ggarrange(p,
                p1,
                ncol=2, nrow=1, labels = "auto")
p4

ggsave("./Results/Figures/Final figures/Figure 1/DietResolutionSummary_long.pdf", p4,
       width = 6, height = 3, units = "in")

table(carnidiet$Sample.Source)
table(carnidiet$Method)
table(carnidiet$TaxonomicPrecision)

# ---- Simulation of prey species richness with random samples ----
# First, descriptive statistics of the number of prey species recorded in diet analyses
x <- carnidiet[carnidiet$TaxonomicPrecision == "Species",]
colnames(x)
prey.info <- x[c(7:13)]

# Proportion of records 
table(prey.info$PreyOrder)
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