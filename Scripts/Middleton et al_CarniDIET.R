#————————————————————————————————————————————————————————————————————————————————————————————####
#
#  ---- CarniDIET V1.0 ----
#   
#   A database on the diets of the world's terrestrial carnivorous mammals
#   
#   Script: Summarise data within CarniDIET for the manuscript "TBC"
#   
#   - 1. Metadata: Information for database metadata
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
# ---- Book keeping ----
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
library(pacman)
pacman::p_load(tidyverse,plyr,ggrepel,sp,sf,viridis,rgdal,
               rgeos,maptools,reshape2,RColorBrewer,gridExtra,
               grid,doParallel,effects,raster,cluster,factoextra,
               vegan,FD,factoextra,pscl,boot,mapproj,ggplot2, dplyr,
               stringr, gridExtra,viridisLite,raster,rgdal,maptools,
               ape, ggtree,ggpubr, update = F)


#---- Load data ----
# Load CarniDIET
carnidiet <- read.csv("./Version 1.0/CarniDIET 1.0.csv")

# Load potential species list
# Includes species selected as primary consumers of mammals from MammalDIET (Kissling et al., 2014)
cd.wos.hits <- read.csv("./Version 1.0/Supplementary data/CarniDiet_PotentialSpecies.csv", header = TRUE)

# Phylacine 
phylacine <- read.csv("../Phylacine/Trait_data.csv")

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- CarniDIET Metadata ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- 1. Carnivore species information ----

# How many ... in CarniDIET
nlevels(carnidiet$CarniBinom) # 104 species
nlevels(carnidiet$Order) # 5 orders
nlevels(carnidiet$Family)  # 15 families

# How many levels of:
nlevels(carnidiet$Age) # 4 ages/life-stages
nlevels(carnidiet$Sex) # Both sexes and both/unknown

# ---- 2. Diet resolution information -----

# How many levels for food-type variables:
# the '-1' accounts for "" entries which could be treated as NA's...
levels(carnidiet$Food.type) # 15 Food types
levels(carnidiet$PreyOrder) # 62 prey orders
levels(carnidiet$PreyFamily) # 130 prey families
levels(carnidiet$PreyGenus) # 490 prey genera
levels(carnidiet$PreySpecies) # 734 prey species
levels(carnidiet$PreyBinom) # 854 prey binomials
levels(carnidiet$Common.name) # 1649 common names

# How many domestic animals
dom <- carnidiet %>% filter(Domestic.Agricultural == 1) %>% droplevels()
levels(dom$PreyBinom)

# Categories of taxonomic precision
levels(carnidiet$TaxonomicPrecision) # 10 categories of taxonomic precision

# Range of quantitative importance of prey species or food type
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
levels(carnidiet$Season) # 56 seasons

# ----- 4. Methods information ----
levels(carnidiet$Method) # 9 quantification methods
levels(carnidiet$Sample.Origin) # 20 Diet Sample Origins

# ----- 5. Spatial information -----
# Altitude
range(carnidiet[!is.na(carnidiet$AltitudeMinimum),]$AltitudeMinimum) # 0 - 4543m
range(carnidiet[!is.na(carnidiet$AltitudeMaximum),]$AltitudeMaximum) # 0 - 8156m

# Latitude and longitude
range(carnidiet[!is.na(carnidiet$Latitude.Centroid.Mean.x),]$Latitude.Centroid.Mean.x) # -55S to 81N
range(carnidiet[!is.na(carnidiet$Longitude.Centroid.Mean.x),]$Longitude.Centroid.Mean.x) # -169W to 178E
levels(carnidiet$Coordinates.Source) # 4 Types of geo-referencing

# Descriptions of study area
levels(carnidiet$SiteName) # 791 sites
levels(carnidiet$Region) # 410 Regions
levels(carnidiet$Country) # 107 "countries"

# Range of study area size
carnidiet$Study.Area.Size.km2 <- as.numeric(as.character(carnidiet$Study.Area.Size.km2))
range(carnidiet[!is.na(carnidiet$Study.Area.Size.km2),]$Study.Area.Size.km2) # 0.03 to 100,000km2


# ---- 6. Bibliographic information ----
levels(carnidiet$PR.Author) # 638
range(carnidiet$PR.Year) # 1952 to 2019
levels(carnidiet$PR.Title) # 749 titles of studies
levels(carnidiet$PR.Journal) # 208 Journals
levels(carnidiet$PrimaryRef) # 749 sources
levels(carnidiet$CollectionRef) # 683 Collection Refs

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- Taxonomy represented  ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# How many potential ...?
length(unique(cd.wos.hits$Bin.)) # 210 species
levels(cd.wos.hits$Order) # 9 orders
levels(cd.wos.hits$Family) # 23 families

# ---- Family representation in diet sources -----
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

# Number of sources included in CarniDIET (by Family)
family.papers <- carnidiet[c(3,44)]
family.papers <- unique(family.papers)
family.papers <- family.papers %>% dplyr::group_by(Family) %>% dplyr::summarise(No.papers = length(Family))
family.papers

# Median and max of number of sources per species (per Family)
family.species.papers <- carnidiet[c(3:4,44)]
family.species.papers <- unique(family.species.papers)
family.species.papers <- family.species.papers %>% dplyr::group_by(Family,CarniBinom) %>% dplyr::summarise(studiesperspecies = length(CarniBinom))
(family.species.papers <- family.species.papers %>% dplyr::group_by(Family) %>% dplyr::summarise(median.papers = median(studiesperspecies),
                                                                                                max.papers = max(studiesperspecies)))

# Building Family-level table
family.info <- merge(family.species,family.CD, by = "Family", all = TRUE)
family.info$Perc.Families <- (family.info$No.Species.y/family.info$No.Species.x)*100
family.info <- merge(family.info, family.wos.hits, by = "Family", all = TRUE)
family.info <- merge(family.info, family.papers, by = "Family", all = TRUE)
family.info$Hits_Papers_Ratio <- (family.info$No.papers/family.info$No.Hits)*100
family.info <- merge(family.info, family.species.papers, by = "Family", all = TRUE)
family.info[is.na(family.info)] <- 0

# How many sources were used compared to actual hits?
sum(family.info$No.Hits)
sum(family.info$No.papers)
median(family.info$Hits_Papers_Ratio)

# Merge in Order information
colnames(phylacine)
order.family <- unique(phylacine[c(2,3)])
colnames(order.family) <- c("Order", "Family")

family.info <- merge(family.info,order.family, by = "Family", all.x = TRUE)
family.info <- family.info[c(10,1:9)]

family.info <- family.info %>% arrange(Order, Family)

# ---- Genus representation in diet sources -----
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

# Number of actual sources per genus
genus.papers <- diet.2[c(3:4,45)]
genus.papers <- unique(genus.papers)
genus.papers <- genus.papers %>% dplyr::group_by(Family,Genus_V_1.2) %>% dplyr::summarise(No.papers = length(Genus_V_1.2))

# Median and max of number of sources per species per Genus
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
colnames(genus.info)[2] <- "Family"
genus.info <- merge(genus.info, order.family, by = "Family", all.x = TRUE)
genus.info <- genus.info[c(11,1:10)]

genus.info <- genus.info %>% arrange(Order, Family, Genus_V_1.2)

# ---- Tables: Family & genus level representation in dietary studies ----
write.csv(family.info, "../Tables/Family.Information.Summary_UPDATE.csv")
write.csv(genus.info, "../Tables/Genus.Information.Summary_UPDATE.csv")

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
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9)))

# Save figure
ggsave(filename = paste0("../Final figures/Figure S1/PercentFamiliesStudiedUpdated.tiff"), plot = plot, device = NULL, path = NULL,
       scale = 1, width = 10.5, height = 10.5, units = c("cm"),
       dpi = 300)



# Supplementary figure: Barplot for top-20 studied species in database ----

# There are 749 sources used. How many feature specific species?
# There are 3160 unique studies. How many are for specific species?

# Number of sources feature each species
carnivore.summary <- carnidiet %>%
                     dplyr::group_by(Family, CarniBinom) %>%
                     dplyr::summarise(papers = length(unique(PrimaryRef)))
 
# Number of studies feature each species
carnivore.summary.2 <- carnidiet %>% 
                       dplyr::group_by(Family, CarniBinom, Age, Sex, # Unique species and demography
                                       Scat.Stomach.Tissue,Prey.Items.Kill.Sample.Size, # Unique methods
                                       Sample.Origin, Method, 
                                       Country, Region, SiteName, AltitudeMaximum, # Unique spatial information
                                       Latitude.Centroid.Mean.x, Longitude.Centroid.Mean.x,
                                       Start.Year,End.Year, Start.Month,End.Month, Season, # Unique temporal information
                                       PR.Title) %>%  # Unique source
                       dplyr::summarise(x = length(CarniBinom))

carnivore.summary.3 <- carnivore.summary.2 %>%
  dplyr::group_by(Family, CarniBinom) %>%
  dplyr::summarise(studies = length(CarniBinom)) # Number of studies

carnivore.summary.master <- merge(carnivore.summary, carnivore.summary.3, by = "CarniBinom")

# How many studies are there
sum(carnivore.summary.master$papers) # 950, unique combinations of sources and carnivore
length(unique(carnidiet$PrimaryRef)) # 749, unique sources

# How many studies are there
sum(carnivore.summary.master$studies) # 3160 - This is all individual studies
3160/749 # On average, there are 4 studies per source

# 475 source = 50% source
# 562 = 75% source

# Tidy table
colnames(carnivore.summary.master)[1] <- "CarniBinom"

# Get top-20 species
carnivore.summary.master$CarniBinom <- reorder(carnivore.summary.master$CarniBinom, carnivore.summary.master$studies)
x <- carnivore.summary.master[order(carnivore.summary.master$studies, decreasing = TRUE),]
y <- x[1:20,]
y$CarniBinom <- factor(y$CarniBinom)
levels(y$CarniBinom)

# How many species need to be included to get 50% of the studies?
(sum(x[1:9,]$studies)/3160)*100 #  9 species includes just over 50% studies
(sum(x[1:20,]$studies)/3160)*100 # 20 species includes 70% studies
(sum(x[1:24,]$studies)/3160)*100 # 24 species includes almost 75% studies

# Find percentage lines
(sum(x[1:1,]$studies)/3160)*100 # 14.8% studies
(sum(x[1:9,]$studies)/3160)*100 # 50% studies

# Rename species factor levels for figure
y$CarniBinom <- revalue(y$CarniBinom, c("Vulpes_vulpes" = "V. vulpes","Lynx_lynx" = "L. lynx",
                                        "Canis_lupus" = "C. lupus","Felis_silvestris" = "F. silvestris",
                                        "Panthera_pardus" = "P. pardus","Neovison_vison" = "N. vison",
                                        "Puma_concolor" = "P. concolor","Martes_foina" = "M. foina",
                                        "Canis_latrans" = "C. latrans","Genetta_genetta" = "G. genetta",
                                        "Ursus_arctos" = "U. arctos","Canis_aureus" = "C. aureus",
                                        "Martes_martes" = "M. martes","Mustela_erminea" = "M. erminea",
                                        "Panthera_tigris" = "P. tigris","Panthera_uncia" = "P. uncia",
                                        "Martes_americana" = "M. americana","Lynx_rufus" = "L. rufus",
                                        "Canis_mesomelas" = "C. mesomelas","Panthera_onca" = "P. onca"))


# Main plot
(p10 <- ggplot(data = y[1:20,], aes(x = reorder(CarniBinom, studies), y = studies, fill = papers)) + # 20 species refers to 66% all papers
  geom_bar(stat = "identity") +
  #annotate("rect", xmin = 19.5, xmax = 20.5, ymin = -Inf, ymax =Inf, fill = "red") +
  #annotate("rect", xmin = 15.5, xmax = 19.5, ymin = -Inf, ymax =Inf, fill = "orange") +
  #annotate("rect", xmin = 7.5, xmax = 15.5, ymin = -Inf, ymax =Inf, fill = "yellow") +
  coord_flip() +
  xlab("") + ylab("Number of studies") +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Sources") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 19.5, colour = "black", linetype = 2, alpha = 0.8) + # 15% studies
  geom_vline(xintercept = 11.5, colour = "black", linetype = 2, alpha = 0.8) + # 50% studies
  annotate("text", x = 19, y = 380, label = "15%", size = 3) +
  annotate("text", x = 11, y = 380, label = "50%", size = 3) +
  ylim(0,500) + theme(axis.text.y = element_text(face = "italic")))

# Save plot
ggsave(filename = paste0("../Final figures/Figure S2/Top-20 studied species.tiff"), plot = p10, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c("cm"),
       dpi = 200)



# Figure: Phylogenetic tree and number of studies per species ----

# Download phylogeny from PHYLACINE
# Save  phylogeny in  repository behind the one with this project
forest <- read.nexus("../Phylogenies/Complete_phylogeny.nex")
names(forest) <- NULL
set.seed(42)
forest <- forest[sample(1:1000, 30)]

# Species list in CarniDIET
(sp <- levels(carnivore.summary.master$CarniBinom))

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
    geom_tippoint(aes(x, y, fill = Family.x, size = studies), pch = 21, colour = "black", alpha = 0.75) +
    scale_fill_manual(values = colours) +
    #xlim(-160, 30) +
    theme(panel.background = element_blank(),
          legend.position = "none"))

# Save phylogenetic tree
ggsave("../Final figures/Figure 2/Phylo tree.pdf", p,
       width = 4, height = 4, units = "in")



# Histogram of studies per species
(hist <- ggplot(data = carnivore.summary.master, aes(x = log10(studies))) +
    geom_histogram(binwidth = 0.2, colour = "black", fill = "lightgrey", alpha = 0.5) +
    scale_x_continuous(labels = c("0" = "1",
                                  "1" = "10",
                                  "2" = "100",
                                  "3" = "1000")) +
    labs(x = "Number of studies", y = "Number of species") +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()))

# Save plot
ggsave("../Final figures/Figure 2/Histogram.pdf", hist,
       width = 2, height = 2, units = "in")

# Alternative: Barplot showing decline in studies and coloured by family

(barplot <- ggplot(data = carnivore.summary.master, aes(x = reorder(Binomial.1.2, -studies),
                                                        fill = Family.x,
                                                     y = studies)) +
    geom_bar(stat = "identity", width = 1) +
    # scale_y_continuous(labels = c("0" = "1",
    #                                "1" = "10",
    #                                "2" = "100",
    #                                "3" = "1000")) +
    scale_fill_manual(values = colours, guide = "none") +
    labs(x = "Species", y = "Number of studies") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))

ggsave("../Final figures/Figure 2/barplot.pdf", barplot,
       width = 2, height = 2, units = "in")



# Add figures together later in Inkscape



# Supplementary figure: Repeat with all mammal-consumers ----

# Reload phylogeny
forest <- read.nexus("../Phylogenies/Complete_phylogeny.nex")
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
    geom_tippoint(aes(x, y, fill = Family, size = studies), pch = 21, colour = "black", alpha = 0.75) +
    scale_fill_manual(values = colours.updated) +
    #xlim(-160, 30) +
    theme(panel.background = element_blank(),
          legend.position = "none"))

# Save plot
ggsave("../Final figures/Figure S3/Phylo tree_ALL Species.pdf", p,
       width = 5, height = 5, units = "in")

# Get histogram of number of studies for all potential species:
tree[is.na(tree$studies),]$studies <- 0 # Turn NA to 0

# Plot
(hist <- ggplot(data = tree[!is.na(tree$label),], aes(x = log10(studies + 1))) +
    geom_histogram(binwidth = 0.1, colour = "black", fill = "lightgrey") +
    # scale_x_continuous(labels = c("0" = "1",
    #                               "1" = "10",
    #                               "2" = "100",
    #                               "3" = "1000")) +
    labs(x = "Number of studies", y = "Number of species") +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()))

# Save histogram
ggsave("../Final figures/Figure S3/Histogram_All potential.pdf", hist,
       width = 2, height = 2, units = "in")


#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- Spatio-temporal distribution of dietary studies ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# Get unique studies ----
# Find the number of studies
temporal.diet <- carnidiet %>% 
  dplyr::group_by(Family, CarniBinom, Age, Sex, # Unique species and demography
                  Scat.Stomach.Tissue,Prey.Items.Kill.Sample.Size, # Unique methods
                  Sample.Origin, Method, 
                  Country, Region, SiteName, AltitudeMaximum, # Unique spatial information
                  Latitude.Centroid.Mean.x, Longitude.Centroid.Mean.x,
                  Start.Year,End.Year, Start.Month,End.Month, Season, # Unique temporal information
                  PR.Title) %>%  # Unique source
  dplyr::summarise(x = length(CarniBinom))

nrow(temporal.diet) # 3160 studies

# Get midpoints for each study:
temporal.diet$Mid.Year <- (temporal.diet$Start.Year + temporal.diet$End.Year)/2

# Figure: Spatio-temporal distribution ----

# Find median and IQR of study mid-year
(mid.years <- temporal.diet$Mid.Year[!is.na(temporal.diet$Mid.Year)])

(median.year <- quantile(mid.years)[3])
(LQ.year <- quantile(mid.years)[2])
(UQ.year <- quantile(mid.years)[4])

# Find median and IQR of study latitudes
latitudes <- temporal.diet$Latitude.Centroid.Mean.x[!is.na(temporal.diet$Latitude.Centroid.Mean.x)]

(median.lat <- quantile(abs(latitudes))[3])
(LQ.lat <- quantile(abs(latitudes))[2])
(UQ.lat <- quantile(abs(latitudes))[4])

# Get colour scheme consistent again:
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

# Plot figure showing temporal distribution and |latitude| of studies
p1 <- ggplot() +
  annotate("rect",xmin =-Inf, xmax = Inf, ymin = LQ.year, ymax =UQ.year, alpha = .2) +
  annotate("rect",xmin = LQ.lat, xmax = UQ.lat, ymin = -Inf, ymax =Inf, alpha = .2) +
  annotate("rect",xmin = LQ.lat, xmax = UQ.lat,ymin = LQ.year, ymax = UQ.year, alpha = .2) +
  ylab("Year") + xlab("|Latitude|") +
  coord_flip() +
  scale_y_continuous(breaks = c(1940,1960,1980,2000,2020)) + 
    geom_errorbar(data = temporal.diet, 
                    aes(ymin = Start.Year, ymax = End.Year, 
                        x = abs(Latitude.Centroid.Mean.x), colour = Family), alpha = .5, width = 0) +
    geom_point(data = temporal.diet, 
               aes(y = Mid.Year, x = abs(Latitude.Centroid.Mean.x),fill = Family),
               pch = 21, alpha = .5, colour = "black", size = 1.5) +
  geom_vline(xintercept = median.lat) + # Median latitude
  geom_vline(xintercept = 23.5, linetype = 2) + # Tropic lines
  annotate("text", x = 22, y = 1936.5, label = "Tropic lines", colour = "black", size = 3) +
  geom_vline(xintercept = 0, linetype = 2) + # Equator
  annotate("text", x = 2.5, y = 1935, label = "Equator", colour = "black", size = 3) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  geom_hline(yintercept = median.year, colour = "black") + # Median year
  theme_bw() + theme(legend.position = "top",
                     legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank())

# Save figure
ggsave("../Final figures/Figure 1/Spatio-temporal information_UPDATED.tiff", p1,
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
temporal.diet.no.nas <- temporal.diet[!is.na(temporal.diet$Latitude.Centroid.Mean.x),]
temporal.diet.no.nas <- temporal.diet[!is.na(temporal.diet$Longitude.Centroid.Mean.x),]

# change to sf object and will change to Mollweide projection...difficult
coords <- as.data.frame(temporal.diet.no.nas[c(14,13)])

# Get Spatial points
temporal.diet.sp <- SpatialPointsDataFrame(coords = coords,
                                  data = temporal.diet.no.nas,
                                  proj4string = CRS(crs(biomes)))
print(temporal.diet.sp)

temporal.diet.sp <- st_as_sf(temporal.diet.sp)

# Need to get colours right again...
levels(factor(temporal.diet$Family))
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

levels(factor(temporal.diet.sp$Family))

colours.updated.again <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
                           "#FFCCCC","#CC6699","#FF99FF","#9999CC","#3399FF","#CCFFFF","#003333")

# Merge taxonomic info 
temporal.diet.sp <- merge(temporal.diet.sp, order.family, by = "Family", all.x = TRUE)
levels(temporal.diet.sp$Order)
temporal.diet.sp$Carnivora <- NA

temporal.diet.sp[temporal.diet.sp$Order == "Carnivora",]$Carnivora <- "Carnivora"
temporal.diet.sp[temporal.diet.sp$Order != "Carnivora",]$Carnivora <- "Non-Carnivora"


# Plot map
p2 <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", alpha = 0.4) +
  geom_sf(data = temporal.diet.sp, aes(fill = Family), 
          alpha = .4, pch = 21, colour = "black", size = 1.5) +
  coord_sf(crs = "+proj=moll") +
  labs(x = "", y = "") +
  scale_fill_manual(values = colours.updated.again) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank())

ggsave("../Final figures/Figure 2/Testfig2.pdf", p2, 
       width = 15, height = 8, units = "in")

p3 <- ggarrange(p1,
                p2,
                ncol=1, nrow=2, labels = "auto",
                common.legend = TRUE, legend="right")
   
ggsave("../Final figures/Figure 1/Combined spatio-temporal_mollweide.pdf", p3, 
       width = 8, height = 8, units = "in")

# ---- Data summary: More information on temporal information ----

# What is the average length of time a study is reported for
x <- temporal.diet[!is.na(temporal.diet$End.Year),]
y <- x[!is.na(x$Start.Year),]

median(y$End.Year - y$Start.Year)
quantile(y$End.Year - y$Start.Year)



# ---- How many time series are there in CarniDIET? ----

# Need to figure out which sources give different  breakdowns with
# these variables: 

# Start.Year,End.Year
# Start.Month,End.Month
# Season

# 1. Inter-annual studies (i.e. seasonal breakdowns)
colnames(temporal.diet)

time.series.diet <- temporal.diet[c(1,2,19,20)]
(time.series.diet <- unique(time.series.diet))
time.series.diet$value <- 1

(time.series.diet <- time.series.diet %>% 
  dplyr::group_by(CarniBinom, PR.Title) %>% 
  dplyr::summarise(num.season = sum(value)))

nrow(time.series.diet[time.series.diet$num.season >1,])
# 223/949 unique combinations of sources and carnivore have
# more than one season
223/949 # 23.5%

# 2. Annual time series (i.e. multiple years)
colnames(temporal.diet)

time.series.diet <- temporal.diet[c(1,2,15,20)]
(time.series.diet <- unique(time.series.diet))
time.series.diet$value <- 1

(time.series.diet <- time.series.diet %>% 
    dplyr::group_by(CarniBinom, PR.Title) %>% 
    dplyr::summarise(num.years = sum(value)))

nrow(time.series.diet[time.series.diet$num.years > 2,])
# Only 45/949 combinations of carnivore and primary referemce
44/949 # 4.6%

# 3. Annual time series that also includes inter-annual breakdowns
colnames(temporal.diet)

time.series.diet <- temporal.diet[c(1,2,15,19,20)]
(time.series.diet <- unique(time.series.diet))
time.series.diet$value <- 1

# First, how many seasons are there per year
(time.series.diet <- time.series.diet %>% 
    dplyr::group_by(CarniBinom, PR.Title, Start.Year) %>% 
    dplyr::summarise(num.season = sum(value)))

# Second, of those years that have multiple seasons, how many sources
# also have multiple years
time.series.diet <- time.series.diet[time.series.diet$num.season >1,]

time.series.diet <- time.series.diet[c(1:3)]
time.series.diet <- unique(time.series.diet)
time.series.diet$value <- 1

(time.series.diet <- time.series.diet %>% 
    dplyr::group_by(CarniBinom, PR.Title) %>% 
    dplyr::summarise(num.year = sum(value)))

nrow(time.series.diet[time.series.diet$num.year >2,])
# Only 14 unique combinations of carnivores and source have diet breakdowns
# across multiple years and multiple seasons
14/949 # 1.5%
 

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
                          dplyr::group_by(CarniBinom, Sample.Source) %>%
                          dplyr::summarise(records = length(Sample.Source))

# Get percentage of records for each sample source
perc <- diet.records_perSource %>% dplyr::group_by(Sample.Source) %>% dplyr::summarise(num = sum(records))
(perc$perc <- (perc$num/sum(perc$num))*100)

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
nlevels(carnidiet$PreyBinom)

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
prey <- carnidiet[c(4,7:12,15)]
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

# Methods of data quantification
data.sources <- diet.summary[c(30:39)]
data.sources <- melt(data.sources)
data.sources <- data.sources %>% dplyr::group_by(variable) %>% dplyr::summarise(sum = sum(value))
data.sources <- data.sources[data.sources$sum > 100,]

data.sources$variable <- factor(data.sources$variable)
data.sources$variable <- reorder(data.sources$variable, -data.sources$sum)



# ---- Figure: Number of records used from each source and methods used in each source ----
diet.method.summary <- carnidiet %>% dplyr::group_by(Sample.Source, Method) %>% dplyr::summarise(count = length(Method))
diet.method.summary <- diet.method.summary[diet.method.summary$count > 100,]

diet.method.summary$Sample.Source <- factor(diet.method.summary$Sample.Source)
diet.method.summary$Sample.Source <- reorder(diet.method.summary$Sample.Source, -diet.method.summary$count)
diet.method.summary$Sample.Source <- factor(diet.method.summary$Sample.Source,
                                            levels(diet.method.summary$Sample.Source)[c(2,1,4,5,3,6:7)])

(p1 <- ggplot(diet.method.summary, aes(x = Sample.Source, y = count, fill = Method)) +
  geom_bar(stat = "identity", colour = "black") +
    theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sample origin") + ylab("Records") +
  scale_fill_brewer(palette = "Greys") +
  theme(axis.text = element_text(size = 10, angle = 90),
        legend.text = element_text(size = 4),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.position=c(.6,.85),
        legend.key.size =  unit(0.1, "in"))) # Change key size in the legend)



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

# Plot figure
(p <- ggplot(diet.resolution, aes(x = Sample.Source, y = count, fill = TaxonomicPrecision)) +
  geom_bar(stat = "identity", colour = "black") +
    theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  scale_fill_brewer(palette = "Greys") +
  scale_y_continuous(limits = c(0,18000), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sample origin") + ylab("Records") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 4),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.position=c(.78,.76),
        legend.key.size =  unit(0.1, "in")))

# Plot both figures together
(p4 <- ggarrange(p + theme(axis.ticks.x = element_blank()),
                p1 + theme(axis.ticks.x = element_blank()),
                      ncol=2, nrow=1, labels = "auto"))

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