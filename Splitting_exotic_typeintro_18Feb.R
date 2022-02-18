alien <- read.csv("Alien Ants_Global.csv")
library(sp)
library(raster)
library(tidyverse) 
library(lubridate)
library(knitr)
library(cowplot)
library(car)
library(lme4)

ant_polygons <- shapefile('GABI_Data_Release1.0_18012020/Bentity2_shapefile_fullres/Bentity2_shapefile_fullres.shp') 
crs(ant_polygons)

ant_polygons$area_sqkm <- area(ant_polygons) / 1000000 #Calculates area of each polygon in square kilometers
polygons_df <- as.data.frame(ant_polygons)[, -2]

antmaps_alien_merged <- merge(alien, polygons_df, by.x = "bentity2_name", by.y = "BENTITY2_N")
#write.csv(antmaps_alien_merged, "antmaps_alien_data_with_area.csv")
antmaps_alien_merged <- read.csv("antmaps_alien_data_with_area.csv")

sum.by.ant.species.intro <- antmaps_alien_merged %>% group_by (Phy, Exotic_status) %>% summarize(area_sqkm = area_sqkm)
sum.by.intro <- sum.by.ant.species.intro %>% group_by(Phy, Exotic_status) %>% summarize(total.area=sum(area_sqkm, na.rm=TRUE))

sum.by.intro <- sum.by.intro[order(sum.by.intro$Exotic_status),]
sum.by.intro.exotic <- sum.by.intro[sum.by.intro$Exotic_status == "Exotic",]
write.csv(sum.by.intro.exotic, "alien_exotic.csv")
sum.by.intro.indoor <- sum.by.intro[sum.by.intro$Exotic_status == "Indoor Introduced",]
write.csv(sum.by.intro.indoor, "alien_indoor.csv")
sum.by.intro.intercepted <- sum.by.intro[sum.by.intro$Exotic_status == "Intercepted",]
write.csv(sum.by.intro.intercepted, "alien_intercepted.csv")
