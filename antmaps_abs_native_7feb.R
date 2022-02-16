antmaps_df <- read_csv("GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020.csv")
#antmaps_df <- data.frame(antmaps_df)[, -c(1, 2, ]
library(rgdal)
library(sp)
library(raster)

ant_polygons <- st_read('GABI_Data_Release1.0_18012020/Bentity2_shapefile_fullres/Bentity2_shapefile_fullres.shp') 
#As a dataframe under the sf package
ant1_polygons <- shapefile('GABI_Data_Release1.0_18012020/Bentity2_shapefile_fullres/Bentity2_shapefile_fullres.shp') 
#As a SpatialPolygonsDataFrame under the raster package
crs(ant1_polygons) #Checking projection

poly1_df <- as.data.frame(as(as(ant1_polygons, "SpatialLinesDataFrame"),"SpatialPointsDataFrame")) 
#Converting to a dataframe With all names and coordinates


poly1_df <- as.data.frame(poly1_df)[, -c(2,3,4,5)] #The result is a dataframe 
#with the name of each polygon (e.g. "India)
#and the list of lat and long coordinates

names <- as.data.frame(ant_polygons)[, -c(2,3)] 
#A vector for all namesnof polygons

#Finding the maximum lat and long for each polygon
maxlat <- vector()
minlat <- vector()

for (i in names) {
  temp <- poly1_df[poly1_df$BENTITY2_N == i, ]
  maxlat <- append(maxlat, max(temp$coords.x2))
  minlat <- append(minlat, min(temp$coords.x2))
                 print(i)}

poly_minmax <- data.frame(names, maxlat, minlat)

#Now merging antmaps data with polygon name + minmax latitudes
antmaps_native_merged <- merge(antmaps_df, poly_minmax, by.x = "bentity2_name", by.y = "names", all.x = "TRUE")
write.csv(antmaps_native_merged, "antmaps_native_lat_2Feb22.csv")

ant_latlong <- antmaps_native_merged %>% group_by (valid_species_name, bentity2_name) %>% 
  summarize(maxlat = maxlat, minlat = minlat)
ant_latlong <- ant_latlong[!duplicated(ant_latlong), ]

#Finding overall max and min latitudes
ant_names <- ant_latlong$valid_species_name[!duplicated(ant_latlong$valid_species_name)] 
#List of all ant species in the dataset
maxlatall <- vector()
minlatall <- vector()

for (i in ant_names) {
  temp <- ant_latlong[ant_latlong$valid_species_name == i, ]
  maxlatall <- append(maxlatall, max(temp$maxlat))
  minlatall <- append(minlatall, min(temp$minlat))
  print(i)}

ant_minmaxall <- data.frame(ant_names, maxlatall, minlatall)
ant_minmaxall$midpoint <- (ant_minmaxall$maxlatall + ant_minmaxall$minlatall)/2
abs_lat_native <- abs(ant_minmaxall$midpoint)

ant_minmaxall <- data.frame(ant_minmaxall, abs_lat_native)
write.csv(ant_minmaxall, "absolute_native_lat_ants7Feb.csv")












