# Functions for obtaining estimates of rockfish abundance and density from hyroacustic data.
# ----
# Author: William Gaeuman
# Contact: william.gaeuman@alaska.gov
# Last updated: 2018-02-01
# ----

get_haver_dist <- function(coords){
# Computes great circle distance (km) between two points given 
# lon, lat in decimal degrees using haversine formula. Coords is
# a vector of form c(lon1, lat1, lon2, lat2). Called by functions 
# display_station_data() and get_station_estimates().
# ----

coords <- pi * coords / 180
lon1 <- coords[1]; lat1 <- coords[2]; lon2 <- coords[3]; lat2 <- coords[4]
dlon <- lon2 - lon1; dlat <- lat2 - lat1
a <- (sin(dlat / 2))^2 + cos(lat1) * cos(lat2) * (sin( dlon / 2))^2
return(2 * asin(sqrt(a)) * 6362.8) # constant for Kodiak latitude

} # END FUNCTION

get_station_estimates <- function(station.name, N, window = 21, grid = TRUE){
# station.name    = station name (character)
# N               = number of transects
# window          = moving-average window size for track position coords
# grid            = Boolean indicator for grid vs star transect design
  
# Automatically loads required packages zoo, tidyverse and rgdal (and sp). Expects files 
# track.csv, fish.csv, boundary.shp, boundary.shx, boundary.prj, boundary.dbf available in 
# current working directory. Results are in terms of km.
# 
# Note: Need to be in proper station directory to run these functions.
#
# Example: 
# > get_station_estimates("NEGrid33", 11)
#
# (Function first plots data. Analyst must then delineate transects using cursor.)
#
# Output:
#    station   area obs.fish  dens se.dens abund se.abund    cv
# 1 NEGrid33 0.3355      398 27640   10970  9275     3681 0.397
# 
# Output also written to file estimates.csv in working (station) directory.
# ----

# Load necessary libraries ----
library(zoo, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(rgdal, quietly = TRUE)

# Get data, process and display ----
suppressMessages(read_csv("fish.csv", na = "-9999") -> fish)
fish <- mutate(fish, time = as.character(time)) %>% 
	arrange(time) %>% 
	fill(depth)

suppressMessages(read_csv("track.csv") -> track)
track <- mutate(track, time = as.character(time)) %>%
	arrange(time) %>%
	mutate(lon = rollmean(lon, window, fill = NA), lat = rollmean(lat, window, fill = NA)) %>% 
	drop_na(lon, lat)

boundary <- readOGR(dsn=".",layer="boundary", verbose = FALSE)
station.area <- round(boundary@polygons[[1]]@Polygons[[1]]@area / 10^6, 5) 
boundary.coords <- project(boundary@polygons[[1]]@Polygons[[1]]@coords, boundary@proj4string@projargs, inv = TRUE)

plot(boundary.coords, type = "l", xlab = "deg lon", ylab ="deg lat", main = paste(station.name," Station Area = ", station.area,"sq km"))
lines(track$lon, track$lat)
points(jitter(fish$lon), fish$lat, col = 2, pch = 16)

# Interactively obtain transect endpoints ----
endpoints <- sort(identify(track$lon, track$lat, n = 2 * N, plot = TRUE, atpen = TRUE, cex = 0.8, col = "purple"))

# Generate grouping variable assigning fish to transects (fish not on an actual transect = 0 by default) ----
trnsct.id <- integer(nrow(fish))
for(k in 1:N){ 
	trnsct.id[fish$time > track$time[endpoints[2 * k - 1]] & fish$time < track$time[endpoints[2 * k]]] <- k
}

# Calculate estimates ----
if(grid){
	
	# Grid methodology ----
  
  # Get length (m) of orthogonal projection of station area onto baseline ----
  L <- read_csv("baseline.csv")[[1]] 
  
  # Calculate grid transect Horwitz-Thompson abundance estimates ----
  estimates <- mutate(fish, trnsct = trnsct.id) %>% 
    filter(trnsct != 0) %>% # Remove fish not on an actual transect ----
    group_by(trnsct) %>% 
    summarise(ht.abund = L * sum(1 / (0.108316 * depth))) # L and depth both in meters
  
	} else {

	# Star methodology ----
	
	# Get center point and make sure coords in proper order for distance calcs ----  
  cntr <- read_csv("center.csv")
  cntr <- c(cntr$lon, cntr$lat)
	
  # Define helper function for computing distance (km) of each fish from star center ----
	d_to_cntr <- function(lon, lat) apply(cbind(lon, lat), 1, function(coords) get_haver_dist(c(cntr, coords)))
  
  # Calculate star transect Horwitz-Thompson abundance estimates (strip.wdth and d.to.cntr both in meters) ---- 
	estimates <- mutate(fish, trnsct = trnsct.id,	
			          strip.wdth = 0.108316 * depth,
				  d.to.cntr = d_to_cntr(lon, lat) * 1000,
				  sine = strip.wdth / (2 * d.to.cntr),
				  fish.wt = ifelse(sine < 1, pi / (2 * asin(sine)), 1)) %>%
    filter(trnsct != 0) %>% # Remove fish not on an actual transect ----
		group_by(trnsct) %>% 
		summarise(ht.abund = sum(fish.wt))
	}

# Include zero estimates for transects with no fish, if any ----
estimates <- left_join(tibble(trnsct = 1:N), estimates, by = "trnsct") %>%
  replace_na(list(ht.abund = 0)) %>%
  
  # Calculate overall estimates of abundance and density, standard errors and cv ----
  summarise(dens = mean(ht.abund) / station.area,
          se.dens = sqrt(var(ht.abund) / N) / station.area,
          abund = station.area * dens,
          se.abund = station.area * se.dens,
          cv = round(se.abund / abund, 3))

# Package results, write to file estimates.csv, and return ----	
cbind(station = station.name, area = station.area, num.trnscts = N, obs.fish = sum(trnsct.id !=0 ), signif(estimates, 4) ) %>% 
	write_csv("estimates.csv") %>% 
	return()

} # END FUNCTION

display_station_data <-function(station.name, grid = TRUE, window = 21){
# This is front end of function get_station_estimates() for data display, except 
# neither number of transects nor estimation method are included in input.
#----

library(zoo, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(rgdal, quietly = TRUE)

# Get data, process and display
suppressMessages(read_csv("fish.csv", na = "-9999") -> fish)
fish <- mutate(fish, time = as.character(time)) %>% 
	arrange(time) %>% 
	fill(depth)

suppressMessages(read_csv("track.csv") -> track)
track <- mutate(track, time = as.character(time)) %>%
	arrange(time) %>%
	mutate(lon = rollmean(lon, window, fill = NA), lat = rollmean(lat, window, fill = NA)) %>% 
	drop_na(lon, lat)

boundary <- readOGR(dsn=".",layer="boundary", verbose = FALSE)
station.area <- round(boundary@polygons[[1]]@Polygons[[1]]@area / 10^6, 5) 
boundary.coords <- project(boundary@polygons[[1]]@Polygons[[1]]@coords, boundary@proj4string@projargs, inv = TRUE)

plot(boundary.coords, type="l", xlab="deg lon", ylab="deg lat", main = paste(station.name," Station Area = ", station.area,"sq km"),
	sub = paste("Total track fish count = ", nrow(fish),"; Average fish depth = ", round(mean(fish$depth), 1),"m"))
lines(track$lon, track$lat)
points(jitter(fish$lon), fish$lat, col=2, pch=16)

# Plot star center for star transect sampling ----
if(!grid) {
	cntr <- read_csv("center.csv")
	points(cntr$lon, cntr$lat, pch = 16, cex = 1.5)
}

} # END FUNCTION