# Functions for obtaining estimates of rockfish abundance and density from hyroacustic data.
# -----------------------------------------------------------------------------------------
#
# Author: William Gaeuman
# Contact: william.gaeuman@alaska.gov
# Last updated: 2017-11-03
#
# Some details have still to be finalized...
# ----

get_haver_dist <- function(coords){
# --------------------------------
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

get_station_estimates <- function(station.dir.path, station.name, N, window = 21, grid = TRUE){
# --------------------------------------------------------------------------------------------
# Automatically loads required packages zoo, tidyverse and rgdal (and sp). Expects files 
# track.csv, fish.csv, boundary.shp, boundary.shx, boundary.prj, boundary.dbf available in 
# station directory. N = number of transects. window = moving-average window size for track 
# position coords. Results are in terms of km.
#
# Example: 
# > get_station_estimates("C:/rockfish/data/NEGrid33", "NEGrid33", 11)
#
# (Function first plots data. Analyst must then delineate transects using cursor.)
#
# Output:
#    station   area obs.fish  dens se.dens abund se.abund    cv
# 1 NEGrid33 0.3355      398 27640   10970  9275     3681 0.397
# ----

library(zoo, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(rgdal, quietly = TRUE)

setwd(station.dir.path)

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

plot(boundary.coords, type = "l", xlab = "deg lon", ylab ="deg lat", main = paste(station.name," Station Area = ", station.area,"sq km"))
lines(track$lon, track$lat)
points(jitter(fish$lon), fish$lat, col = 2, pch = 16)

# Get transect endpoints and generate transect grouping variables
endpoints <- sort(identify(track$lon, track$lat, n = 2 * N, plot = TRUE, atpen = TRUE, cex = 0.8, col = "purple"))

trnsct.id <- integer(nrow(track))
fish.grp <- integer(nrow(fish))
for(k in 1:N){ 
	trnsct.id[ endpoints[2 * k - 1]:endpoints[2 * k] ] <- k
	fish.grp[ fish$time > track$time[ endpoints[2 * k - 1] ] & fish$time < track$time[ endpoints[2 * k] ] ] <- k
	}

# Calculate estimates
if(grid){
	
	# Grid methodology
	get_trnsct_lngth <- function(lon, lat){
		n <- length(lon)
		coords <- cbind(lon[-n], lat[-n], lon[-1], lat[-1])
		return(sum(apply(coords, 1, get_haver_dist)))
		}
	
	trnsct.lngths <- mutate(track, trnsct = trnsct.id) %>% 
		filter(trnsct != 0) %>% 
		slice(seq(1, nrow(track), by = 20)) %>% 
		group_by(trnsct) %>% 
		summarise(lngth = get_trnsct_lngth(lon, lat))

	wtd.fish.cnts <- mutate(fish, trnsct = fish.grp) %>% 
		filter(trnsct != 0) %>% 
		group_by(trnsct) %>% 
		summarise(cnt = 1000 * sum(1 / (0.108316 * depth)))
		
	estimates <- left_join(trnsct.lngths, wtd.fish.cnts, by = "trnsct") %>%
		mutate(cnt = replace(cnt, is.na(cnt), 0), trnsct.den = cnt / lngth * station.area) %>% 
		summarise( dens = mean(trnsct.den), 
		   	   se.dens = sqrt(var(trnsct.den) / N),  
                           abund = station.area * dens, 
                           se.abund = station.area * se.dens,
		           cv = round(se.dens / dens, 3))
	} else {

	# Star methodology
	cntr <- c(-152.150, 57.804) # SURROGATE cntr FOR STATION NEG33; TBD HOW TO COME UP WITH THIS POINT IN GENERAL

	d_to_cntr <- function(lon, lat) apply(cbind(lon, lat), 1, function(coords) get_haver_dist(c(cntr, coords)))
	
	suppressWarnings(trnsct.abund.estimates <- mutate(fish, trnsct = fish.grp,
			           	                        strip.wdth = 0.108316 * depth / 1000,
				                                d.to.cntr = d_to_cntr(lon, lat),
				                                sine = strip.wdth / (2 * d.to.cntr),
				                                fish.wt = ifelse(sine < 1, pi / (2 * asin(sine)), 1))) %>%
		filter(trnsct != 0) %>% 
		group_by(trnsct) %>% 
		summarise(abund = sum(fish.wt))
	
	estimates <- left_join(tibble(trnsct = 1:N), trnsct.abund.estimates, by = "trnsct") %>%
		mutate(abund = replace(abund, is.na(abund), 0)) %>%
		summarise( dens = mean(abund) / station.area,
			   se.dens = sqrt(var(abund) / N) / station.area,
			   abund = mean(abund),
			   se.abund = station.area * se.dens,
			   cv = round(se.dens / dens, 3))
	}

# Package results, store in a file, and return	
cbind(station = station.name, area = station.area, obs.fish = sum(fish.grp !=0 ), signif(estimates, 4) ) %>% 
	write_csv("estimates.csv") %>% 
	return()

} # END FUNCTION

display_station_data <-function(station.dir.path, station.name, window = 21){
# --------------------------------------------------------------------------
# This is front end of function get_station_estimates() for data display, except 
# neither number of transects nor estimation method are included in input.
#----

library(zoo, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(rgdal, quietly = TRUE)

setwd(station.dir.path)

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

} # END FUNCTION