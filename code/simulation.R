# Simulation code for comparing star and grid transect abundance estimators
# ------------------------------------------------------------------------
#
# Author: william gaeuman
# Contact: william.gaeuman@alaska.gov
# Last update: 2017-13-06
# 
# Generates fish distribution in cylindrical coordinates (r, theta, d), where
# d is depth (+) in water column and allows repeated transect sampling of 
# distribution using both Star and grid transect estimators. Both are 
# Horvitz-Thompson estimators of abundance equal to sum of reciprocals of 
# inclusion probabilities of detected fish. We expect star transect estimator 
# to perform well for radially symmetric distributions and much better than 
# grid transect estimator when distribution is highly concentrated at center,
# like the default example below, which assumes a 600 m x 600 m study area 
# containing 1,500 fish. 
#
# Yields station width W, tibble 'fish' and vectors 'star.hat' and 'grid.hat' of 
# simulated estimates.
# ----

library(tidyverse)
set.seed(20171106)


# Create distribution of N fish on square region of width W (meters) centered at (0, 0)
N <- 1500
W <- 600
r <- rbeta(N, 1, 3) * W/2
theta <- rbeta(N, 1, 1) * 2 * pi
d <- rbeta(N, 2, 2) * 40 + 10

# x, y from r, theta: need x for grid estimator and x and y for easy plotting
x <- r * cos(theta)
y <- r * sin(theta)

# Display fish distribution
plot(x, y, xlim = c(-W/2, W/2), ylim = c(-W/2, W/2), pch = 16, cex = 0.5)
abline(h=0, lty=2)
abline(v=0, lty=2)

fish <- tibble(r, theta, x, y, d)

# Calcuate Horvitz-Thompson weights as inverse of inclusion probabilities
fish <- suppressWarnings(mutate(fish, w = 0.108316 * d,
				        grid.wt = W / w, 			
				        sine = (w / 2) / r, 				
				        delta = ifelse(sine < 1, asin(sine), pi / 2),
				        star.wt = pi / (2 * delta)))  

# Functions for determining abundance estimates given grid transect location x0 and star transect angle theta0
get.grid.hat <- function(x0, W){
	smpl <- abs(fish$x - x0) < fish$w / 2
	return(sum(fish$grid.wt[smpl]))
	}

get.star.hat <- function(theta0){
	a1 <- abs(theta0 - fish$theta)
	a2 <- abs((theta0 + pi) %% (2 * pi) - fish$theta)
	smpl <- a1 < fish$delta | a1 > (2 * pi - fish$delta) | a2 < fish$delta | a2 > (2 * pi - fish$delta)
	return(sum(fish$star.wt[smpl]))
	}

# Replicate star and grid transect sampling for uniform random x0 and theta0
R <- 10000
grid.hat <- replicate(R, get.grid.hat(runif(1, -1, 1) * W / 2, W), simplify = TRUE)
star.hat <- replicate(R, get.star.hat(runif(1) * 2 * pi), simplify = TRUE) 

# Clean up
rm(N, r, theta, d, x, y, R, get.grid.hat, get.star.hat)