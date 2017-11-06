# Simulation code for comparing star and grid transect abundace estimators
# ------------------------------------------------------------------------
#
# Author: william gaeuman
# Contact: william.gaeuman@alaska.gov
# Last update: 2017-11-06
# 
# Generates fish distribution in either cylindrical coordinates (r, theta, d) or
# Cartesian coordinates (x, y, d), where d is depth (+) in water column and allows
# repeated transect sampling of distribution using both Star and grid transect
# estimators. Both are Horvitz-Thompson estimators of abundance equal to sum of
# reciprocals of inclusion probabilities of detected fish. We expect star transect 
# estimator to perform well for radially symmetric distributions and much better
# than grid transect estimator when distribution is highly concentrated at center,
# like the default example below. Set up to assume a 600 m x 600 m study area 
# containing 1,500 fish. 
#
# Yields tibble 'fish' and vectors 'star.hat' and 'grid.hat' of simulated estimates.
# ----

library(tidyverse)
set.seed(20171106)


# Create distribution of N fish
N <- 1500
r <- rbeta(N, 1, 3) * 300
theta <- rbeta(N, 1, 1) * 2 * pi
d <- rbeta(N, 2, 2) * 40 + 10

# x, y from r, theta: need x for grid estimator and x and y for easy plotting
x <- r * cos(theta)
y <- r * sin(theta)

# Display fish distribution
plot(x, y, xlim = c(-300, 300), ylim = c(-300, 300), pch = 16, cex = 0.5)
abline(h=0, lty=2)
abline(v=0, lty=2)


fish <- tibble(r, theta, x, y, d)
fish <-	suppressWarnings(mutate(fish, w = 0.108316 * d, 			
				      sine = (w / 2) / r, 				
				      delta = ifelse(sine < 1, asin(sine), pi / 2)))  

# Replicate star and grid transect sampling
R <- 5000
star.hat <- grid.hat <- numeric()
for(k in 1:R){
	
	# Draw random star transect angle
	theta.t <- runif(1) * 2 * pi
	
	# Draw random grid transect location along base line
	x0 <- runif(1, -1, 1) * 300

	fish <- mutate(fish, a1 = abs(theta.t - theta), 
		             a2 = abs((theta.t + pi) %% (2 * pi) - theta),
		             star.smpl = a1 < delta | a1 > (2 * pi - delta) | a2 < delta | a2 > (2 * pi - delta),
			     grid.smpl = abs(x - x0) < w / 2)
	
	# Calculate Horvitz_Thompson abundance estimates
	star.hat <- c(star.hat, sum(pi / (2 * fish$delta[fish$star.smpl])))
	grid.hat <- c(grid.hat, sum(600 / fish$w[fish$grid.smpl])) 
	}

rm(N, r, theta, d, x, y, R, k, theta.t, x0)
