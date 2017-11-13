# ======================================
# Load Packages & Source Files
# ======================================

library(rootSolve)
library(tidyverse)
library(deSolve)
if (! 'R6' %in% installed.packages()[,'Package']) install.packages('R6')
library(R6)


# The sourced file specifies a `web` class that can implement either of two different 
# versions of the model (A and B).
# Upon creating a `web` object, it solves for selected unknown parameters, given 
# the 'known' parameters and equilibria.
# It also includes an internal class function `ode_solve` that solves the ODEs. 

source('aaa-class.R')



# Both models A and B are very similar to the translations of Tony's discrete time 
# formulation into ODEs that we worked with last week. The key difference is that 
# nutrients can be lost from the system directly from each pool, rather than only being 
# lost from the nutrients. This made the model much more stable, and is more realistic. 
# It also matches the way in which nutrient-dynamic models have been formulated by others.  

# Version A retains the original structure of direct density dependence for P, V, H, 
# and R. Version B only includes direct density dependence for P, while V, H, and R 
# have saturating functional responses. As you'll see, A is very stable, while B is 
# rather unstable (although better than before). That said, I think B makes more 
# mechanistic sense (and is more conventional for consumer-resource type models). 
# Something like A is probably more suitable for our purposes, although there may be 
# some modifications we can make to B to balance plausibility and stability

# Version A is very nicely behaved. For a variety of pulse values, most pools increase 
# slightly in response to midge inputs, although R increases a lot. Everything then 
# returns to equilibrium. One interesting question is whether changing the 'internal' 
# parameters of the system (e.g. attack rates or loss rates) can shift which trophic 
# levels respond most strongly.

# Version B is quite unstable, although better behaved than the original model. With 
# the midge pulse set to 500, R explodes and drives H, V, and eventually itself extinct. 
# If you set pulse to 1, you get oscillations that appear to expand in amplitude 
# indefinitely. If you set pulse to 100, then you get something in between and rather 
# unpredictable.





# ======================================
# Model A
# ======================================

#Initialize model (A is default)
foodweb_A = web$new()

# Viewing class:
foodweb_A

# Create function for time-varying midge pulse (make sure to wrap it in a list)
# pulse is the instantaneous rate of midge input over the specified time frame
# pulse_tmin and pulse_tmax define the duration over which the midge pulse occurs
foodweb_A$iM_func = list(
    function(t) {
        pulse=500; pulse_tmin=100; pulse_tmax=150
        ifelse(t > pulse_tmin & t < pulse_tmax, pulse, 0)
    }
)

#Solve ODEs
output_A = foodweb_A$ode_solve(tmin = 0, tmax = 1000, tstep = 1)

# Plot
# Absolute biomass
output_A  %>%
    gather('pool', 'biomass', -time) %>%
    group_by(pool) %>%
    # Define 'minb' to set the minimum value for the y-axis. 
    # This allows different y-scales for different facets, with the ymin set to 0 
    mutate(minb = 0) %>%
    ggplot(aes(time, biomass)) +
    facet_wrap(~pool, scales="free_y") + 
    # The horizontal lines show the initial states
    geom_hline(data = foodweb_A$initial_states, 
               aes(yintercept=biomass), color="firebrick") +
    geom_line(size = 1) +
    geom_point(aes(time, minb), shape="") +
    theme_classic()



# Relative to Equilibrium
# Note that this gives weird results when there is no deviation from equilibrium
# This is probably due to small numerical errors, but is not a big issue
output_A %>%
    gather('pool', 'biomass', -time) %>%
    filter(pool!="M") %>%
    group_by(pool) %>%
    # Scale the biomass relative to the initial state
    mutate(biomass_scale = biomass/biomass[1]) %>%
    ggplot(aes(time, biomass_scale)) +
    facet_wrap(~pool) + 
    geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1) +
    theme_classic()





# ======================================
# Model B
# ======================================

#Initialize model (specify model B)
foodweb_B = web$new(model='B')

# Viewing class:
foodweb_B

# Create function for time-varying midge pulse (make sure to wrap it in a list)
# Changed midge pulse relative to above to give less volatile dynamics
foodweb_B$iM_func = list(
    function(t) {
        pulse=0.1; pulse_tmin=100; pulse_tmax=150
        ifelse(t > pulse_tmin & t < pulse_tmax, pulse, 0)
    }
)

#Solve ODEs
output_B = foodweb_B$ode_solve(tmin = 0, tmax = 1000, tstep = 1)

# Plot
# Absolute biomass
output_B  %>%
    gather('pool', 'biomass', -time) %>%
    group_by(pool) %>%
    # Define 'minb' to set the minimum value for the y-axis. 
    # This allows different y-scales for different facets, with the ymin set to 0 
    mutate(minb = 0) %>%
    ggplot(aes(time, biomass)) +
    facet_wrap(~pool, scales="free_y") + 
    # The horizontal lines show the initial states
    geom_hline(data = foodweb_A$initial_states, 
               aes(yintercept=biomass), color="firebrick") +
    geom_line(size = 1) +
    geom_point(aes(time, minb), shape="") +
    theme_classic()



# Relative to Equilibrium
# Note that this gives weird results when there is no deviation from equilibrium
# This is probably due to small numerical errors, but is not a big issue
output_B %>%
    gather('pool', 'biomass', -time) %>%
    filter(pool!="M") %>%
    group_by(pool) %>%
    # Scale the biomass relative to the initial state
    mutate(biomass_scale = biomass/biomass[1]) %>%
    ggplot(aes(time, biomass_scale)) +
    facet_wrap(~pool) + 
    geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1) +
    theme_classic()

