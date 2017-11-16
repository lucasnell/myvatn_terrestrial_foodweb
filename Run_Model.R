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

iM_func = function(t) {
        pulse=500; pulse_tmin=100; pulse_tmax=150 
        ifelse(t > pulse_tmin & t < pulse_tmax, pulse, 0) 
    }

many_webs <- multi_web(list(Neq = seq(343000, 343500, 250), 
                            Deq = seq(114000, 114100, 50), 
                            iM_func = list(iM_func)), 
                       expand = TRUE)
many_outputs <- many_webs %>% 
    lapply(function(w) w$ode_solve(tmin = 0, tmax = 1000, tstep = 1) %>% 
               mutate(Neq = w$Neq, Deq = w$Deq)) %>% 
    bind_rows

many_outputs %>% 
    gather('pool', 'biomass', -time, -Neq, -Deq) %>%
    filter(pool != "M") %>%
    group_by(pool, Neq, Deq) %>%
    # Scale the biomass relative to the initial state
    mutate(biomass_scale = biomass/biomass[1]) %>%
    ungroup %>% 
    ggplot(aes(time, biomass_scale, color = pool)) +
    facet_grid(Deq ~ Neq, labeller = label_both) + 
    geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1) +
    theme_classic()



# ======================================
# Model A
# ======================================

#Initialize model (A is default)
foodweb_A = web$new(model="A")

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
        pulse=0.01; pulse_tmin=100; pulse_tmax=150
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

