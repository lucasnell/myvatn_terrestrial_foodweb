# ======================================
# Load Packages & Source Files
# ======================================

library(rootSolve)
library(tidyverse)
library(deSolve)



# source('Define_Model_A.R')
# source('Define_Model_B.R')
source('aaa-class.R')

# (Model A is the default)
foodweb_A = web$new()
# Viewing class:
foodweb_A

# How you change midge function (make sure to wrap it in a list):
foodweb_A$iM_func = list(
    function(t) {
        pulse=500; pulse_tmin=100; pulse_tmax=150
        ifelse(t > pulse_tmin & t < pulse_tmax, pulse, 0)
    }
)
output_A = foodweb_A$ode_solve(tmin = 0, tmax = 1000, tstep = 1)


foodweb_B = foodweb_A$clone()  # use clone to keep from modifying foodweb_A
foodweb_B$model = 'B'
foodweb_B$iN = 10  # this parameter was set differently from model A
foodweb_B$re_solve()  # re-solve for unknown parameters
output_B = foodweb_B$ode_solve(tmin = 0, tmax = 1000, tstep = 1)


# options("device" = "quartz"); graphics.off()

# The source files specify functions to implement two different versions of the model 
# (A and B). The objects defined in the source files have the same names, so only one 
# can be loaded at a time.

# The source files include code that solves for selected unknown paramters, given 
# the 'known' paramters and equilibria. It also includes the code that Lucas and Eric 
# developed for running the ODEs, with a wrapper to make it a bit easier to run with 
# different midge scenarios. 

# Both models A and B are very similar to the translations of Tony's discrete time 
# formulation into ODEs that we worked with last week. The key difference is that 
# nutrients can be lost from the system directly from each pool, rather than only being 
# lost from the nutrients. This made the model much more stable, and is more realistic. 
# It also matches the way in which nutrient-dynamic models have been formulated by others.  

# Version A retain's the original structure of direct density dependence for P, V, H, 
# and R. Version B only includes direct density dependence for P, while V, H, and R 
# have saturating functional responses. As you'll see, A is very stable, while B is 
# rather unstable (although better than before). That said, I think B makes more 
# mechanistic sense (and is more conventional for consumer-resource type models). 
# Something like A is probably more suitable for our purposes, although there may be 
# some modifications we can make to B to balance plausibility and stability





# ======================================
# Define parameters
# ======================================
# There are two parameter sets, "params_solve" and "equil_states" that are sourced from 
# the model files. To facilitate modfications later, but without overwriting the 
# original parameters, define new corresponding parameter sets:
params_final = params_solve
initial_states = as.numeric(equil_states)
names(initial_states) = names(equil_states)

# Add midge parameters and initial states
params_final$mM = 0.5
params_final$lM = 0.1
initial_states['M'] = 0


# ======================================
# Run Model
# ======================================
# Solve ODEs
# tmin and tmax specify the duration over which to run the model
# tstep specifies the step size
# pulse is the instantaneous rate of midge input over the specified time frame
# pulse_tmin and pulse_tmax define the duration over which the midge pulse occurs
# parms and init give the paramter values and initial states
output = 
  ode_solve(tmin=0, tmax=1000, tstep=1, pulse=500, pulse_tmin=100, pulse_tmax=150, 
            parms=params_final,init=initial_states) %>%
  as.data.frame %>%
  as_tibble

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

# Plot
# Absolute biomass
output  %>%
  gather('pool', 'biomass', -time) %>%
  group_by(pool) %>%
  # Define 'minb' to set the minimum value for the y-axis. 
  # This allows different y-scales for different facets, with the ymin set to 0 
  mutate(minb = 0) %>%
  ggplot(aes(time, biomass)) +
  facet_wrap(~pool, scales="free_y") + 
  # The horizontal lines show the initial states
  geom_hline(data = initial_states %>% t %>% tbl_df %>% 
               gather('pool','biomass'),
             aes(yintercept=biomass), color="firebrick") +
  geom_line(size = 1) +
  geom_point(aes(time, minb), shape="") +
  theme_classic()



# Relative to Equilibrium
# Note that this gives funny results when there is supposed to be no change, due to 
# tiny fluctuations about equilibrium, presumably due to numerical errors. I bet this 
# would go away if the step size was diminished, although I don't think its all that 
# important.
output %>%
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

# Comparing results
output_A %>%
    gather('pool', 'biomass', -time) %>%
    filter(pool!="M") %>%
    group_by(pool) %>%
    # Scale the biomass relative to the initial state
    mutate(biomass_scale = biomass/biomass[1]) %>%
    ggplot(aes(time, biomass_scale)) +
    facet_wrap(~pool) + 
    geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1, color = 'dodgerblue') +
    theme_classic()
