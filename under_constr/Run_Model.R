# ======================================
# Load Package
# ======================================

library(mtf)


# The package has a `web` class that can implement either of two different
# versions of the model (A and B).
# Upon creating a `web` object, it solves for selected unknown parameters, given
# the 'known' parameters and equilibria.
# It also includes an internal class function `ode_solve` that solves the ODEs.



# ======================================
# Model A
# ======================================

# Initialize model (A is default)
foodweb_A = web$new(model = "A")

# Solve for unknown values
foodweb_A$eq_solve()

foodweb_A$R0 = 0
foodweb_A$eq_solve(solve_pars = c("N0", "D0", "P0", "V0", "H0", "M0"))

# Viewing class:
foodweb_A

# Outputting class values:
foodweb_A$values()



# Test midge pulse:
plot(foodweb_A$test_midges(1000, a=1e9, b=10, r=400, w=40, d=1), type = 'l', ylab = "iM")

# Solve ODEs
output_A = foodweb_A$ode_solve(tmax = 1000, a=1e9, b=0, r=400, w=40, d=1)  %>%
    gather('pool', 'biomass', -time)  # %>%
    # mutate(pool = factor(pool, levels = c('D', 'H', 'N', 'P', 'R', 'V', 'M'),
    #                      labels = c('detritus', 'herbivores', 'N pool', 'plants',
    #                                 'predators', 'detritivores', 'midges')))

# Plot absolute biomass
output_A  %>%
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



modA_out <- lapply(c(0, 10, 20, 40),
       function(.b) {
           foodweb_A$ode_solve(tmax = 1000, a=1e9, b=.b, r=400, w=40, d=1)  %>%
               gather('pool', 'biomass', -time) %>%
               mutate(b = .b)
       }) %>%
    bind_rows()

modA_out %>%
    mutate(b = factor(b)) %>%
    filter(pool!="M") %>%
    group_by(pool, b) %>%
    # Scale the biomass relative to the initial state
    mutate(biomass_scale = biomass/biomass[1]) %>%
    ggplot(aes(time, biomass_scale, color = b)) +
    facet_wrap(~pool, scales = 'free_y') +
    # geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1) +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 500, 1000))



# Relative to Equilibrium
# Note that this gives weird results when there is no deviation from equilibrium
# This is probably due to small numerical errors, but is not a big issue
output_A %>%
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

# Initialize model (specify model B)
foodweb_B = web$new(model='B')

# Solve for unknown values
foodweb_B$eq_solve()

# Viewing class:
foodweb_B

# Outputting class values:
foodweb_B$values()


# Create function for time-varying midge pulse (make sure to wrap it in a list)
# Changed midge pulse relative to above to give less volatile dynamics
foodweb_B$iM_func = list(
    function(t) {
        pulse=0.01; pulse_tmin=100; pulse_tmax=150
        ifelse(t > pulse_tmin & t < pulse_tmax, pulse, 0)
    }
)

# Solve ODEs
output_B = foodweb_B$ode_solve(tmax = 1000)

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







# =========
# Multiple foodwebs
# =========


# Define midge function
iM_func = function(t) {
    pulse=500; pulse_tmin=100; pulse_tmax=150
    ifelse(t > pulse_tmin & t < pulse_tmax, pulse, 0)
}

# Input a list containing vectors of values you want to use for webs
# If `expand = FALSE`, then it provides `n` webs, where `n` is the length of each entry
# in the input list (all vectors in the list must be the same length)
# If `expand = TRUE`, then it provides a web for each combination of entries in the
# parameter list
many_webs <- multi_web(list(Neq = seq(343e3, 343e3 + 50e3, 25e3),
                            Deq = seq(114e3, 114e3 + 10e3, 5e3),
                            iM_func = list(iM_func)),
                       expand = TRUE)
# Solve for unknown values
for (w in many_webs) w$eq_solve()

# ODE solve, add columns for the parameters that differ, then combine
many_outputs <- many_webs %>%
    lapply(function(w) w$ode_solve(tmax = 1000) %>%
               mutate(Neq = w$Neq, Deq = w$Deq)) %>%
    bind_rows

# Plot
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
