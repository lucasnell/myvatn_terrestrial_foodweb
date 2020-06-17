
Code for the following paper:

Transient top-down and bottom-up effects of resources pulsed to multiple trophic levels
-------

by 

Matthew A. McCary, Joseph S. Phillips, Tanjona Ramiadantsoa, Lucas A. Nell,
Amanda R. McCormick, and Jamieson C. Botsch



## File organization:


We've organized much of the code into an R package, `mtf`.
Inside the `R` folder, you'll find...

- `data.R`: a description of the default-parameter estimates present in the `mtf::par_estimates` object
- `equil_pools.R`: The `equil_pools` function that re-estimates equilibrium pool sizes, for use when you adjust parameter values that affect equilibria
- `food_web.R`: The `food_web` function (and required inner functions) that simulates the food web through time, given a set of parameter values
- `zzz.R`: The miscellaneous functions `.onLoad` and `color_pal` that set the `ggplot2` theme and create the color palette for figures, respectively

Obviously `man` contains the documentation for the package, and
`data` contains the data files (the only data object in `mtf` is `par_estimates`).

The folder `data-raw` contains the pre-processed data for the parameter estimates, plus
the code to clean it for the package.
The script `data-raw/par_estimates.R` cleans the data.
The folder `data-raw/parameters` contains the following pre-processed data files:

- `full_pars.csv`: The "final" parameter values that are solved using the specified equilibrium pool sizes and the initially specified parameter values
- `structure_eqs.csv`: The equilibrium pool sizes for different model configurations, based on `full_pars.csv`
- `iN_eqs.csv` The equilibrium pool sizes for different values of the nutrient input rate i_I (with all other values based on `full_pars.csv`)


The `analysis` folder contains the files that use the `mtf` package to produce the
figures and supplemental information file.

- `fig_2-4.R`: Creates figures 2, 3, and 4
- `fig_5.R`: Creates figure 5
- `AppendixS1.Rmd` Creates Appendix S1
- `AppendixS2.Rmd` Creates Appendix S2
- `template.tex`: LaTeX template that helps generate the supplemental information PDF
