# Mývatn terrestrial foodweb

Matt McCary

## Original description of model:

This is a nutrient-balance model which was much more dynamically stable than a hybrid 
nutrient/population model that I tried first. There is mass balance, with nutrients 
"leaking" from the system through leaching and entering the system through external 
inputs and midges that enter the detritus pool.

Note that this is an example in which a nutrient balance model is inherently more 
dyanmically stable than a population-level model. This topic (relative stability of 
mass-balance vs. population models) is something Joe and I talked about looking into 
a long time ago.


The model tracks the change of 6 components of a food web through time.
At time t, the total nitrogen pool is denoted by N(t),
the total amount of detritus is denoted by D(t),
the abundance of detritivores by V(t), plants by P(t), herbivores by H(t),
and the predators by P(t).

We assume that at each time step, except for the detritus,
each component has specific mortality rate denoted by mX where X indicate
with component it is. For instance, mD denotes the per unit mortality rate
of the detritivores. A consumed unit of a component X is converted into a
consumer component Y at a rate aXY. For instance, one unit of nutrient is
converted into aND unit of detritus.

Finally, we assume that consumptions by plants and predators follow a
Holling type II functional response with an attack rate of bP and bR, respectively.
At each time step, we add nitrogen and midges denoted by G and M, respectively.
A fraction f of the midge goes to the detritivores whereas the remaining fraction (1-f)
goes to the plant. Under these assumptions, the "abundances" at time t +1 are...
[it then gets into the code]
