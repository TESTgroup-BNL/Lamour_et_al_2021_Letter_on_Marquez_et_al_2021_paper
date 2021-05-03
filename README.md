# Marquez_et_al_2021_New_Gasex_theory
 The aim of this repository is to implement the Marquez et al. 2021 new gas exchange theory.
 The values of intercellular CO2 (Ci), CO2 at the surface of the leaf (Cs) and the stomatal conductance (gsw) are recalculated from the output of the gas exchange instruments.
 The calculation of those variables using the Marquez et al. 2021 theory uses several variables :
 - Boundary layer conductance to water
 - Leaf temperature
 - Atmospheric pressure
 - Chamber overpressure
 - Sample cell H2O concentration
 - Sample cell CO2 concentration
 - Leaf transpiration
 - Leaf CO2 assimilation rate
It also uses two parameters, the leaf cuticular conductance to water and Beta, the ratio between the cuticular conductance to CO2 and to water. 
Functions are given to recompute the gas exchange variables from the outputs of the LICOR6800 and LICOR6400.

Several exemples of the effect of this new theory are given for Aci data and conductance data.

Márquez, D.A., Stuart-Williams, H. & Farquhar, G.D. An improved theory for calculating leaf gas exchange more precisely accounting for small fluxes. Nat. Plants 7, 317–326 (2021). https://doi.org/10.1038/s41477-021-00861-w
