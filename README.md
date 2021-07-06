# Aim of the repository
 
 The aim of this repository is to implement the Marquez et al. 2021 new gas transport theory between the leaf and the atmosphere, hereafter called the M2021 theory.
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

Functions are given to recompute the gas exchange variables from the outputs of the LICOR6800 and LICOR6400 [0_LICOR_Recalculations_functions_M2021.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/0_LICOR_Recalculations_functions_M2021.R). Examples or recalculations are made using the script [2_Recalculation_ACi_using_M2021_model.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/2_Recalculation_ACi_using_M2021_model.R).

The effect of the MSWF theory on the ACi parameters are obtained using also the script [2_Recalculation_ACi_using_M2021_model.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/2_Recalculation_ACi_using_M2021_model.R).

Simulations of leaf gas exchange for a tropical species are made using the script [3_Simulations_LeafGasEx_Fick_vCF1981_M2021.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/3_Simulations_LeafGasEx_Fick_vCF1981_M2021.R). For those simulations the Farquhar et al. 1980 photosynthesis model is coupled with a stomatal conductance model (Medlyn et al. 2011). The M2021 theory is used to calculate the gas transport between the leaf and the atmosphere and is compared with the vCF theory (von Caemmerer and Farquhar, 1981) and the Fick's law of diffusion which is often implemented in crops and terrestrial biosphere models. The equations used to make those simulations are given in a PDF file [Equations for the FvCB USO M2021 model.pdf](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/Equations%20for%20the%20FvCB%20USO%20M2021%20model.pdf).

Note that the code uses the [LeafGasExchange](https://github.com/TESTgroup-BNL/LeafGasExchange) package, available on github. https://github.com/TESTgroup-BNL/LeafGasExchange
It also uses the packages 'here' and 'cowplot'.


## References
von Caemmerer S, Farquhar GD. 1981. Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153: 376–387.

Farquhar, G. D., von Caemmerer, S. V., & Berry, J. A. (1980). A biochemical model of photosynthetic CO 2 assimilation in leaves of C 3 species. Planta, 149(1), 78-90.

Márquez, D.A., Stuart-Williams, H. & Farquhar, G.D. An improved theory for calculating leaf gas exchange more precisely accounting for small fluxes. Nat. Plants 7, 317–326 (2021). https://doi.org/10.1038/s41477-021-00861-w

Medlyn BE, Duursma RA, Eamus D, Ellsworth DS, Prentice IC, Barton CVM, Crous KY, Angelis PD, Freeman M, Wingate L. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology 17: 2134–2144.
