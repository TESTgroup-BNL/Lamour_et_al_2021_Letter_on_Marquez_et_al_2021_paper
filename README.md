## Repository aim

This code was used to analyze data and generate figures associated with Lamour et al. (2021).

Lamour, J., Davidson, K.J., Ely, K.S., Li, Q., Serbin, S.P. and Rogers, A. (2022), New calculations for photosynthesis measurement systems: what's the impact for physiologists and modelers?. New Phytol, 233: 592-598. https://doi.org/10.1111/nph.17762

## Márquez et al. 2021 gas transport model
 
 Márquez et al. 2021 proposed a new gas transport theory between the leaf and the atmosphere, hereafter called the M2021 theory.
 This new theory can be used to re-estimate variables measured by gas exchange instruments (LI-6400, LI-6800,...) such as Ci, while accounting for cuticular conductance.
 This new theory can also be used as a replacement for Fick's law of diffusion for simulating gas exchanges in vegetation models. Usually, leaf gas exchanges are simulated by coupling models of photosynthesis (e.g. Farquhar et al. 1980), stomatal conductance (e.g., Medlyn et al. 2011), and gas transport (e.g., Fick's Law, or Farquhar et al. 1981). However, Fick's Law and Farquhar et al. 1981 fail at representing gas transport through the cuticle.

 ## Repository organization
 
Functions are given to recompute the gas exchange variables from the outputs of the LICOR6800 and LICOR6400 [0_LICOR_Recalculations_functions_M2021.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/0_LICOR_Recalculations_functions_M2021.R). Examples or recalculations of A-Ci curves are made using the script [2_Recalculation_ACi_using_M2021_model.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/2_Recalculation_ACi_using_M2021_model.R).

The effect of the M2021 theory on the A-Ci parameters of the FvCB model (Vcmax, Jmax, Tp, Rd, Farquhar et al. 1980) are obtained using also the script [2_Recalculation_ACi_using_M2021_model.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/2_Recalculation_ACi_using_M2021_model.R).

Simulations of leaf gas exchange for a tropical species are made using the script [3_Simulations_LeafGasEx_Fick_vCF1981_M2021.R](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/3_Simulations_LeafGasEx_Fick_vCF1981_M2021.R). For those simulations, the FvCB photosynthesis model is coupled with a stomatal conductance model (Medlyn et al. 2011). The M2021 theory is used to calculate the gas transport between the leaf and the atmosphere and is compared with the vCF1981 theory (von Caemmerer and Farquhar, 1981) and the Fick's law of diffusion which is often implemented in crops and terrestrial biosphere models. The equations used to make those simulations are given in a PDF file [Equations for the FvCB USO M2021 model.pdf](https://github.com/TESTgroup-BNL/Marquez_et_al_2021_New_Gasex_theory/blob/main/Equations%20for%20the%20FvCB%20USO%20M2021%20model.pdf).

Note that the code uses the [LeafGasExchange](https://github.com/TESTgroup-BNL/LeafGasExchange) package, available on github. https://github.com/TESTgroup-BNL/LeafGasExchange.
The packages 'here' and 'cowplot' are also needed.


## References
von Caemmerer S, Farquhar GD. 1981. Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta 153: 376–387.

Farquhar, G. D., von Caemmerer, S. V., & Berry, J. A. (1980). A biochemical model of photosynthetic CO 2 assimilation in leaves of C 3 species. Planta, 149(1), 78-90.

Márquez, D.A., Stuart-Williams, H. & Farquhar, G.D. An improved theory for calculating leaf gas exchange more precisely accounting for small fluxes. Nat. Plants 7, 317–326 (2021). https://doi.org/10.1038/s41477-021-00861-w

Medlyn BE, Duursma RA, Eamus D, Ellsworth DS, Prentice IC, Barton CVM, Crous KY, Angelis PD, Freeman M, Wingate L. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology 17: 2134–2144.

## DOI

[![DOI](https://zenodo.org/badge/364026620.svg)](https://zenodo.org/badge/latestdoi/364026620)
