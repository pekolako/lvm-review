# **A review of generalized linear latent variable models and related computational approaches**

This project contains code to reproduce the results of the paper with the same name by Korhonen, P., Nordhausen, K. and Taskinen, S.

The goal of the project was to review and compare different algorithms in the context of generalized linear latent variable models (GLLVMs).

**Summary**

All code here is written in R and requires the packages ltm, boral, gllvm,
glmmTMB, gmf, vegan, MASS, doMPI, xtable, corrplot and gcus.

The file 'simulate_bernoulli.R' reproduces the simulation setup 1 in Section 4.1, while 'simulate_poisson.R' reproduces the simulation setup 2.

'example.R' reproduces the example in Section 4.2, including Figures 3, 4 and 5.

**Authors**

Korhonen, P., Nordhausen, K. and Taskinen, S.


**Further reading**

Note that more details about fitting GLLVMs and worked examples using the package 'gllvm' are given in Niku et al. (2019) and in the package vignettes.

**License**

GNU GPLv3

**References**

Niku, J., Hui, F. K. C.,  Taskinen, S.and D. I. Warton D. I. (2019): gllvm: Fast analysis of multivariate abundance data with generalized
linear latent variable models in R. Methods in Ecology and Evolution, 10, 2173-2182.
