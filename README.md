# PHM Users' Guide

## Introduction
Piecewise helical model (PHM) is a parsimonious, easy to interpret, and robust model for inferring 
three-dimensional (3D) chromosomal structure from Hi-C data. PHM takes the Hi-C contact matrix and 
local genomic features (restriction enzyme cutting frequencies, GC content and sequence uniqueness) 
as input and produces, via MCMC computation, the posterior distribution of three-dimensional (3D)
chromosomal structure.

In Piecewise Helical Model, we assume that chromatin within a topologically associating domain exhibits 
a consensus spatial organization among the cell population. Additionally, it is known from geometry 
that any 3D curve can be uniquely determined by its local curvature and torsion. As a special case, 
a constant curvature and constant torsion lead to a helical curve. Analogous to the fact that any 
continuous function can be approximated by a constant at each point, the curvature and torsion of an 
arbitrary 3D curve can be approximated by piecewise constant functions. Therefore, any continuous 3D 
curve can be approximated by several well-connected helixes, which we refer to as a piecewise helical 
curve. Based on this, we further assume that chromatin folds like a piecewise helical curve in three 
dimension space.

## How to run PHM
The full command is:


```{r, engine='bash', count_lines}
./phm -i heatmap_filename -v covariates_filename -NP number_of_helixpieces -NG number_of_iteration -NT tune_interval -SEED seed 
```

heatmap_filename: string of characters, file name of the input Hi-C contact matrix.

covariates_filename: string of characters, file name of the input local genomic features.

number_of_helixpieces: integers, number of helixes within the piecewise helical curve. Default value is 2.

number_of_iteration: integers, number of Gibbs sampler iterations. Default value is 5,000.

tune_interval: integers, length of tune interval in HMC. Default value is 50.

seed: integers, seed for gsl random number generator. Default value is 1.


Example:
```{r, engine='bash', count_lines}
./PHM -i heatmap.txt -v cov.txt -NP 2 -NG 5000 -NT 50 -SEED 1 
```
