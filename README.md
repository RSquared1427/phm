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

## Input Files
### Format of input Hi-C contact matrix

Assume the genomic region of interest contains ![](http://latex.codecogs.com/gif.latex?N) loci. The input file of Hi-C
contact matrix is a ![](http://latex.codecogs.com/gif.latex?N%5Ctimes%20N) symmetric matrix separated by the tab delimiter.
All off-diagonal numbers should be non-negative integers. All diagonal numbers
should be zero. The number in the (i, j) cell is the total number of Hi-C
reads spanning the i-th locus and the j-th locus.

Example:

![](http://latex.codecogs.com/gif.latex?%5Cbegin%7Bmatrix%7D%200%20%26%20197%20%26%20175%20%26154%20%26140%20%26147%20%26102%20%26122%26%20%5Ccdots%5C%5C%20197%26%200%20%26210%20%26138%20%26124%20%2698%20%2684%20%26102%26%20%5Ccdots%5C%5C%20175%26%20210%26%200%26%20348%26%20143%20%26110%20%26115%20%26130%26%20%5Ccdots%5C%5C%20154%20%26138%26%20348%26%200%20%26176%20%26171%20%26202%20%26167%26%20%5Ccdots%5C%5C%20140%26%20124%26%20143%20%26176%26%200%20%26448%20%26248%20%26153%20%26%5Ccdots%5C%5C%20147%26%2098%26%20110%26%20171%26%20448%20%260%20%26303%20%26180%20%26%5Ccdots%5C%5C%20102%26%2084%26%20115%26%20202%26%20248%20%26303%20%260%26%20243%26%20%5Ccdots%5C%5C%20122%26%20102%26%20130%26%20167%26%20153%20%26180%20%26243%20%260%26%20%5Ccdots%5C%5C%20%5Ccdots%20%26%5Ccdots%20%26%5Ccdots%20%26%5Ccdots%26%20%5Ccdots%26%20%5Ccdots%26%20%5Ccdots%26%20%5Ccdots%26%20%5Ccdots%20%5Cend%7Bmatrix%7D)


### Format of input local genomic features

The input file of local genomic features is a ![](http://latex.codecogs.com/gif.latex?N%5Ctimes%206) matrix separated by the table delimiter. For the i-th row:

Column 1: chromosome name for the i-th locus.

Column 2: start position for the i-th locus.

Column 3: end position for the i-th locus.

Column 4: number of restriction enzyme cut fragment ends in the i-th locus (positive integer).

Column 5: mean GC content in the i-th locus (positive real number).

Column 6: mean mappability score in the i-th locus (positive real number).

Example:

![](http://latex.codecogs.com/gif.latex?%5Cbegin%7Bmatrix%7D%2022%20%26%201.6e&plus;07%20%26%201.7e&plus;07%26%20389%26%200.4511%26%200.8831%5C%5C%2022%26%201.7e&plus;07%26%201.8e&plus;07%26%20272%26%200.4534%26%200.8951%5C%5C%2022%26%201.8e&plus;07%26%201.9e&plus;07%26%20218%26%200.5365%26%200.9065%5C%5C%2022%26%201.9e&plus;07%26%202.0e&plus;07%26%20230%26%200.4726%26%200.8704%5C%5C%2022%26%202.0e&plus;07%26%202.1e&plus;07%26%20235%26%200.4366%26%200.9304%5C%5C%2022%26%202.1e&plus;07%26%202.2e&plus;07%26%20400%26%200.4562%26%200.9118%5C%5C%2022%26%202.2e&plus;07%26%202.3e&plus;07%26%20246%26%200.4928%26%200.8788%5C%5C%2022%26%202.3e&plus;07%26%202.4e&plus;07%26%20336%26%200.4595%26%200.9067%5C%5C%20%5Ccdots%20%26%5Ccdots%20%26%5Ccdots%20%26%5Ccdots%26%20%5Ccdots%26%20%5Ccdots%20%5Cend%7Bmatrix%7D)

## Output Files

time.txt: start, end and duration of the PHM running time.

mode_loglike.txt: the posterior model of the log likelihood.

record_loglike.txt: the log likelihood in each iteration of the Gibbs sampler
(a vector with size no of iteration). It can be used to check the convergence of
MCMC chain.

mode_p.txt: the posterior mode of the 3D coordinates (![](http://latex.codecogs.com/gif.latex?N%5Ctimes%203) matrix).

record p.txt: the posterior samplers of the 3D coordinates (number of ![](http://latex.codecogs.com/gif.latex?%5Ctext%7Biteration%7D%5Ctimes%20%023N) matrix). In each row, the 3D coordinates are in the order of: ![](http://latex.codecogs.com/gif.latex?x_1%2C%20y_1%2C%20z_1%2C%20x_2%2C%20y_2%2C%20z_2%2C%5Ccdots%2C%20x_N%2C%20y_N%2C%20z_N).

mode_kappa_tau.txt: the posterior mode of kappas and taus of helixes in the piecewise helical curve ![](http://latex.codecogs.com/gif.latex?%28%5Ckappa%2C%5Ctau%29).

    Column 1: number of helix.

    Column 2: posterior mode of the curvature of the helix(![](http://latex.codecogs.com/gif.latex?%5Ckappa)).
    
    Column 3: posterior mode of the torsion of the helix(![](http://latex.codecogs.com/gif.latex?%5Ctau)).
    
mode_tnb.txt: the posterior mode of t, n, b vectors in Frenet framework at the change points.
