
Be Bayesian my friend

================
This is a course for everyone who want to start learning Bayesian Inference. 

# Program

The course is divided in two parts combining theoretical and practical aspects, and the idea is to teach it in two hours.

**PART I: Introduction to Bayesian inference**. History of Bayes theorem. Bayes theorem. Bayesian inference. Posterior distribution. Credible intervals. Predictive distribution.
 
**PART II: Hierarchical Bayesian models**.  Latent Gaussian models (LGMs). Laplace approximation. Gaussian Markov random fields (GMRFs). Fitting GLMMs using INLA. Structured temporal and spatial random effects.
 
**class2-PART 3: Geostatistics using INLA and SPDE**. Geostatistics in the context of LGMs. The Stochastic partial differential equation (SPDE). 
 

# Software

To take full advantage of the course, it is necessary that everyone has the following programs installed:

- version 4.0.0 of [R](https://cran.r-project.org/) or posterior, and
- [RStudio](https://www.rstudio.com/products/rstudio/download/).

# R packages

This will be the packages required for the course

```r
install.packages(pkgs = c("ggplot2", "gridExtra", "maptools", "rgdal", "spdep", "lattice", "latticeExtra", "viridis", "splancs", "lattice", "fields", "plotKML", "raster", "sp"))

```

The R-INLA package can be downloaded directly from the webpage https://www.r-inla.org/download-install

```r
### --- INLA --- ###
install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
```

Also, other packages from Bioconductor
```r
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
```
