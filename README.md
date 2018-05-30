# gentleHMC

This repository contains materials that appeared in the EMAC 2018 Quantitative Marketing SIG: 
A Gentle Introduction to Estimation Bayesian Models Using Stan. 

The code is organized into an R package, but we haven't tested whether it works as such. 
To install the code as an R package, do the following:

```{r}
# install.packages('devtools')
devtools::install_github('jasonmtroos/gentleHMC')
```

A list of relevant demo code and/or slides can be found via

```{r}
vignette(package = 'gentleHMC')
```

To view one of these "vignettes" (e.g., the vignette called `part_1_slides`), do this:

```{r}
vignette('part_1_slides', package = 'gentleHMC')
```

