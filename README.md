# Reductionism in Rasch Modeling
This R package is associated with a talk given at [ECPA13](http://www.ecpa13.com/) by Peter Edelsbrunner and Fabian Dablander. Idea and major theoretical work was done by Peter Edelsbrunner. For questions about the code, please email me.

## Installation
```{r}
devtools::install_github('dostodabsi/simrasch', build_vignettes = TRUE)
```

## Talk & Simulation Notes
```{r}
library('simrasch')

vignette('simulations') # for discussion of simulation results
vignette('talk') # for the presentation given at ECPA13
```

If you calling the respective vignette doesn't work, simply clone this repository:

```{r}
git clone https://github.com/dostodabsi/simrasch
```

or download the repo manually using the button on the right side of the webpage. In the folder **inst/doc** you will find **talk.html** and **simulations.html**. To view them, simply open them with your browser.
