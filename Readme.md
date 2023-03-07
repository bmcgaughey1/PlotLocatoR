
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PlotLocatoR

<!-- badges: start -->
<!-- badges: end -->

The PlotLocatoR package is designed to help georeference forestry field
plot data. These data should include stem-mapped trees with either XY
locations for trees or a distance and azimuth from plot center to each
tree. The main objetive of teh package is to align stem-mapped tree data
with a lidar-derived anopty height model (CHM). The initial plot
locations need to be “close”. That is, within \~10m of true locations.
Locations can be farther off but the brute force search method
implemented in the package will be extremely slow when searching for a
pattern match over a large area.

## Installation

This package is only distributed from my GitHub account. It may make it
to CRAN at some point but no promises.

~~You can install the released version of PlotLocatoR from
[CRAN](https://CRAN.R-project.org) with:~~

``` r
#install.packages("PlotLocatoR")
```

**PlotLocatoR** is currently available as a development version only.
The **devtools** package is required to install **PlotLocatoR**. If you
have not previously used **devtools**, use the commented line of code to
install the package. Note that this will also install several additional
packages needed for devtools. If you do not want the vignettes, set
*build_vignettes = FALSE*.

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
devtools::install_github("bmcgaughey1/PlotLocatoR", build_vignettes = TRUE)
```
