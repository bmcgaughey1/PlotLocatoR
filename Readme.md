
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PlotLocatoR

<!-- badges: start -->
<!-- badges: end -->

The PlotLocatoR package is designed to help georeference forestry field
plot data. These data should include stem-mapped trees with either XY
locations for trees or a distance and azimuth from plot center to each
tree. The main objetive of teh package is to align stem-mapped tree data
with a lidar-derived anopty height model (CHM). The initial plot
locations need to be “close”. That is, within ~10m of true locations.
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

## Example using brute force approach

The PlotLocatoR package includes sample data for a plot and a canopy
height model (CHM) clipped to cover the area containing the “best” plot
location. The CHM was created from high density lidar data (1000+
pulses/m<sup>2</sup>) collected from a UAS platform. CHM resolution is
0.5m. The initial plot location was collected with a survey-grade GNSS
receiver (Javad Triumph 2). Under test conditions, this receiver can
collect location data accurate to ~1m in forest conditions with a
15-minute occupation time and 1-second epochs. However, for this plot,
the GNSS location was off by several meters.

The objective of this example is to improve the plot location. For the
field data, distance was measured to the center of the stem and azimuth
included local declination. Tree heights were predicted using an
equation developed from field DBH and lidar-derived tree heights. In
practice, dominant and co-dominant trees should produce the best match
but it really depends on the forest type, slope and tree lean. The field
protocol used to collect the field data included a special attribute for
each tree indicating whether or not the tree was visible from above
(LiDAR_visible). This attribute was used to select the trees used for
the location search.

``` r
library(PlotLocatoR)

# load the plot data
P13Field <- read.csv(system.file("extdata", "plot13.csv", package = "PlotLocatoR", mustWork = TRUE), 
                     header = TRUE,
                     stringsAsFactors = FALSE)

# this is the plot location from GNSS in UTM NAD83 zone 10 (EPSG:26910)
plotX <- 418620.0
plotY <- 5276654.0

# create the stem map from the field data and the GNSS location
P13_map <- PlotLocatoR::computeTreePositions(P13Field,
                                                  plotX,
                                                  plotY,
                                                  azLabel = "Azimuth",
                                                  distLabel = "Distance",
                                                  dbhLabel = "DBH_cm"
                                                  )

# add XY locations to plot data
P13 <- data.frame(P13Field, P13_map)

# sort by decreasing height...not really needed unless 
# you want to use the n tallest trees
P13 <- P13[order(-P13$Ht_m),]

# select 10 tallest trees
##P13FieldSubset <- P13[1:10,]

# select trees based on the LiDAR_visible attribute
P13FieldSubset <- P13[P13$LiDAR_visible == "Y",]
#P13FieldSubset <- P13

# read the CHM into a terra SpatRaster object
CHM <- terra::rast(system.file("extdata", "plot13_CHM.tif", package = "PlotLocatoR", mustWork = TRUE))
#CHM <- terra::rast("data/Plot13_CHM_4m_stems.tif")

# set tree heights to 4m plus random variation...necessary or correlation is undefined
# since all pixel values will be equal to the mean value
#P13FieldSubset$truncHt <- 4 + runif(nrow(P13FieldSubset), -0.1, 0.1)

# plot the CHM with the original stem map and plot location
terra::plot(CHM, main = "Plot 13--Original GNSS location")
points(plotX, plotY, pch = 4, col = "red")
points(P13$Xfield, P13$Yfield, pch = 3, col = "black")
points(P13FieldSubset$Xfield, P13FieldSubset$Yfield, pch = 19, col = "black")
legend("bottom"
       , c("GNSS location", "plot trees", "trees used for alignment")
       , col = c("red", "black", "black")
       , pch = c(4, 3, 19)
       )

# perform a brute-force assessment of alternate plot locations
# grid resolution for possible locations will match the cell size of the CHM
sr <- testPlotLocations(
  P13FieldSubset,
  coords = c("Xfield", "Yfield"),
  htField = "Ht_m",
#  htField = "truncHt",
  plotX = plotX,
  plotY = plotY,
  plotRadius = 17.68,
  treeBufferSize = 1,
  crs = 26910,
  CHM = CHM,
  searchRadius = 7,
  progress = FALSE
)
```

![](Readme_files/figure-gfm/BruteForce-1.png)<!-- -->

``` r

# get the index of the best location
bestIndex <- findBestPlotLocation(sr, rule = "minheighterror")

# compute new plot location
plotXAdj <- plotX + sr$offsetX[bestIndex]
plotYAdj <- plotY + sr$offsetY[bestIndex]

# compute new tree locations using adjusted plot location
P13_mapAdj <- PlotLocatoR::computeTreePositions(P13Field,
                                                  plotXAdj,
                                                  plotYAdj,
                                                  azLabel = "Azimuth",
                                                  distLabel = "Distance",
                                                  dbhLabel = "DBH_cm"
                                                  )
# bind XY locations to plot data
P13Adj <- data.frame(P13Field, P13_mapAdj)

# plot everything
terra::plot(CHM, main = "Plot 13--Adjusted location (brute force)")
points(P13Adj$Xfield, P13Adj$Yfield, pch = 19, col = "magenta")
points(plotXAdj, plotYAdj, pch = 9, col = "black")
legend("bottom"
       , c("adjusted trees", "adjusted plot")
       , col = c("magenta", "black")
       , pch = c(19, 9))
```

![](Readme_files/figure-gfm/BruteForce-2.png)<!-- -->

## Example using FFT cross correlation

This example uses cross correlation in the frequency domain using
Fourier transforms. This method is commonly used in image registration
but can be adapted to work with patterns of tree locations. In
operation, a stem map for a large area is created by segmenting
individual trees from a canopy height model to serve as the reference
image. Then a stem map is created using plot measurements (stem-mapped
trees). The two images are aligned by finding the location with the
highest correlation in the frequency domain.

The cross correlation method is MUCH faster than the brute force method.
However, it does not always produce the best location for a plot. This
usually happens when you have dense stands with trees about the same
height as those on the plot. The correlation in these areas is high so
you end up with what appears to be a good location when, in fact the
location is not correct. This problem stems from the fact that the
algorithm is designed for aligning images where you have a fairly
complete picture so lots of features to raise the correlation. the stem
map generated for the plot only includes small areas around each tree so
many fewer features to match.

``` r
newLocation <- findBestPlotLocationCorrelation(
    CHM,
    P13FieldSubset,
    initialX = plotX,
    initialY = plotY,
    searchRadius = 100,
    stemLocationFields = c("Xfield", "Yfield"),
    stemHeightField = "Ht_m",
    CHMbuffer = 1.0,
    stemMapBuffer = 1.0,
    cropCHM = TRUE,
    method = "mask",
    addNoise = FALSE,
    noiseMultiplier = 0.1,
    addSelectiveNoise = "chm",
    corrMethod = "cross",
    includeRasters = FALSE
)
#> The legacy packages maptools, rgdal, and rgeos, underpinning this package
#> will retire shortly. Please refer to R-spatial evolution reports on
#> https://r-spatial.org/r/2023/05/15/evolution4.html for details.
#> This package is now running under evolution status 0

adjP13 <- moveTreesToPlotXY(P13FieldSubset, newLocation[[1]], newLocation[[2]])

terra::plot(CHM, main = "Plot 13--Adjusted location (cross correlation)")
points(plotX, plotY, pch = 4, col = "red")
points(newLocation[[3]], newLocation[[4]], pch = 19, col = "cyan")
points(adjP13$Xfield, adjP13$Yfield, pch = 3, col = "black")
legend("bottom", c("original plot", "adjusted plot", "trees"), col = c("red", "cyan", "black"), 
       pch = c(4, 19, 3))
```

![](Readme_files/figure-gfm/CrossCorrelation-1.png)<!-- -->
