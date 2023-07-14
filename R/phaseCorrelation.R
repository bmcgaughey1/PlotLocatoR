# useful links for the phase correlation approach
#
# Code to do image correlation to help adjust plot locations
# look for Fourier Mellin transform
# https://stackoverflow.com/questions/57801071/get-rotational-shift-using-phase-correlation-and-log-polar-transform
# https://sthoduka.github.io/imreg_fmt/docs/phase-correlation/
# https://sthoduka.github.io/imreg_fmt/docs/overall-pipeline/
# https://stackoverflow.com/questions/30630632/performing-a-phase-correlation-with-fft-in-r
# https://en.wikipedia.org/wiki/Phase_correlation



#'
#' xcoor3d: perform image registration using phase correlation. This function
#' was taken from the imagefx package and modified to better handle NA values
#' in the image matrices.
#'
#' @param img1 matrix: image representing reference set of objects
#' @param img2  matrix: image representing movable set of objects
#'
#' @return \code{list} containing the following items:
#'    max.shifts shifts to align img2 with img1
#'    max.corr correlation values
#'    max.count
#'    corr.mat
#'
# ' @examples
xcorr3d <- function(
    img1,
    img2
) {
  ## normalize by subtracting the mean
  img1 <- img1-mean(img1, na.rm = TRUE)
  img2 <- img2-mean(img2, na.rm = TRUE)

  ## go to the frequency domain and take conjugate of first image
  IMG1c <- Conj(stats::fft(img1))
  IMG2 <- stats::fft(img2)

  ## calculate the cross power spectrum
  R <- (IMG1c * IMG2)

  ## back to the spatial domain
  r.shift <- stats::fft(R,inverse=TRUE)/length(R)

  ## rearrange so zero frequency is in middle
  r <- imagefx::fftshift(Re(r.shift))

  ## normalize the r matrix
  r.norm <- r/(length(img1))/(stats::sd(as.vector(img1), na.rm = TRUE)*stats::sd(as.vector(img2), na.rm = TRUE))

  ## what is the correlation value associated with best shift
  max.cor <- max(r, na.rm = TRUE)

  max.count <- sum(r == max.cor)

  ## normalize the correlation value
  max.cor.norm <- max.cor/(length(img1))/(stats::sd(as.vector(img1), na.rm = TRUE)*stats::sd(as.vector(img2), na.rm = TRUE))

  ## find where the zero frequency is
  if(nrow(r)%%2 == 1) { zero.freq.x = -1.5 }
  if(nrow(r)%%2 == 0) { zero.freq.x = -1   }
  if(ncol(r)%%2 == 1) { zero.freq.y = -1.5 }
  if(ncol(r)%%2 == 0) { zero.freq.y = -1   }


  ## find the maximum value relative to the middle
  max.inds.abs <- which(r==max.cor,arr.ind=TRUE)-(dim(r)/2)

  ## adjust the maximum indices according to the true zero frequency
  max.inds <- max.inds.abs + c(zero.freq.x,zero.freq.y)

  ## if the input images are made of 0s, the r matrix could be filled
  ## with NaN and there is no max
  if(length(max.inds)==0){max.inds=matrix(c(0,0),nrow=1)}

  ## sometimes there is more than one max (normally in blank images)
  ## so take only the first max ind
  max.inds <- max.inds[1,]

  ## create a list to hold the max correlation value, its indicies after
  ## shifting according to the zero frequency, and the original correlation matrix
  return.list = list()
  return.list$max.shifts <- max.inds
  return.list$max.corr <- max.cor.norm
  return.list$max.count <- max.count
  return.list$corr.mat <- r.norm

  ## return the above list
  return(return.list)
}


#' Compute the best location for the stem map relative to the CHM using phase correlation.
#'
#' findBestPlotLocationPhaseCorrelation identifies tree approximate objects (TAOs) in the \code{CHM},
#' buffers the TAOs using \code{CHMbuffer}, and creates a raster image of the objects. \code{method}
#' controls the way the CHM is used for the image. When \code{method="buffer"}, the circular area
#' associated with each TAO is assigned the same height as the TAO. When \code{method="mask"}, the \code{CHM}
#' is masked using the TAO circles. The \code{stemMap} is used to create a second raster image representing
#' the trees on a plot. Each tree location is buffered using \code{stemMapBuffer} and all pixels in
#' the buffer are assigned the tree height.
#'
#' The shift of the stem map image that best aligns with the CHM image is determined using phase correlation
#' and converted into a new coordinate for the reference point. This new coordinate can be used to translate
#' the trees in the \code{stemMap}.
#'
#' @param CHM SpatRaster representing the canopy height surface
#' @param stemMap sf point object with tree positions based on \code{initialX,initialY} location
#' @param initialX initial reference point location
#' @param initialY initial reference point location
#' @param searchRadius radius defining area for possible locations. \code{searchRadius} will be used
#'  to crop the \code{CHM} to limit the area for the search.
#' @param stemLocationFields character list with the field names in \code{stemMap} containing the X
#'  and Y values for the tree locations.
#' @param stemHeightField character name of the field in \code{stemMap} containing the tree height.
#' @param CHMbuffer radius for each TAO detected in \code{CHM}
#' @param stemMapBuffer radius for each tree in the \code{stemMap}
#' @param cropCHM boolean to indicate whether or not the CHM should be cropped to the search area
#'  defined by the \code{(initialX,initialY)} and \code{searchRadius}.
#' @param method character string defining the method used to handle the \code{CHM}. Possible
#'  values are "buffer" and "mask".
#'
#' @return invisible list containing the (X,Y) for the best plot location
#' @export
#'
# ' @examples
findBestPlotLocationPhaseCorrelation <- function(
    CHM,
    stemMap,
    initialX,
    initialY,
    searchRadius,
    stemLocationFields = c("X", "Y"),
    stemHeightField = "Ht",
    CHMbuffer = 1.0,
    stemMapBuffer = 1.0,
    cropCHM = TRUE,
    method = "buffer"
) {
  NAReplaceValue <- 0

  # crop CHM...if requested
  if (cropCHM) {
    CHM <- terra::crop(CHM
                , terra::ext(initialX - searchRadius, initialX + searchRadius, initialY - searchRadius, initialY + searchRadius)
                , snap = "near")
  }

  # replace NA values
  CHM[is.na(CHM)] <- NAReplaceValue

  # find TAOs
  # generate local maxima
  TAO <- lidR::locate_trees(CHM, lidR::lmf(ws = 5))
  TAO$Ht <- TAO$Z

  TAOxy <- sf::st_coordinates(TAO)

  TAO$X <- TAOxy[,1]
  TAO$Y <- TAOxy[,2]

  TAOtrees <- sf::st_as_sf(TAO,
                           coords = c("X", "Y"),
                           remove = FALSE,
                           crs = terra::crs(CHM))

  # buffer trees extracted from the CHM
  TAOtreesBuf <- sf::st_buffer(TAOtrees, CHMbuffer)

  # rasterize the buffered TAOs
  r <- terra::rast(terra::vect(TAOtreesBuf)
                   , extent = terra::ext(CHM)
                   , resolution = terra::res(CHM)[1]
  )

  TAOPM <- terra::rasterize(terra::vect(TAOtreesBuf)
                            , r
                            , field = "Ht"
                            , extent = terra::ext(CHM)
                            , resolution = terra::res(CHM)[1]
                            , background = 0
  )

  if (method == "mask") {
    # create CHM masked using TAOs extracted from CHM...this will keep the top of the crowns
    referenceCHM <- terra::mask(CHM, TAOPM, maskvalues = 0)
    TAOCHM[is.na(TAOCHM)] <- NAReplaceValue
  } else {
    # keep original CHM values within TAO buffers
    referenceCHM <- TAOPM
  }

  # create stem map object
  trees <- sf::st_as_sf(stemMap
                        , coords = stemLocationFields
                        , remove = FALSE
                        , crs = terra::crs(stemMap))

  # buffer trees...entire buffer around tree will have the same height
  treesBuf <- sf::st_buffer(trees, stemMapBuffer)

  r <- terra::rast(terra::vect(treesBuf)
                   , extent = terra::ext(CHM)
                   , resolution = terra::res(CHM)[1]
  )

  PM <- terra::rasterize(terra::vect(treesBuf)
                         , r
                         , field = stemHeightField
                         , extent = terra::ext(CHM)
                         , resolution = terra::res(CHM)[1]
                         , background = 0
  )

  # create matrix objects
  i1 <- terra::as.matrix(TAOCHM, wide = TRUE)
  i2 <- terra::as.matrix(PM, wide = TRUE)

  # do correlation to get offset for stem map
  shifts <- xcorr3d(i1,i2)

  # shift vector is rotated 90 degree clockwise
  offsetX <- -shifts$max.shifts[2] * terra::res(CHM)[1]
  offsetY <- shifts$max.shifts[1] * terra::res(CHM)[1]

  invisible(return(list(initialX + offsetX, initialY + offsetY)))
}