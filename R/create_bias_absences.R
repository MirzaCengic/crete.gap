# TODO: Make it so simple feature species data can be passed
# Add argument checks
# Add option to return dataframe with presences/absences or a spatial file with the same
# Make processing boundary optional
# Remove package loading after imports are defined


# Code is inspired by Fitzpatrick et al. 2013 R code (DOI: 10.1890/ES13-00066.1)

#' Function to create absences that account for sampling bias
#' Create biased pseudo-absence points.
#'
#' @param species_data SpatialPointsDataFrame with species data.
#' @param raster_mask Raster mask with same processing resolution.
#' @param process_boundary SpatialPolygonsDataFrame with the processing extent.
#' @param abs_number Number of created pseudo-absences.
#'
#' @return
#' @export
#'
#' @examples
create_bias_absences <- function(species_data, raster_mask, process_boundary, abs_number, return_spatial = FALSE)
{

  stopifnot(abs_number > 1)

  # load needed packages
  pacman::p_load(sm, PresenceAbsence, sp, raster)

  # Create processing mask
  extent_mask <- raster_mask[[1]] >- 1000 # Check this line
  extent_mask <- mask(extent_mask, process_boundary)

  ## Process species data

  # Get raster cells that have points
  bias <- cellFromXY(extent_mask, species_data)
  # Get unique cell IDs
  cells <- unique(sort(bias))

  kernelXY <- xyFromCell(extent_mask, cells) #Get XY coordinates for the cells, from the mask
  samps <- as.numeric(table(bias)) #Count how many locality points were in the cells

  # Create kernel density surface
  density_surface <- sm.density(kernelXY, weights = samps, display = "none")
  density_points_tmp <- SpatialPoints(expand.grid(x = density_surface$eval.points[, 1], y = density_surface$eval.points[, 2]))
  density_points_tmp <- SpatialPixelsDataFrame(density_points_tmp, data.frame(kde = array(density_surface$estimate,
                                                                                          length(density_surface$estimate))))
  # Convert to raster
  density_raster <- raster(density_points)
  density_raster <- resample(density_raster, extent_mask) # Check this
  density_raster_masked <- density_raster * extent_mask #Excludes areas outside the study area based on the Mask layer
  density_points <- rasterToPoints(density_raster_masked, spatial = FALSE) #(Potential) PseudoAbsences are created here

  #Now to sample for Pseuooabsences Presence Data
  species_absences <- density_points[sample(seq(1:nrow(density_points)), size = abs_number, replace = TRUE,
                                            prob = density_points[, "layer"]), 1:2]

  species_absences_pts <- SpatialPointsDataFrame(species_absences, proj4string = crs(raster_mask), data = as.data.frame(rep(0, abs_number)))
  names(species_absences_pts) <- "PA"
  # species_absences <- SpatialPointsDataFrame(a.sp, data = as.data.frame(rep(1, nrow(species_data))))
  # species_absences$PresAbs <- 0 #Set values for PresAbs column of Absences to 0
  #Single Dataset

  species_data$PA <- rep(1, nrow(species_data))
  species_presences_pts <- species_data[, "PA"]

  species_PA_data <- rbind(species_presences_pts, species_absences_pts)


  if (isTRUE(return_spatial))
  {
    return(species_PA_data)
  } else {
    species_PA_data_df <- as.data.frame(cbind(coordinates(species_PA_data), species_PA_data@data$PA))
    names(species_PA_data_df)[1:2] <- c("X", "Y")
    return(species_PA_data_df)
  }

}
