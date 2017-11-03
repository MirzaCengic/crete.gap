#' Rarify points
#' Rarify points according to the raster grid mask. Points will be rarified so only one
#' point per grid cell is preserved. Purposfully made for rarifying points in presence-absence modeling. Useful for SDM.
#' @param point_data Input point data to be rarified. SpatialPointsDataFrame.
#' @param raster_data Input grid data to which point_data is rarified. RasterLayer. Has to have overlapping extent.
#'
#' @return SpatialPointsDataFrame with rarified points.
#' @export
#'
#' @examples
rarify_points <- function(point_data, raster_data)
{
  if (missing(raster_data) | missing(point_data))
    {
    stop("Input data missing.")
    }
  # Convert raster to a spatial grid
  raster_grid <- as(raster_data, "SpatialGrid")
  # Get point data
  initial_pt_data <- point_data
  if (
    !identical(raster::crs(raster_data, asText = TRUE), sp::proj4string(initial_pt_data))
  ) {
    initial_pt_data <- sp::spTransform(initial_pt_data, raster::crs(raster_data, asText = TRUE))
  }
  initial_pt_data$grid <- sp::over(initial_pt_data, raster_grid)
  gridlist <- split(initial_pt_data, initial_pt_data$grid)
  # Take one point per grid
  samples <- lapply(gridlist, function(x) x[sample(1:nrow(x), 1, FALSE),])
  # Bind those rows back together in a new data frame
  sampledgrid <- do.call(rbind, samples)
  return(sampledgrid)
}
