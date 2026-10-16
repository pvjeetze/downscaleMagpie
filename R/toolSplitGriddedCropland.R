#' Split gridded cropland
#'
#' Splits the cropland of each grid cell by the cluster shares of a MAgPIE run: into crop area,
#' fallow land and tree cover, and the crop area into crops by water supply. Every grid cell of a
#' cluster gets the same composition; nothing is allocated beyond the cropland of the grid cell.
#'
#' @details The land split has the names "land.type.water": cropland as crop area of food and feed
#' crops (kfo) and of bioenergy crops (kbe, betr and begr), each rainfed and irrigated, fallow land and
#' tree cover; the other land types with their sub-types from \code{landHr} and water "total".
#' Clusters without cropland in a year take the shares of all clusters (at the start, cropland can
#' exist in clusters that have none in the first year). The start of \code{landHr} (y1985) takes the
#' shares of the first year.
#'
#' @param landHr magclass object of grid cell land (Mha) by land type and sub-type, as \code{landHr}
#' of \code{\link{downscaleLand}}.
#' @param cropSplitLr magclass object of cluster cropland by subcategory (crop_area, crop_fallow,
#' crop_treecover), as from \code{magpie4::land(gdx, level = "cell", subcategories = "crop")}.
#' @param cropAreaLr magclass object of cluster crop area by crop and water supply, as from
#' \code{magpie4::croparea(gdx, level = "cell", products = "kcr", product_aggr = FALSE,
#' water_aggr = FALSE)}.
#' @inheritParams luscale::toolAllocateByCrossEntropy
#' @return List with
#' \item{landSplitHr}{magclass object of the land split (Mha, names "land.type.water")}
#' \item{landHr}{magclass object of the land by land type (Mha), the sum of the land split}
#' \item{cropAreaShareHr}{magclass object of the crop area by crop and water supply as share of the
#' grid cell area, in the years of \code{cropAreaLr}}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolDownscalingInputs}}, \code{\link{toolWriteDownscaledOutputs}}
#' @export
toolSplitGriddedCropland <- function(landHr, cropSplitLr, cropAreaLr, mapping) {
  if (is.character(mapping)) mapping <- toolGetMapping(mapping, where = "local")

  # Shares of the items in each cluster and year; clusters without area take the shares of all
  # clusters.
  .clusterShares <- function(x) {
    filled <- x + (dimSums(x, dim = 3) == 0) * dimSums(x, dim = 1)
    shares <- filled / dimSums(filled, dim = 3)
    shares[is.nan(shares)] <- 0
    return(shares)
  }

  # Cluster shares of the model years, the first year also used for the start (y1985), copied to the
  # grid cells of each cluster.
  .gridShares <- function(x) {
    shares <- .clusterShares(x)
    return(toolAggregate(mbind(setYears(shares[, 1, ], "y1985"), shares), mapping, from = "cluster", to = "cell"))
  }

  message("Disaggregating cropland")
  cropHr <- setNames(dimSums(landHr[, , "crop"], dim = 3), NULL)
  cropSplitHr <- .gridShares(cropSplitLr) * cropHr

  message("Disaggregating MAgPIE crop types")
  cropAreaHr <- .gridShares(cropAreaLr) * setNames(cropSplitHr[, , "crop_area"], NULL)
  cropAreaShareHr <- cropAreaHr / setNames(dimSums(landHr, dim = 3), NULL)

  # crop area of food and feed crops and of bioenergy crops, rainfed and irrigated
  kbe <- c("betr", "begr")
  kfo <- setdiff(getNames(cropAreaHr, dim = 1), kbe)
  cropTypesHr <- mbind(
    setNames(dimSums(cropAreaHr[, , kfo][, , "rainfed"], dim = 3), "crop.kfo.rainfed"),
    setNames(dimSums(cropAreaHr[, , kfo][, , "irrigated"], dim = 3), "crop.kfo.irrigated"),
    setNames(dimSums(cropAreaHr[, , kbe][, , "rainfed"], dim = 3), "crop.kbe.rainfed"),
    setNames(dimSums(cropAreaHr[, , kbe][, , "irrigated"], dim = 3), "crop.kbe.irrigated"),
    setNames(cropSplitHr[, , "crop_fallow"], "crop.fallow.total"),
    setNames(cropSplitHr[, , "crop_treecover"], "crop.treecover.total")
  )

  # the land types other than cropland, by sub-type
  nonCropHr <- landHr[, , "crop", invert = TRUE]
  getNames(nonCropHr) <- paste0(getNames(nonCropHr), ".total")
  landSplitHr <- mbind(cropTypesHr, nonCropHr)
  getSets(landSplitHr, fulldim = FALSE)[3] <- "land.type.water"

  # the land split covers all land of the grid cells
  splitDeviation <- max(abs(dimSums(landSplitHr, dim = 3) - dimSums(landHr, dim = 3)))
  if (splitDeviation > 1e-6) {
    warning("land split deviates from the disaggregated land by up to ", signif(splitDeviation, 3),
            " Mha per grid cell")
  }

  return(list(landSplitHr = landSplitHr, landHr = dimSums(landSplitHr, dim = c(3.2, 3.3)),
              cropAreaShareHr = cropAreaShareHr[, getYears(cropAreaLr), ]))
}
