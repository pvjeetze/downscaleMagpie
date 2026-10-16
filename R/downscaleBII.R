#' Downscale MAgPIE BII
#'
#' Biodiversity Intactness Index (BII) of each grid cell from the BII coefficients of the MAgPIE
#' clusters and the land shares of the grid cell.
#'
#' @details The land classes of the BII coefficients split pasture into managed pasture and
#' rangeland, other land into primary and secondary other land, and each land class into forested
#' and non-forested land (side layers). Pasture follows the managed pasture and rangeland shares of
#' the side layers; primary and secondary other land follow the initial shares
#' (\code{magpie4::PrimSecdOtherLand()}). Primary other land has BII = 1.
#'
#' @param biiLr magclass object of the BII coefficients of the clusters, as from
#' \code{magpie4::BII(gdx, level = "cell", mode = "auto", landClass = "all", bii_coeff = NULL,
#' side_layers = NULL)}.
#' @param landHr magclass object of grid cell land (Mha) by land type, as \code{landHr} of
#' \code{\link{toolSplitGriddedCropland}}.
#' @param sideLayersHr magclass object of the side layers (\code{luh2_side_layers_0.5.mz}): shares of
#' managed pasture and rangeland, and of forested and non-forested land.
#' @param landIniFile Path of the land-use initialisation of the run (\code{avl_land_full_t_0.5.mz}),
#' for the initial shares of primary and secondary other land.
#' @inheritParams luscale::toolAllocateByCrossEntropy
#' @return magclass object of the BII of each grid cell (unitless).
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolDownscalingInputs}}
#' @export
downscaleBII <- function(biiLr, landHr, sideLayersHr, landIniFile, mapping) {
  if (is.character(mapping)) mapping <- toolGetMapping(mapping, where = "local")
  message("Disaggregating BII values")

  # BII coefficients of the clusters; primary other land has BII = 1
  biiLr <- mbind(biiLr[, , "other", invert = TRUE],
                 setNames(biiLr[, , "other"], c("primother.forested", "primother.nonforested")),
                 setNames(biiLr[, , "other"], c("secdother.forested", "secdother.nonforested")))
  biiLr[, , c("primother.forested", "primother.nonforested")] <- 1
  biiCoeffHr <- toolAggregate(biiLr, mapping, from = "cluster", to = "cell")

  # land shares of the grid cells in the BII land classes, pasture split into managed pasture and rangeland
  landBiiHr <- mbind(landHr[, , "past", invert = TRUE],
                     collapseNames(landHr[, , "past"]) * sideLayersHr[, , c("manpast", "rangeland")])
  landBiiHr <- landBiiHr / (dimSums(landBiiHr, dim = 3) + 1e-10)
  # primary and secondary other land, and potential natural vegetation
  landBiiHr <- magpie4::PrimSecdOtherLand(landBiiHr, landIniFile)
  landBiiHr <- landBiiHr * sideLayersHr[, , c("forested", "nonforested")]

  return(dimSums(landBiiHr * biiCoeffHr, dim = 3, na.rm = TRUE))
}
