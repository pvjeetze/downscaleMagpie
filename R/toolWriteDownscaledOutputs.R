#' Write downscaled MAgPIE outputs
#'
#' Writes the gridded outputs of the MAgPIE downscaling to the output folder of the run, with the
#' file names and comments of the MAgPIE output scripts. Each product is given as its own argument,
#' so that a script can write and free one product before computing the next.
#'
#' @param inputs List from \code{\link{toolDownscalingInputs}} or \code{\link{toolAgeClassInputs}}
#' (output folder of the run; for peatland its categories and first year).
#' @param land Optional result of \code{\link{downscaleLand}}: writes
#' \code{cell.land_transitions_0.5.mz}, \code{cell.land_diagnostics.csv},
#' \code{cell.land_constraint_excess.csv}, \code{cell.land_transition_change.csv} and, with peatland,
#' \code{cell.peatland_0.5.mz}, \code{cell.peatland_0.5_share.mz} and \code{cell.peatland_excess.csv}.
#' @param cropland Optional result of \code{\link{toolSplitGriddedCropland}}: writes
#' \code{cell.croparea_0.5_share.mz}, \code{cell.land_split_0.5.mz}, \code{cell.land_0.5.mz} and their
#' share files.
#' @param conservation Optional result of \code{\link{downscaleLandConservation}}: writes
#' \code{cell.conservation_land_0.5.mz}.
#' @param bii Optional result of \code{\link{downscaleBII}}: writes \code{cell.bii_0.5.mz}.
#' @param ageClasses Optional result of \code{\link{downscaleAgeClasses}}: writes
#' \code{cell.ageclass_<type>_0.5_share.mz} for secdforest, forestry and other land, and
#' \code{cell.ageclass_diagnostics.csv}.
#' @return \code{NULL}, invisibly.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolDownscalingInputs}}
#' @export
toolWriteDownscaledOutputs <- function(inputs, land = NULL, cropland = NULL, conservation = NULL, bii = NULL,
                                       ageClasses = NULL) {
  outFile <- function(name) file.path(inputs[["outputdir"]], name)

  # Writes a 0.5 degree output with the spatial sets of the outputs.
  .write <- function(x, file, comment, message) {
    base::message(message)
    getSets(x, fulldim = FALSE)[1] <- "x.y.iso"
    write.magpie(x, outFile(file), comment = comment)
  }

  if (!is.null(land)) {
    message("Maximum deviation of the cluster states from the model: ", signif(land[["deviationLr"]], 3), " Mha")
    if (land[["deviationLr"]] > 1e-5) {
      warning("cluster states deviate from the model by up to ", signif(land[["deviationLr"]], 3), " Mha")
    }
    .write(land[["transitionsHr"]], "cell.land_transitions_0.5.mz",
           comment = "unit: Mha per grid-cell and time step, gross transitions by stage of the step",
           message = "Write outputs cell.land_transitions")
    utils::write.csv(land[["diagnostics"]], outFile("cell.land_diagnostics.csv"), row.names = FALSE)
    if (!is.null(land[["excess"]])) {
      utils::write.csv(land[["excess"]], outFile("cell.land_constraint_excess.csv"), row.names = FALSE)
    }
    if (!is.null(land[["transitionChange"]])) {
      utils::write.csv(land[["transitionChange"]], outFile("cell.land_transition_change.csv"), row.names = FALSE)
    }

    if (!is.null(land[["peatlandHr"]])) {
      # peatland from the year in which the model starts to change it
      peatHr <- land[["peatlandHr"]]
      peatHr <- peatHr[, getYears(peatHr, as.integer = TRUE) >= inputs[["peatlandFirstYear"]],
                       inputs[["peatlandCategories"]]]
      utils::write.csv(land[["peatlandExcess"]], outFile("cell.peatland_excess.csv"), row.names = FALSE)
      .write(peatHr, "cell.peatland_0.5.mz", comment = "unit: Mha per grid-cell",
             message = "Write outputs peatland Mha")
      # share of the grid cell area (area of a 0.5 degree grid cell in Mha)
      coords <- mstools::toolGetMappingCoord2Country(pretty = TRUE)
      cellAreaHr <- new.magpie(cells_and_regions = coords[["coords"]],
                               fill = 1e-10 * (111.263 * 1000 * 0.5)^2 * cos(coords[["lat"]] * pi / 180))
      peatShareHr <- peatHr / cellAreaHr
      peatShareHr[is.nan(peatShareHr) | is.infinite(peatShareHr)] <- 0
      .write(peatShareHr, "cell.peatland_0.5_share.mz", comment = "unit: grid-cell area fraction",
             message = "Write outputs peatland share")
    }
  }

  if (!is.null(cropland)) {
    .write(cropland[["cropAreaShareHr"]], "cell.croparea_0.5_share.mz",
           comment = "unit: croparea fractions of total grid-cell", message = "Write outputs cell.croparea_share")
    landSplitHr <- cropland[["landSplitHr"]]
    landHr <- cropland[["landHr"]]
    .write(landSplitHr, "cell.land_split_0.5.mz", comment = "unit: Mha per grid-cell",
           message = "Write outputs cell.land_split")
    .write(landSplitHr / dimSums(landSplitHr, dim = 3), "cell.land_split_0.5_share.mz",
           comment = "unit: grid-cell land area fraction", message = "Write outputs cell.land_split_share")
    .write(landHr, "cell.land_0.5.mz", comment = "unit: Mha per grid-cell", message = "Write outputs cell.land")
    .write(landHr / dimSums(landHr, dim = 3), "cell.land_0.5_share.mz",
           comment = "unit: grid-cell land area fraction", message = "Write outputs cell.land_share")
  }

  if (!is.null(conservation)) {
    .write(conservation, "cell.conservation_land_0.5.mz", comment = "unit: Mha per grid-cell",
           message = "Write outputs cell.conservation_land")
  }

  if (!is.null(bii)) {
    .write(bii, "cell.bii_0.5.mz", comment = "unitless", message = "Write output BII at 0.5\u00b0")
  }

  if (!is.null(ageClasses)) {
    diagnostics <- ageClasses[["diagnostics"]]
    utils::write.csv(diagnostics, outFile("cell.ageclass_diagnostics.csv"), row.names = FALSE)
    message("Largest deviation of the grid cell area from the grid cell land: ",
            signif(max(diagnostics[["gridCellDeviation"]]), 3),
            " Mha; summed deviation of the cluster age classes from the model: ",
            signif(max(diagnostics[["clusterDeviation"]]), 3), " Mha in a year and land type")
    # one file per land type: the memory use of write.magpie() grows with the number of items
    longName <- c(secdforest = "secondary forest", forestry = "forestry", other = "other land")
    for (type in names(longName)) {
      share <- ageClasses[[type]] / dimSums(ageClasses[[type]], dim = 3)
      share[is.nan(share)] <- NA
      getSets(share, fulldim = FALSE) <- c("x.y.iso", "year", "ac")
      file <- outFile(paste0("cell.ageclass_", type, "_0.5_share.mz"))
      message("Write outputs ", basename(file))
      write.magpie(share, file, comment = paste0("unit: share of ", longName[[type]], " area in each age class"))
      rm(share)
      gc(verbose = FALSE)
    }
  }
  return(invisible(NULL))
}
