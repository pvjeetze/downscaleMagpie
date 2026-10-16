#' Downscale MAgPIE land conservation
#'
#' Fits gridded land conservation data (0.5 degree) to the land conservation of a MAgPIE run in each
#' cluster.
#'
#' @details The land conservation of the data (for example WDPA and conservation priority areas,
#' as selected by the protection settings of the run) is fitted in each cluster to the land
#' conservation of the model (\code{pm_land_conservation}, protection and restoration), separately
#' for natural vegetation (primary forest, secondary forest and other land) and pasture.
#'
#' Without \code{landHr}, the data are scaled with one factor per cluster.
#'
#' With \code{landHr}, the land conservation of the model in each cluster is allocated to the grid
#' cells with \code{\link[luscale]{toolAllocateByFilling}}, natural vegetation as one group and pasture: the
#' protected area of the data within the land of its type is the prior, and the land of its type in
#' the grid cell is the cap, because the model can only protect land that exists. Where the model
#' protects less than the data within the land, the protected area is scaled down. Where it protects
#' more (because the data put protection on land of other types and the model keeps the natural
#' vegetation elsewhere in the cluster, or because of restoration), the rest first fills the
#' unprotected land of the grid cells with protection in the data, in proportion to their protected
#' area, and only what their land cannot hold goes to the other grid cells, in proportion to their
#' land. Protected natural vegetation is split over primary forest, secondary forest and other land as
#' the natural vegetation of the grid cell.
#'
#' @param landConservationLr magclass object with the land conservation of each cluster, year and
#' land type by conservation type (\code{pm_land_conservation} of the fulldata.gdx, names
#' "land.consv_type"), in Mha.
#' @param conservationDataHr magclass object with the land conservation of the data for each grid
#' cell, year and land type (Mha), in the years of \code{landConservationLr}.
#' @inheritParams luscale::toolAllocateByCrossEntropy
#' @param landHr Optional magclass object of grid cell land (Mha) by land type, as \code{landHr} of
#' \code{\link{downscaleLand}}. If given, the result is the land conservation that the model
#' achieves on this land (see details).
#' @return magclass object of the land conservation of each grid cell, year and land type (Mha).
#' @author Patrick v. Jeetze, Pascal Sauer
#' @seealso \code{\link{downscaleLand}}
#' @export
downscaleLandConservation <- function(landConservationLr, conservationDataHr, mapping, landHr = NULL) {
  if (is.character(mapping)) mapping <- toolGetMapping(mapping, where = "local")
  cells <- as.character(mapping[["cell"]])
  years <- getYears(landConservationLr)
  stopifnot("conservationDataHr lacks grid cells of mapping" = all(cells %in% getCells(conservationDataHr)),
            "conservationDataHr lacks years of landConservationLr" = all(years %in% getYears(conservationDataHr)))
  dataHr <- conservationDataHr[cells, years, ]
  # land conservation of the model: protection and restoration
  modelLr <- dimSums(landConservationLr, dim = 3.2)

  # the groups natural vegetation and pasture, in the model and in the data
  natveg <- c("primforest", "secdforest", "other")
  modelGroupLr <- mbind(modelLr[, , "past"], setNames(dimSums(modelLr[, , natveg], dim = 3), "natveg"))
  dataGroupHr <- mbind(dataHr[, , "past"], setNames(dimSums(dataHr[, , natveg], dim = 3), "natveg"))

  conservationHr <- dataHr
  if (is.null(landHr)) {
    # the data scaled with one factor per cluster
    scalingHr <- modelGroupLr / toolAggregate(dataGroupHr, mapping, from = "cell", to = "cluster")
    scalingHr[is.na(scalingHr) | is.infinite(scalingHr)] <- 1
    scalingHr <- toolAggregate(scalingHr, mapping, from = "cluster", to = "cell")
    conservationHr[, , "past"] <- scalingHr[, , "past"] * dataHr[, , "past"]
    conservationHr[, , natveg] <- scalingHr[, , "natveg"] * dataHr[, , natveg]
  } else {
    # land conservation achieved on the given land
    stopifnot("landHr lacks grid cells of mapping" = all(cells %in% getCells(landHr)))
    landHr <- landHr[cells, years, ]
    landGroupHr <- mbind(landHr[, , "past"], setNames(dimSums(landHr[, , natveg], dim = 3), "natveg"))
    getSets(landGroupHr) <- getSets(dataGroupHr)

    groups <- getItems(dataGroupHr, dim = 3)
    protectedGroupHr <- dataGroupHr
    for (year in years) {
      landYear <- matrix(as.array(landGroupHr[, year, groups]), length(cells), dimnames = list(cells, groups))
      dataYear <- matrix(as.array(dataGroupHr[, year, groups]), length(cells), dimnames = list(cells, groups))
      modelCluster <- matrix(as.array(modelGroupLr[, year, groups]), ncol = length(groups),
                             dimnames = list(getItems(modelGroupLr, dim = 1), groups))
      # the model protection of each cluster: first in the grid cells with protection in the data, in
      # proportion to it, up to their land, then the rest in the other grid cells, in proportion to their land
      protectedGroupHr[, year, groups] <- luscale::toolAllocateByFilling(modelCluster, pmin(dataYear, landYear),
                                                                         landYear, mapping)
    }

    # protected natural vegetation split as the natural vegetation of the grid cell
    natvegShare <- landHr[, , natveg] / dimSums(landHr[, , natveg], dim = 3)
    natvegShare[is.na(natvegShare)] <- 0
    conservationHr[, , ] <- 0
    conservationHr[, , "past"] <- protectedGroupHr[, , "past"]
    conservationHr[, , natveg] <- collapseDim(protectedGroupHr[, , "natveg"], dim = 3) * natvegShare
  }
  return(conservationHr)
}
