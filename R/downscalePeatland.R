#' Downscale MAgPIE peatland
#'
#' Disaggregates the peatland categories of a MAgPIE run (module 58) from clusters to grid cells on a
#' given grid cell land.
#'
#' @details
#' In each year, the categories follow the cluster values: each grid cell keeps its peat area, categories
#' that shrink in the cluster shrink in each grid cell, and the freed area goes to the categories that
#' grow (as in \code{\link[luscale]{interpolate2}}). Where \code{rewettingHr} gives a land reduction on
#' drained peat in a year, the drained peat of each grid cell shrinks by the reduction placed there
#' instead. The result is fitted to the cluster values with
#' \code{\link[luscale]{toolAllocateByCrossEntropy}}, with the peat area of
#' each grid cell as cap and the drained peat under cropland, pasture and forestry within the grid cell
#' area of that land type (soft cap).
#'
#' Each year depends only on the peat before it, so \code{\link{downscaleLand}} calls this
#' function year by year together with the land, and the same call on its outputs gives the same
#' peatland.
#'
#' @param landHr magclass object of grid cell land (Mha) by land type, with the years of
#' \code{peatlandLr}. Sub-types (for example "forestry.aff" in \code{landHr} of
#' \code{\link{downscaleLand}}) are summed to their land type.
#' @param peatlandLr magclass object with the peatland area of each cluster, year and category
#' (as from \code{magpie4::PeatlandArea(gdx, level = "cell", sum = FALSE)}, which is cluster level).
#' @param peatlandStartHr magclass object with the peatland area of each grid cell by category before
#' the first year; missing categories are set to 0.
#' @inheritParams luscale::toolAllocateByCrossEntropy
#' @param rewettingHr Optional magclass object with the land reduction on drained peat of each grid cell
#' and year by land type (crop, past, forestry). From \code{transitionsHr} of
#' \code{\link{downscaleLand}}: \code{dimSums(transitionsHr[, , "peatRewetting"], dim =
#' c("stage", "fromType", "to", "toType"))}.
#' @return List with
#' \item{peatlandHr}{magclass object of grid cell peatland by category (Mha)}
#' \item{excess}{data frame of drained peat above the grid cell land, by year and category (Mha)}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{downscaleLand}}
#' @examples
#' peatlandLr <- magclass::new.magpie("A.1", "y2000", c("crop", "rewetted"), fill = c(0.5, 0.5))
#' peatlandStartHr <- magclass::new.magpie(c("c1", "c2"), "y1995", "crop", fill = 0.5)
#' landHr <- magclass::new.magpie(c("c1", "c2"), "y2000", "crop", fill = 1)
#' rewettingHr <- magclass::new.magpie(c("c1", "c2"), "y2000", "crop", fill = c(0, 0.5))
#' # the cropland on peat of grid cell c2 is reduced, so its drained peat is rewetted
#' mapping <- data.frame(cell = c("c1", "c2"), cluster = "A.1")
#' downscalePeatland(landHr, peatlandLr, peatlandStartHr, mapping, rewettingHr)
#' @export
downscalePeatland <- function(landHr, peatlandLr, peatlandStartHr, mapping, rewettingHr = NULL) {
  cells <- getCells(landHr)
  if (is.character(mapping)) mapping <- toolGetMapping(mapping, where = "local")
  cluster <- as.character(mapping[["cluster"]][match(cells, mapping[["cell"]])])
  stopifnot("mapping lacks grid cells of landHr" = !anyNA(cluster))
  clusters <- getItems(peatlandLr, 1)
  years <- getYears(peatlandLr)
  categories <- getNames(peatlandLr)
  drained <- intersect(categories, c("crop", "past", "forestry"))
  landType <- getItems(landHr, dim = 3.1, full = TRUE)

  # peat of the grid cells before the first year (grid cells x categories); magclass inputs matched by grid cell
  current <- matrix(0, length(cells), length(categories), dimnames = list(cells, categories))
  start <- as.array(peatlandStartHr)[, 1, , drop = FALSE]
  common <- intersect(categories, dimnames(start)[[3]])
  current[, common] <- start[cells, 1, common]
  loss <- array(0, c(length(cells), length(years), length(drained)), dimnames = list(cells, years, drained))
  if (!is.null(rewettingHr)) {
    rewetting <- as.array(rewettingHr)
    common <- intersect(drained, dimnames(rewetting)[[3]])
    loss[, , common] <- rewetting[cells, years, common]
  }

  out <- array(0, c(length(cells), length(years), length(categories)), dimnames = list(cells, years, categories))
  excess <- list()
  for (year in years) {
    # land of the drained categories (crop, past, forestry)
    land <- vapply(drained, function(k) rowSums(as.array(landHr[, year, landType == k])[, 1, , drop = FALSE]),
                   numeric(length(cells)))
    land <- matrix(land, length(cells), length(drained), dimnames = list(cells, drained))

    # grid cells with peat, and all grid cells of clusters without such grid cells
    hasPeat <- rowSums(current) > 0
    active <- hasPeat | !(cluster %in% cluster[hasPeat])
    peat <- current[active, , drop = FALSE]
    total <- as.array(peatlandLr[, year, ])[, 1, , drop = FALSE]
    total <- matrix(total, nrow = dim(total)[1], dimnames = dimnames(total)[c(1, 3)])
    total <- total[clusters, categories, drop = FALSE]

    # one step of interpolate2 as prior: in each grid cell, categories shrink by the share by which they
    # shrink in the cluster (drained peat by the land reduction on it, where given), and the freed
    # area goes to the growing categories in proportion to their cluster growth
    previousTotal <- rowsum(peat, cluster[active])[clusters, , drop = FALSE]
    change <- total - previousTotal
    shrinkShare <- pmax(-change, 0) / ifelse(previousTotal > 0, previousTotal, 1)
    shrinkSum <- rowSums(pmax(-change, 0))
    growShare <- pmax(change, 0) / ifelse(shrinkSum > 0, shrinkSum, 1)
    shrink <- peat * shrinkShare[cluster[active], , drop = FALSE]
    if (sum(loss[active, year, ]) > 0) {
      shrink[, drained] <- pmin(loss[active, year, ], peat[, drained, drop = FALSE])
    }
    prior <- peat - shrink + rowSums(shrink) * growShare[cluster[active], , drop = FALSE]

    # caps: the peat area of each grid cell, and drained peat within the land of its type (soft)
    area <- rowSums(peat)
    caps <- list(list(group = rep("peat", length(categories)),
                      cap = matrix(area, ncol = 1, dimnames = list(NULL, "peat"))))
    if (length(drained) > 0) {
      caps[[2]] <- list(group = ifelse(categories %in% drained, categories, NA),
                        cap = land[active, , drop = FALSE], minFactor = 1e-6)
    }
    fallback <- matrix(area, length(area), length(categories))
    result <- suppressWarnings(luscale::toolAllocateByCrossEntropy(total, prior, mapping, caps, fallbackHr = fallback,
                                                                   tol = 1e-6, maxIter = 100))

    current[active, ] <- result[["flowsHr"]]
    out[, year, ] <- current
    excess[[year]] <- data.frame(year = rep(year, length(drained)), category = drained,
                                 area = colSums(pmax(current[, drained, drop = FALSE] - land, 0)))
  }

  peatlandHr <- as.magpie(out, spatial = 1, temporal = 2)
  getSets(peatlandHr, fulldim = FALSE) <- c("x.y.iso", "year", "data")
  return(list(peatlandHr = peatlandHr, excess = do.call(rbind, excess)))
}
