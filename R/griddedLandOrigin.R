#' Track the origin of gridded land
#'
#' Area of land types in each grid cell by the type it came from, along gross grid cell transitions (as from
#' \code{\link{downscaleLand}}). Area that enters a tracked type from another type takes that type
#' as origin, area that moves between tracked types keeps its origin, and the initial area of a type
#' has the type itself as origin. Outflow from a tracked type takes the area that entered last first
#' (\code{order = "newest"}), or takes from all area in proportion (\code{order = "proportional"}).
#'
#' @param transitionsHr magclass object of the gross grid cell transitions per year with names
#' "stage.from.to". The stages of a year are replayed in the order of their items. The outflow of a
#' stage comes from the area at the start of the stage, and only beyond it from the inflow of the
#' stage (area passing through).
#' @param startHr magclass object of the grid cell area by type in the year before the first year of
#' \code{transitionsHr} (one year, with the grid cells of \code{transitionsHr}).
#' @param types Types to track (items of "from" and "to" in \code{transitionsHr}).
#' @param order "newest" or "proportional", see Description.
#' @return magclass object of the grid cell area by type and origin (names "type.origin") in the year of
#' \code{startHr} and the years of \code{transitionsHr}.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{downscaleLand}}
#' @examples
#' initial <- magclass::new.magpie("A.1", "y2000", c("crop", "other"), c(2, 1))
#' transitions <- magclass::new.magpie("A.1", c("y2010", "y2020"),
#'                                     c("change.crop.other", "change.other.crop"), c(0.5, 0, 0, 0.7))
#' magclass::getSets(transitions, fulldim = FALSE)[3] <- "stage.from.to"
#' griddedLandOrigin(transitions, initial, "other")
#' @export
griddedLandOrigin <- function(transitionsHr, startHr, types, order = "newest") {
  # One stage of griddedLandOrigin (area: list by tracked type; flows: grid cells x transitions of the stage). The
  # outflow comes from the area at the start of the stage, and beyond it from the inflow of the stage
  # (area passing through); it is carried to tracked types with its origin. Inflow enters at entry.
  .originStage <- function(area, flows, from, to, origins, order, entry) {
    types <- names(area)
    nCell <- nrow(flows)

    # inflow from types that are not tracked
    inflow <- list()
    carried <- list()
    for (j in types) {
      inflow[[j]] <- matrix(0, nCell, length(origins))
      carried[[j]] <- matrix(0, nCell, length(origins))
      for (m in which(to == j & !from %in% types)) {
        origin <- match(from[m], origins)
        inflow[[j]][, origin] <- inflow[[j]][, origin] + flows[, m]
      }
    }

    shortfall <- 0
    for (k in intersect(types, from)) {
      leaving <- which(from == k)
      total <- rowSums(flows[, leaving, drop = FALSE])
      taken <- .takeArea(area[[k]], total, order)
      area[[k]] <- taken[["area"]]
      inflowSum <- rowSums(inflow[[k]])
      passing <- inflow[[k]] * pmin(taken[["shortfall"]] / ifelse(inflowSum > 0, inflowSum, 1), 1)
      inflow[[k]] <- inflow[[k]] - passing
      shortfall <- max(shortfall, taken[["shortfall"]] - rowSums(passing))
      for (m in leaving[to[leaving] %in% types]) {
        carried[[to[m]]] <- carried[[to[m]]] + (taken[["mix"]] + passing) * (flows[, m] / ifelse(total > 0, total, 1))
      }
    }

    for (j in types) {
      area[[j]][, , entry] <- area[[j]][, , entry] + inflow[[j]] + carried[[j]]
    }
    return(list(area = area, shortfall = shortfall))
  }

  # Takes the grid cell outflow total from area (grid cells x origins x years of entry), the newest entries first
  # or in proportion; returns the remaining area, the origin mix of the outflow (grid cells x origins) and
  # the outflow above the area of each grid cell.
  .takeArea <- function(area, total, order) {
    stock <- rowSums(area)
    if (order == "proportional") {
      taken <- area * pmin(total / ifelse(stock > 0, stock, 1), 1)
    } else {
      taken <- array(0, dim(area))
      remaining <- total
      for (e in rev(seq_len(dim(area)[3]))) {
        entry <- rowSums(area[, , e, drop = FALSE])
        take <- pmin(remaining, entry)
        taken[, , e] <- area[, , e] * (take / ifelse(entry > 0, entry, 1))
        remaining <- remaining - take
      }
    }
    return(list(area = area - taken, mix = rowSums(taken, dims = 2), shortfall = pmax(total - stock, 0)))
  }

  # Area by grid cell and origin of each tracked type (summed over the years of entry), as one matrix
  # grid cells x (types and origins).
  .originArea <- function(area) {
    return(do.call(cbind, lapply(area, rowSums, dims = 2)))
  }

  stopifnot("order must be \"newest\" or \"proportional\"" = order %in% c("newest", "proportional"),
            "startHr must have one year" = nyears(startHr) == 1,
            "grid cells of startHr and transitionsHr differ" = identical(getCells(startHr), getCells(transitionsHr)))

  stage <- getItems(transitionsHr, "stage", full = TRUE)
  from <- getItems(transitionsHr, "from", full = TRUE)
  to <- getItems(transitionsHr, "to", full = TRUE)
  years <- getYears(transitionsHr)
  origins <- unique(c(types, from[to %in% types & !from %in% types]))
  nCell <- ncells(startHr)

  # area by grid cell, origin and year of entry, for each tracked type
  area <- lapply(stats::setNames(types, types), function(k) {
    x <- array(0, c(nCell, length(origins), length(years) + 1))
    x[, match(k, origins), 1] <- as.vector(startHr[, , k])
    return(x)
  })

  out <- array(0, c(nCell, length(years) + 1, length(types) * length(origins)),
               dimnames = list(getCells(startHr), c(getYears(startHr), years),
                               paste(rep(types, each = length(origins)), origins, sep = ".")))
  out[, 1, ] <- .originArea(area)

  shortfall <- 0
  for (i in seq_along(years)) {
    flowsYear <- matrix(transitionsHr[, years[i], ], nrow = nCell)
    for (s in unique(stage)) {
      n <- which(stage == s & from != to & (from %in% types | to %in% types))
      step <- .originStage(area, flowsYear[, n, drop = FALSE], from[n], to[n], origins, order, i + 1)
      area <- step[["area"]]
      shortfall <- max(shortfall, step[["shortfall"]])
    }
    out[, i + 1, ] <- .originArea(area)
  }

  if (shortfall > 1e-6) {
    warning("outflow exceeds the tracked area by up to ", signif(shortfall, 3), " in a grid cell")
  }

  result <- as.magpie(out, spatial = 1, temporal = 2)
  getSets(result, fulldim = FALSE) <- c(getSets(startHr, fulldim = FALSE)[1], "year", "type.origin")
  return(result)
}
