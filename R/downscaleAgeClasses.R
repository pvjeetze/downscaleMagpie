#' Downscale MAgPIE age classes
#'
#' Disaggregates the age classes of secondary forest, forestry (by type of module 32) and other
#' land (othernat and youngsecdf of module 35) of a MAgPIE run from clusters to grid cells, along the
#' grid cell land and the grid cell transitions of \code{\link{downscaleLand}}.
#'
#' @details
#' In the first year, each grid cell gets the age composition of its cluster. Each later time step
#' follows the age-class bookkeeping of modules 32 and 35, in the order of the model:
#' \enumerate{
#'   \item disturbance: disturbed secondary forest and afforestation regrow in the establishment
#'   age classes (\code{ac_est}, as many as age classes are shifted in the time step) of the same
#'   grid cell; primary forest lost to disturbance (stage "disturbance") enters secondary forest in
#'   \code{ac_est};
#'   \item recovery (stages "recoveryForestry" and "recovery"): othernat becomes youngsecdf at its
#'   age;
#'   \item ageing: all age classes move on by the length of the time step
#'   (\code{\link[luscale]{toolShiftCohorts}});
#'   \item maturation (stage "maturation"): youngsecdf becomes secondary forest at its age;
#'   \item the other stages: the outflow of each land type to other land types leaves its age
#'   classes; harvested secondary forest and othernat and replanted forestry regrow in
#'   \code{ac_est} of the same grid cell; inflow from other land types (including harvested youngsecdf,
#'   which becomes othernat) enters \code{ac_est}.
#' }
#' Recovery, maturation and outflow are land transitions by age class. Each has two margins, the
#' grid cell total from the grid cell transitions and the cluster flow by age class from the model, and is
#' allocated with \code{\link[luscale]{toolAllocateCohortFlows}}. Disturbance, harvest and replanting stay within a
#' land type and follow the grid cell stock by age class. Outflow of a grid cell beyond its stock (area that
#' enters and leaves a land type in the same time step) passes through without age.
#'
#' The grid cell stocks by age class can then deviate from the age classes of the model in the
#' clusters: where the grid cell transitions differ from the model transitions (urban land from gridded
#' data, see \code{\link[luscale]{toolAllocateTransitions}}), and where the grid cells that the land disaggregation
#' chose for a transition lack the age classes that the model removes. At the end of each time
#' step, the grid cell stocks are therefore fitted to the age classes of the model with
#' \code{\link[luscale]{toolFitCohortStocks}}, keeping the area of each land type in each grid cell. The area
#' of each land type in each grid cell equals the grid cell land, and the age classes of the clusters equal those of the
#' model; the area that the fit moves between age classes is reported.
#'
#' @param dataLr Named list of magclass objects read from the fulldata.gdx, at cluster level:
#' \code{secdforest} (level of ov35_secdforest), \code{forestry} (ov32_land), \code{other}
#' (ov_land_other), \code{secdforestDisturbance} (p35_disturbance_loss_secdf), \code{recovery}
#' (p35_forest_recovery_area), \code{maturation} (p35_maturesecdf), \code{secdforestReduction}
#' (ov35_secdforest_reduction), \code{secdforestHarvest} (ov35_hvarea_secdforest),
#' \code{otherReduction} (ov35_other_reduction), \code{otherHarvest} (ov35_hvarea_other),
#' \code{forestryReduction} (ov32_land_reduction), \code{forestryReplant} (ov32_land_replant) and
#' optionally \code{forestryDisturbance} (p32_disturbance_loss_ftype32). Variables are given as
#' levels.
#' @param landHr magclass object of grid cell land (Mha) by land type and sub-type (names
#' "land.type", optionally followed by further sub-dimensions that are summed), with the years of
#' \code{dataLr[["secdforest"]]}, as \code{landHr} of \code{\link{downscaleLand}} or the land
#' split \code{cell.land_split_0.5.mz}.
#' @param transitionsHr magclass object of the grid cell transitions with names
#' "stage.from.fromType.to.toType", as \code{transitionsHr} of \code{\link{downscaleLand}}.
#' @inheritParams luscale::toolAllocateByCrossEntropy
#' @return List with
#' \item{secdforest, forestry, other}{magclass objects of the grid cell area (Mha) by age class
#' (grid cells x years x age classes); forestry is the sum of its types, other the sum of othernat and
#' youngsecdf}
#' \item{diagnostics}{data frame by year and land type with the largest deviation of the grid cell area
#' from the grid cell land (\code{gridCellDeviation}), the summed absolute deviation of the clusters from the
#' model by age class (\code{clusterDeviation}), the land type area of the model
#' (\code{clusterArea}), the area passing through (\code{passThrough}) and the area that the fit
#' moves between age classes (\code{movedByFit}), all in Mha}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{downscaleLand}}, \code{\link[luscale]{toolAllocateCohortFlows}},
#' \code{\link[luscale]{toolFitCohortStocks}}
#' @export
downscaleAgeClasses <- function(dataLr, landHr, transitionsHr, mapping) {
  # Values of a magclass object by age class in a year (clusters x age classes), for one type where
  # the object has types (names "type.ac").
  .clusterAges <- function(x, year, clusters, ageClasses, type = NULL) {
    if (!is.null(type)) x <- collapseDim(x[, , type], dim = 3.1)
    out <- as.array(x[, year, ])[, 1, , drop = FALSE]
    out <- matrix(out, nrow = dim(out)[1], dimnames = dimnames(out)[c(1, 3)])
    return(out[clusters, ageClasses, drop = FALSE])
  }

  # Cluster stock of a land type by age class in a year (clusters x age classes).
  .clusterAgeStock <- function(dataLr, pool, year, clusters, ageClasses) {
    return(switch(pool,
                  secdforest = .clusterAges(dataLr[["secdforest"]], year, clusters, ageClasses),
                  othernat = , youngsecdf = .clusterAges(dataLr[["other"]], year, clusters, ageClasses, pool),
                  .clusterAges(dataLr[["forestry"]], year, clusters, ageClasses, pool)))
  }

  # Model flows of a land type by age class in a year (clusters x age classes): out, the outflow to
  # other land types (reduction minus harvest or replanting), and within, the flow within the land
  # type to ac_est (harvested secondary forest and othernat, replanted forestry; zero for youngsecdf,
  # whose harvested area becomes othernat, a land transition).
  .clusterAgeOutflow <- function(dataLr, pool, year, clusters, ageClasses) {
    reduction <- switch(pool,
                        secdforest = .clusterAges(dataLr[["secdforestReduction"]], year, clusters, ageClasses),
                        othernat = ,
                        youngsecdf = .clusterAges(dataLr[["otherReduction"]], year, clusters, ageClasses, pool),
                        .clusterAges(dataLr[["forestryReduction"]], year, clusters, ageClasses, pool))
    within <- switch(pool,
                     secdforest = .clusterAges(dataLr[["secdforestHarvest"]], year, clusters, ageClasses),
                     othernat = .clusterAges(dataLr[["otherHarvest"]], year, clusters, ageClasses, pool),
                     youngsecdf = reduction * 0,
                     .replanted(dataLr, pool, year, clusters, reduction))
    return(list(out = pmax(reduction - within, 0), within = within))
  }

  # Replanted forestry area by age class in a year (clusters x age classes): the replanted part of
  # the harvest, in proportion to the reduction by age class.
  .replanted <- function(dataLr, pool, year, clusters, reduction) {
    replanted <- rep(0, length(clusters))
    if (pool %in% getItems(dataLr[["forestryReplant"]], dim = 3)) {
      replanted <- as.array(dataLr[["forestryReplant"]][, year, pool])[clusters, 1, 1]
    }
    reductionSum <- rowSums(reduction)
    return(reduction * (replanted / ifelse(reductionSum > 0, reductionSum, 1)))
  }

  # Adds the area of each grid cell (amount: vector by grid cell) to the first nShift age classes of stock, in
  # equal parts (the establishment age classes ac_est of the time step).
  .toEntry <- function(stock, amount, nShift) {
    entry <- seq_len(nShift)
    stock[, entry] <- stock[, entry] + amount / nShift
    return(stock)
  }

  clusters <- getItems(dataLr[["secdforest"]], 1)
  cells <- getCells(landHr)
  if (is.character(mapping)) mapping <- toolGetMapping(mapping, where = "local")
  cluster <- as.character(mapping[["cluster"]][match(cells, mapping[["cell"]])])
  stopifnot("mapping lacks grid cells of landHr" = !anyNA(cluster))
  ageClasses <- getItems(dataLr[["secdforest"]], dim = 3)
  years <- getYears(dataLr[["secdforest"]])
  yearNumbers <- as.numeric(sub("y", "", years))
  forestryTypes <- getItems(dataLr[["forestry"]], dim = 3.1)

  # land types with age classes, with their label "land.type" in landHr and transitionsHr
  pools <- c("secdforest", forestryTypes, "othernat", "youngsecdf")
  label <- stats::setNames(c("secdforest.total", paste0("forestry.", forestryTypes), "other.othernat",
                             "other.youngsecdf"), pools)

  # grid cell land of each land type (grid cells x years)
  landLabel <- paste(getItems(landHr, dim = 3.1, full = TRUE), getItems(landHr, dim = 3.2, full = TRUE),
                     sep = ".")
  cellArea <- lapply(label, function(x) {
    matrix(as.array(dimSums(landHr[, years, landLabel == x], dim = 3)), ncol = length(years),
           dimnames = list(NULL, years))
  })

  # labels of the grid cell transitions
  stage <- getItems(transitionsHr, dim = 3.1, full = TRUE)
  from <- paste(getItems(transitionsHr, dim = 3.2, full = TRUE), getItems(transitionsHr, dim = 3.3, full = TRUE),
                sep = ".")
  to <- paste(getItems(transitionsHr, dim = 3.4, full = TRUE), getItems(transitionsHr, dim = 3.5, full = TRUE),
              sep = ".")
  presolve <- c("disturbance", "recoveryForestry", "recovery", "maturation")

  # first year: the age composition of the cluster in each grid cell
  stock <- lapply(stats::setNames(pools, pools), function(p) {
    clusterStock <- .clusterAgeStock(dataLr, p, years[1], clusters, ageClasses)
    stockSum <- rowSums(clusterStock)
    share <- clusterStock / ifelse(stockSum > 0, stockSum, 1)
    share[stockSum == 0, ] <- 1 / length(ageClasses)
    x <- share[cluster, , drop = FALSE] * cellArea[[p]][, 1]
    rownames(x) <- cells
    return(x)
  })

  out <- list()
  out[c("secdforest", "forestry", "other")] <- list(array(0, dim = c(length(cells), length(years), length(ageClasses)),
                                                          dimnames = list(cells, years, ageClasses)))
  diagnostics <- list()

  for (i in seq_along(years)) {
    year <- years[i]
    passThrough <- stats::setNames(rep(0, length(pools)), pools)

    if (i > 1) {
      nShift <- (yearNumbers[i] - yearNumbers[i - 1]) / 5
      stopifnot("time steps must be multiples of 5 years" = nShift == round(nShift))
      # grid cell transitions of the year, in the grid cell order of landHr
      flows <- matrix(as.array(transitionsHr[, year, ])[cells, 1, , drop = FALSE], nrow = length(cells))
      .cellFlow <- function(select) rowSums(flows[, select, drop = FALSE])

      # ---- presolve: disturbance, recovery, ageing, maturation ----------------------------------

      disturbance <- .clusterAges(dataLr[["secdforestDisturbance"]], year, clusters, ageClasses)
      disturbed <- luscale::toolAllocateCohortFlows(stock[["secdforest"]], disturbance, mapping)
      stock[["secdforest"]] <- .toEntry(stock[["secdforest"]] - disturbed, rowSums(disturbed) +
                                          .cellFlow(stage == "disturbance" & to == label[["secdforest"]]), nShift)
      if (!is.null(dataLr[["forestryDisturbance"]])) {
        for (k in forestryTypes) {
          disturbance <- .clusterAges(dataLr[["forestryDisturbance"]], year, clusters, ageClasses, k)
          disturbed <- luscale::toolAllocateCohortFlows(stock[[k]], disturbance, mapping)
          stock[[k]] <- .toEntry(stock[[k]] - disturbed, rowSums(disturbed), nShift)
        }
      }

      recovery <- .clusterAges(dataLr[["recovery"]], year, clusters, ageClasses)
      recoveryHr <- .cellFlow(stage %in% c("recoveryForestry", "recovery") &
                                from == label[["othernat"]] & to == label[["youngsecdf"]])
      recovered <- luscale::toolAllocateCohortFlows(stock[["othernat"]], recovery, mapping, recoveryHr, tol = 1e-6)
      stock[["othernat"]] <- stock[["othernat"]] - recovered
      stock[["youngsecdf"]] <- stock[["youngsecdf"]] + recovered

      stock <- lapply(stock, luscale::toolShiftCohorts, nShift = nShift)

      maturation <- .clusterAges(dataLr[["maturation"]], year, clusters, ageClasses)
      matured <- luscale::toolAllocateCohortFlows(stock[["youngsecdf"]], maturation, mapping,
                                                  .cellFlow(stage == "maturation" & from == label[["youngsecdf"]] &
                                                              to == label[["secdforest"]]), tol = 1e-6)
      stock[["youngsecdf"]] <- stock[["youngsecdf"]] - matured
      stock[["secdforest"]] <- stock[["secdforest"]] + matured

      # ---- optimisation: outflow, harvest and replanting, inflow --------------------------------

      for (p in pools) {
        outflow <- .cellFlow(!stage %in% presolve & from == label[[p]] & to != label[[p]])
        inflow <- .cellFlow(!stage %in% presolve & to == label[[p]] & from != label[[p]])
        # area that enters and leaves in the same time step passes through without age
        through <- pmin(pmax(outflow - rowSums(stock[[p]]), 0), inflow)
        passThrough[[p]] <- sum(through)

        model <- .clusterAgeOutflow(dataLr, p, year, clusters, ageClasses)
        left <- luscale::toolAllocateCohortFlows(stock[[p]], model[["out"]], mapping, outflow - through, tol = 1e-6)
        stock[[p]] <- stock[[p]] - left
        regrowth <- luscale::toolAllocateCohortFlows(stock[[p]], model[["within"]], mapping)
        stock[[p]] <- .toEntry(stock[[p]] - regrowth, rowSums(regrowth) + inflow - through, nShift)
      }
    }

    # ---- fit to the age classes of the model, store and diagnose -------------------------------

    for (p in pools) {
      clusterStock <- .clusterAgeStock(dataLr, p, year, clusters, ageClasses)
      tracked <- stock[[p]]
      stock[[p]] <- luscale::toolFitCohortStocks(tracked, clusterStock, mapping, tol = 1e-6)
      diagnostics[[length(diagnostics) + 1]] <- data.frame(
        year = year, type = p,
        gridCellDeviation = max(abs(rowSums(stock[[p]]) - cellArea[[p]][, i])),
        clusterDeviation = sum(abs(rowsum(stock[[p]], cluster)[clusters, , drop = FALSE] - clusterStock)),
        clusterArea = sum(clusterStock), passThrough = passThrough[[p]],
        movedByFit = sum(abs(stock[[p]] - tracked)) / 2
      )
    }
    out[["secdforest"]][, i, ] <- stock[["secdforest"]]
    out[["forestry"]][, i, ] <- Reduce(`+`, stock[forestryTypes])
    out[["other"]][, i, ] <- stock[["othernat"]] + stock[["youngsecdf"]]
  }

  result <- list()
  for (group in names(out)) {
    result[[group]] <- as.magpie(out[[group]], spatial = 1, temporal = 2)
    getSets(result[[group]], fulldim = FALSE) <- c("x.y.iso", "year", "ac")
    out[[group]] <- NULL
  }
  result[["diagnostics"]] <- do.call(rbind, diagnostics)
  return(result)
}
