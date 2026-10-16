#' Downscale MAgPIE land
#'
#' Disaggregates the land of a MAgPIE run from clusters to grid cells with
#' \code{\link[luscale]{toolAllocateTransitions}}. The cluster transitions of each time step are allocated
#' stage by stage in the order of the model: the presolve shifts of module 35 (disturbance,
#' recovery, maturation), area passing through forestry, the cropland relocation of the SNV policy,
#' the reduction of managed land on drained peat and the transitions of the optimisation. Forestry
#' and other land are split into the sub-types of modules 32 and 35.
#'
#' @details
#' In the optimisation, urban land follows the 0.5 degree data, scaled to the cluster urban area of
#' the model; the other transitions are changed so that the cluster states stay those of the model
#' (see \code{\link[luscale]{toolAllocateTransitions}}). The \code{constraints} (for example potential forest,
#' available cropland and protected land) apply in all stages.
#' Recovery is placed on the othernat inflow of the previous step (optimisation, SNV relocation and
#' land reduction on drained peat), forestry abandonment first.
#' The cropland relocation of the SNV policy (\code{p29_snv_relocation}, part of the transitions
#' from cropland to secondary forest and other land, \code{q29_land_snv_trans}) is placed in
#' proportion to the part of \code{snvTargetHr} not relocated in earlier steps.
#'
#' Peatland is updated with the land year by year by \code{\link{downscalePeatland}}. The
#' decline of drained peat under cropland, pasture and forestry in a cluster
#' (\code{q58_peatlandMan}: managed land reduction on drained peat, followed by rewetting or
#' abandonment) is a stage of land reduction (peatRewetting) placed in proportion to the drained
#' peat of the grid cells in the previous year; in each grid cell, the drained peat then shrinks by the land
#' reduction placed there.
#'
#' @param dataLr Named list of magclass objects read from the fulldata.gdx, at cluster level:
#' \code{transitions} (level of ov_lu_transitions), \code{land} (ov_land), \code{forestry}
#' (ov_land_forestry), \code{other} (ov_land_other), \code{disturbance}
#' (p35_disturbance_loss_primf), \code{recovery} (p35_forest_recovery_area), \code{maturation}
#' (p35_maturesecdf), \code{forestryExpansion} (ov_landexpansion_forestry),
#' \code{forestryReduction} (ov_landreduction_forestry), \code{otherReduction}
#' (ov35_other_reduction), \code{otherHarvest} (ov35_hvarea_other), \code{transitionsUpper} (upper bound
#' of ov_lu_transitions: the transitions that the model allows) and optionally
#' \code{snvRelocation} (p29_snv_relocation). Variables are given as levels.
#' @param startHr magclass object with the land of each grid cell before the first year (one year) in
#' the types of the land-use initialisation (crop, past, range, forestry, primforest, secdforest,
#' urban, primother, secdother).
#' @inheritParams luscale::toolAllocateByCrossEntropy
#' @param urbanHr Optional magclass object with the urban land of each grid cell and year (one data
#' column).
#' @param constraints Optional named list of constraints as in \code{\link[luscale]{toolAllocateTransitions}}.
#' \code{types} can name the main land types of MAgPIE (forestry and other stand for their
#' sub-types) and sub-types (for example youngsecdf); bounds can also be magclass objects (grid cells x
#' years, one data column).
#' @param snvTargetHr Optional numeric vector with the cropland of each grid cell that requires
#' relocation under the SNV policy (in the order of the grid cells of \code{startHr}); without it, the
#' relocation follows the cropland area.
#' @param peatlandLr Optional magclass object with the peatland area of each cluster, year and
#' category (as from \code{magpie4::PeatlandArea(gdx, level = "cell", sum = FALSE)}, which is cluster
#' level), in the years of \code{dataLr[["transitions"]]}.
#' @param peatlandStartHr magclass object with the peatland area of each grid cell by category before
#' the first year, needed with \code{peatlandLr}; missing categories are set to 0.
#' @return List with
#' \item{landHr}{magclass object of grid cell land (Mha) by land type and sub-type (names "land.type", type
#' "total" for land types without sub-types), the state before the first year labelled y1985}
#' \item{transitionsHr}{magclass object of the gross grid cell transitions per year (Mha) by stage, land
#' type and sub-type (names "stage.from.fromType.to.toType", type "total" as in \code{land}, stages
#' in the order of the land accounting); the transitions between the main types are its sum over
#' "stage", "fromType" and "toType"}
#' \item{diagnostics, excess, transitionChange}{as returned by \code{\link[luscale]{toolAllocateTransitions}}}
#' \item{deviationLr}{maximum deviation of the cluster sums from the model states (Mha)}
#' \item{peatlandHr}{magclass object of grid cell peatland by category (Mha), or \code{NULL}}
#' \item{peatlandExcess}{area of drained peat above the grid cell land, by year and category (Mha)}
#' @author Patrick v. Jeetze
#' @seealso \code{\link[luscale]{toolAllocateTransitions}}, \code{\link{downscalePeatland}}
#' @export
downscaleLand <- function(dataLr, startHr, mapping, urbanHr = NULL, constraints = NULL, # nolint: cyclocomp_linter
                          snvTargetHr = NULL, peatlandLr = NULL, peatlandStartHr = NULL) {
  # =================================================================================================
  # Land types and cluster transitions of the run
  # =================================================================================================

  # Transitions between sub-types that the model allows: those between main types with a positive upper
  # bound in the run (modules/10_land/landmatrix_dec18/presolve.gms); youngsecdf receives area only in
  # the presolve of module 35.
  .allowedTransitions <- function(mainOf, upper) {
    types <- names(mainOf)
    mainAllowed <- apply(as.array(upper) > 0, 3, any)
    allowed <- outer(types, types, "!=") & matrix(mainAllowed[outer(mainOf, mainOf, paste, sep = ".")], length(types))
    dimnames(allowed) <- list(types, types)
    allowed[, "youngsecdf"] <- FALSE
    return(allowed)
  }

  # Cluster x item matrix of one year of a magclass object.
  .yearMatrix <- function(x, year, clusters) {
    out <- as.array(x[, year, ])[, 1, , drop = FALSE]
    out <- matrix(out, nrow = dim(out)[1], dimnames = dimnames(out)[c(1, 3)])
    return(out[clusters, , drop = FALSE])
  }

  # Cluster transitions by year and stage (each a matrix clusters x "from.to"), in the order of the
  # model. peatLoss is the decline of drained peat by land type (clusters x years x types) or NULL.
  .stageTransitions <- function(dataLr, mainOf, allowed, clusters, statesLr, peatLoss = NULL) {
    types <- names(mainOf)
    forestryTypes <- types[mainOf == "forestry"]
    years <- getYears(dataLr[["transitions"]])

    recovery <- dimSums(dataLr[["recovery"]], dim = 3)
    maturation <- dimSums(dataLr[["maturation"]], dim = 3)
    # forestry abandonment of the previous year (0 in the first year)
    abandoned <- matrix(as.array(dataLr[["transitions"]])[clusters, , "forestry.other"], nrow = length(clusters))
    abandoned <- cbind(0, abandoned[, -length(years), drop = FALSE])

    stages <- list()
    for (i in seq_along(years)) {
      t <- years[i]
      mainTransitions <- .yearMatrix(dataLr[["transitions"]], t, clusters)

      # transitions of the optimisation between sub-types, with the area passing through forestry
      # as a stage of its own
      forestryArea <- matrix(statesLr[, forestryTypes, max(i - 1, 1)], nrow = length(clusters),
                             dimnames = list(clusters, forestryTypes))
      split <- .optimisationTransitions(dataLr, t, mainTransitions, forestryArea, mainOf, allowed, clusters)
      optimisation <- split[["optimisation"]]

      # SNV relocation (q29_land_snv_trans), taken from the transitions from cropland to secondary
      # forest and other land in proportion to them
      snvRelocation <- NULL
      if (!is.null(dataLr[["snvRelocation"]])) {
        toSnv <- c("crop.secdforest", "crop.othernat")
        available <- rowSums(optimisation[, toSnv, drop = FALSE])
        amount <- pmin(as.vector(dataLr[["snvRelocation"]][clusters, t, ]), available)
        snvRelocation <- optimisation[, toSnv, drop = FALSE] * (amount / ifelse(available > 0, available, 1))
        optimisation[, toSnv] <- pmax(optimisation[, toSnv] - snvRelocation, 0)
      }

      # land reduction on drained peat (q58_peatlandMan), taken from the transitions out of each
      # managed land type in proportion to them
      peatRewetting <- NULL
      if (!is.null(peatLoss)) {
        fromMain <- mainOf[sub("[.].*$", "", colnames(optimisation))]
        out <- fromMain %in% dimnames(peatLoss)[[3]] & fromMain != mainOf[sub("^.*[.]", "", colnames(optimisation))]
        available <- t(rowsum(t(optimisation[, out, drop = FALSE]), fromMain[out]))
        loss <- matrix(peatLoss[, t, colnames(available)], nrow(available))
        share <- pmin(available, loss) / ifelse(available > 0, available, 1)
        peatRewetting <- optimisation[, out, drop = FALSE] * share[, fromMain[out], drop = FALSE]
        optimisation[, out] <- optimisation[, out] - peatRewetting
        optimisation <- pmax(optimisation, 0)
        peatRewetting <- peatRewetting[, colSums(peatRewetting) > 0, drop = FALSE]
      }

      # recovery from the forestry abandonment of the previous step first, then from other othernat
      recoveryTotal <- as.vector(recovery[clusters, t, ])
      recoveryForestry <- pmin(recoveryTotal, abandoned[, i])

      stages[[t]] <- list(
        disturbance = .stageMatrix(dataLr[["disturbance"]][clusters, t, ], "primforest.secdforest", clusters),
        recoveryForestry = .stageMatrix(recoveryForestry, "othernat.youngsecdf", clusters),
        recovery = .stageMatrix(recoveryTotal - recoveryForestry, "othernat.youngsecdf", clusters),
        maturation = .stageMatrix(maturation[clusters, t, ], "youngsecdf.secdforest", clusters),
        forestryPassThrough = split[["passThrough"]]
      )
      stages[[t]][["snvRelocation"]] <- snvRelocation
      stages[[t]][["peatRewetting"]] <- peatRewetting
      stages[[t]][["optimisation"]] <- optimisation
    }
    return(stages)
  }

  # Stage matrix (clusters x one transition) of cluster values.
  .stageMatrix <- function(values, transition, clusters) {
    return(matrix(as.vector(values), ncol = 1, dimnames = list(clusters, transition)))
  }

  # Cluster inflow into a main type, from the transitions between main types.
  .mainInflow <- function(mainTransitions, type, mainTypes) {
    return(rowSums(mainTransitions[, paste(setdiff(mainTypes, type), type, sep = "."), drop = FALSE]))
  }

  # Transitions of the optimisation of year t between sub-types (clusters x "from.to"). The model
  # gives the outflow and inflow of each forestry type (expansion and reduction) and of each other
  # land type (reduction and harvest), but not their sources and destinations; the split fits them
  # to the transitions between main types. Area passing through forestry (where the transitions into
  # forestry exceed the expansion) is split over the forestry types in proportion to their area and
  # returned as a stage of its own (passThrough), because it enters forestry before it leaves.
  .optimisationTransitions <- function(dataLr, t, mainTransitions, forestryArea, mainOf, allowed, clusters) {
    types <- names(mainOf)
    mainTypes <- unique(mainOf)
    forestryTypes <- types[mainOf == "forestry"]
    otherTypes <- types[mainOf == "other"]

    # forestry types: expansion and reduction, plus the area passing through
    expansion <- .yearMatrix(dataLr[["forestryExpansion"]], t, clusters)[, forestryTypes, drop = FALSE]
    reduction <- .yearMatrix(dataLr[["forestryReduction"]], t, clusters)[, forestryTypes, drop = FALSE]
    reshuffle <- rowSums(expansion) - .mainInflow(mainTransitions, "forestry", mainTypes)
    passThrough <- pmax(-reshuffle, 0)

    forestrySum <- rowSums(forestryArea)
    areaShare <- forestryArea / ifelse(forestrySum > 0, forestrySum, 1)
    areaShare[forestrySum == 0, ] <- 1 / length(forestryTypes)
    expansion <- expansion + passThrough * areaShare
    reduction <- reduction + passThrough * areaShare

    # other land types: reduction and harvest; harvested youngsecdf regrows as othernat
    otherReduction <- .yearMatrix(dimSums(dataLr[["otherReduction"]], dim = "ac"), t, clusters)
    otherReduction <- otherReduction[, otherTypes, drop = FALSE]
    otherHarvest <- .yearMatrix(dimSums(dataLr[["otherHarvest"]], dim = "ac"), t, clusters)
    otherHarvest <- otherHarvest[, otherTypes, drop = FALSE]

    # outflow and inflow of each sub-type
    outflow <- matrix(0, length(clusters), length(types), dimnames = list(clusters, types))
    inflow <- outflow
    for (k in setdiff(types, c(forestryTypes, otherTypes))) {
      outflow[, k] <- rowSums(mainTransitions[, paste(k, setdiff(mainTypes, k), sep = "."), drop = FALSE])
      inflow[, k] <- .mainInflow(mainTransitions, k, mainTypes)
    }
    outflow[, forestryTypes] <- reduction
    inflow[, forestryTypes] <- expansion
    outflow[, "othernat"] <- otherReduction[, "othernat"] - otherHarvest[, "othernat"]
    outflow[, "youngsecdf"] <- otherReduction[, "youngsecdf"]
    inflow[, "othernat"] <- .mainInflow(mainTransitions, "other", mainTypes) + otherHarvest[, "youngsecdf"]
    inflow[, "youngsecdf"] <- 0
    outflow <- pmax(outflow, 0)

    # transitions between sub-types of one main type: changes of forestry type, harvest of youngsecdf
    internal <- matrix(0, length(clusters), length(mainTypes), dimnames = list(clusters, mainTypes))
    internal[, "forestry"] <- pmax(reshuffle, 0)
    internal[, "other"] <- otherHarvest[, "youngsecdf"]

    # a forestry type can expand and reduce in the same cluster and year; the model nets the
    # transitions through forestry, so the split needs transitions from a type to itself
    splitAllowed <- allowed
    splitAllowed[cbind(forestryTypes, forestryTypes)] <- TRUE
    # other land types whose whole reduction is harvest (regrowing as othernat)
    noExternal <- (otherReduction - otherHarvest) < 1e-9
    optimisation <- luscale::toolFitSubtypeTransitions(mainTransitions, outflow, inflow, mainOf, splitAllowed, internal,
                                                       noExternal)

    # the share of the inflow to the forestry types that passes through forestry
    from <- sub("[.].*$", "", colnames(optimisation))
    to <- sub("^.*[.]", "", colnames(optimisation))
    intoForestry <- to %in% forestryTypes & !(from %in% forestryTypes)
    passShare <- passThrough * areaShare / ifelse(expansion > 0, expansion, 1)
    passIn <- optimisation[, intoForestry, drop = FALSE] * passShare[, to[intoForestry], drop = FALSE]
    optimisation[, intoForestry] <- optimisation[, intoForestry] - passIn

    return(list(optimisation = optimisation, passThrough = passIn))
  }

  # Initial grid cell state by sub-type (grid cells x types): pasture and rangeland merged, forestry and other
  # land split into equal parts; toolAllocateTransitions then fits each sub-type to its cluster start total.
  .initialSubtypes <- function(initial, mainOf) {
    types <- names(mainOf)
    x <- as.array(initial)[, 1, ]
    x[x < 0] <- 0

    out <- matrix(0, nrow(x), length(types), dimnames = list(NULL, types))
    out[, c("crop", "primforest", "secdforest", "urban")] <- x[, c("crop", "primforest", "secdforest", "urban")]
    out[, "past"] <- x[, "past"] + x[, "range"]

    .split <- function(total, subTypes) {
      return(matrix(total / length(subTypes), length(total), length(subTypes)))
    }

    out[, types[mainOf == "forestry"]] <- .split(x[, "forestry"], types[mainOf == "forestry"])
    out[, types[mainOf == "other"]] <- .split(x[, "primother"] + x[, "secdother"], types[mainOf == "other"])
    return(out)
  }

  # =================================================================================================
  # Stage weights
  # =================================================================================================

  # Weight of the othernat outflow for recovery: the othernat inflow of the previous year (from
  # forestry for recoveryForestry) per othernat area. Othernat that regrows after harvest is not a
  # grid cell transition, so grid cells without inflow take the rest (weight 1e-6) in proportion to their
  # othernat area.
  .recoveryWeight <- function(stage, state, previousFlows, mainOf) {
    inflow <- 0
    for (flows in Filter(Negate(is.null), previousFlows)) {
      toOthernat <- sub("^.*[.]", "", colnames(flows)) == "othernat" &
        (stage != "recoveryForestry" | mainOf[sub("[.].*$", "", colnames(flows))] == "forestry")
      inflow <- inflow + rowSums(flows[, toOthernat, drop = FALSE])
    }
    area <- state[, "othernat"]
    weight <- ifelse(area > 0, inflow / area, 0) + 1e-6
    return(matrix(weight, ncol = 1, dimnames = list(NULL, "othernat")))
  }

  # Grid cell area of the land types with drained peat (grid cells x types crop, past, forestry).
  .drainedLand <- function(state, drained) {
    land <- vapply(drained, function(k) rowSums(state[, names(mainOf)[mainOf == k], drop = FALSE]),
                   numeric(nrow(state)))
    return(matrix(land, nrow(state), length(drained), dimnames = list(rownames(state), drained)))
  }

  # Weight of the outflow of cropland, pasture and forestry for land reduction on drained peat: the
  # drained peat under the land type per area of the type in the grid cell (peat: grid cells x categories);
  # grid cells without drained peat take the rest (weight 1e-6) where the peat grid cells lack land.
  .peatWeight <- function(state, peat, drained) {
    area <- .drainedLand(state, drained)
    weight <- ifelse(area > 0, peat[, drained, drop = FALSE] / area, 0) + 1e-6
    subTypes <- names(mainOf)[mainOf %in% drained]
    weight <- weight[, mainOf[subTypes], drop = FALSE]
    colnames(weight) <- subTypes
    return(weight)
  }

  # =================================================================================================
  # Outputs
  # =================================================================================================

  # Grid cell state by sub-type as magclass object with names "land.type" (type "total" for land types
  # without sub-types), the start labelled y1985.
  .landOutput <- function(state, mainOf, cells) {
    out <- aperm(state, c(1, 3, 2))
    dimnames(out) <- list(cells, c("y1985", dimnames(state)[[3]][-1]), unname(.typeLabels(mainOf)))
    return(.asCellMagpie(out, "land.type"))
  }

  # Gross grid cell transitions per year by stage and sub-type, as magclass object with the names
  # "stage.from.fromType.to.toType" in the order of the stages; transitions without flow in any grid cell
  # and year are left out.
  .transitionOutput <- function(transitions, mainOf, cells) {
    years <- names(transitions)
    label <- .typeLabels(mainOf)

    # flows of each year (grid cells x "stage.from.fromType.to.toType"), without columns of zeros
    byYear <- lapply(transitions, function(stages) {
      do.call(cbind, lapply(names(stages), function(stage) {
        flows <- stages[[stage]]
        colnames(flows) <- paste(stage, label[sub("[.].*$", "", colnames(flows))],
                                 label[sub("^.*[.]", "", colnames(flows))], sep = ".", recycle0 = TRUE)
        return(flows[, colSums(flows) > 0, drop = FALSE])
      }))
    })

    # items in the order of the stages, so that readers can replay them
    itemNames <- unique(unlist(lapply(byYear, colnames)))
    stageNames <- unique(unlist(lapply(transitions, names)))
    itemNames <- itemNames[order(match(sub("\\..*", "", itemNames), stageNames))]

    out <- array(0, dim = c(length(cells), length(years), length(itemNames)),
                 dimnames = list(cells, years, itemNames))
    for (t in years) {
      out[, t, colnames(byYear[[t]])] <- byYear[[t]]
    }

    return(.asCellMagpie(out, "stage.from.fromType.to.toType"))
  }

  # Labels "land.type" of the sub-types (type "total" for land types without sub-types).
  .typeLabels <- function(mainOf) {
    subType <- ifelse(names(mainOf) == mainOf, "total", names(mainOf))
    return(stats::setNames(paste(mainOf, subType, sep = "."), names(mainOf)))
  }

  # magclass object of an array grid cells x years x items, with the sets of the 0.5 degree outputs.
  .asCellMagpie <- function(x, itemSet) {
    out <- as.magpie(x, spatial = 1, temporal = 2)
    getSets(out, fulldim = FALSE) <- c("x.y.iso", "year", itemSet)
    return(out)
  }


  clusters <- getItems(dataLr[["transitions"]], 1)
  cells <- getCells(startHr)
  if (is.character(mapping)) mapping <- toolGetMapping(mapping, where = "local")
  cluster <- as.character(mapping[["cluster"]][match(cells, mapping[["cell"]])])
  stopifnot("mapping lacks grid cells of startHr" = !anyNA(cluster),
            "clusters of dataLr and mapping differ" = setequal(clusters, unique(cluster)),
            "snvTargetHr needs one value per grid cell of startHr" =
              is.null(snvTargetHr) | length(snvTargetHr) == length(cells))
  years <- getYears(dataLr[["transitions"]])

  # main land type of each sub-type of the run: forestry types of module 32 and other land types of
  # module 35; the other land types are their own sub-type
  forestryTypes <- getItems(dataLr[["forestry"]], "type32")
  otherTypes <- getItems(dataLr[["other"]], "othertype35")
  mainOf <- c(crop = "crop", past = "past",
              stats::setNames(rep("forestry", length(forestryTypes)), forestryTypes),
              primforest = "primforest", secdforest = "secdforest", urban = "urban",
              stats::setNames(rep("other", length(otherTypes)), otherTypes))
  types <- names(mainOf)
  allowed <- .allowedTransitions(mainOf, dataLr[["transitionsUpper"]])

  # cluster state at the end of each year by sub-type (clusters x types x years)
  otherLr <- dimSums(dataLr[["other"]], dim = "ac")
  statesLr <- array(0, dim = c(length(clusters), length(types), length(years)),
                    dimnames = list(clusters, types, years))
  for (t in years) {
    land <- .yearMatrix(dataLr[["land"]], t, clusters)
    state <- cbind(land[, c("crop", "past", "primforest", "secdforest", "urban"), drop = FALSE],
                   .yearMatrix(dataLr[["forestry"]], t, clusters),
                   .yearMatrix(otherLr, t, clusters))
    statesLr[, , t] <- state[, types]
  }

  # peat of the grid cells (grid cells x categories) as layer of toolAllocateTransitions: updated at the end of
  # each year with downscalePeatland on the land and the peatRewetting stage of the year
  peat <- NULL
  updatePeat <- NULL
  peatLoss <- NULL
  if (!is.null(peatlandLr)) {
    categories <- getNames(peatlandLr)
    drained <- intersect(categories, c("crop", "past", "forestry"))
    peatStart <- as.array(peatlandStartHr)[, 1, , drop = FALSE]
    common <- intersect(categories, dimnames(peatStart)[[3]])
    peat <- list(current = matrix(0, length(cells), length(categories), dimnames = list(cells, categories)),
                 out = list(), excess = list())
    peat[["current"]][, common] <- peatStart[cells, 1, common]

    # decline of drained peat by land type (clusters x years x types), 0 in the first year
    drainedLr <- as.array(peatlandLr)[clusters, years, drained, drop = FALSE]
    peatLoss <- pmax(drainedLr[, c(1, seq_along(years)[-length(years)]), , drop = FALSE] - drainedLr, 0)
    dimnames(peatLoss) <- dimnames(drainedLr)

    updatePeat <- function(year, state, flows, peat) {
      # land of the drained types and its reduction in the peatRewetting stage (grid cells x types), and
      # the peat before the year, as magclass objects of the year
      .asYear <- function(x) {
        .asCellMagpie(array(x, c(length(cells), 1, ncol(x)), list(cells, year, colnames(x))), "data")
      }
      rewetting <- flows[["peatRewetting"]]
      fromMain <- mainOf[sub("[.].*$", "", colnames(rewetting))]
      loss <- vapply(drained, function(k) rowSums(rewetting[, fromMain == k, drop = FALSE]), numeric(length(cells)))
      loss <- matrix(loss, length(cells), length(drained), dimnames = list(cells, drained))

      result <- downscalePeatland(.asYear(.drainedLand(state, drained)), peatlandLr[, year, ],
                                  .asYear(peat[["current"]]), mapping, .asYear(loss))
      peat[["current"]][] <- as.array(result[["peatlandHr"]])[, 1, ]
      peat[["out"]][[year]] <- result[["peatlandHr"]]
      peat[["excess"]][[year]] <- result[["excess"]]
      return(peat)
    }
  }

  # cluster transitions by stage, start state of the clusters and initial state of the grid cells
  stages <- .stageTransitions(dataLr, mainOf, allowed, clusters, statesLr, peatLoss)
  initialSub <- .initialSubtypes(startHr, mainOf)
  rownames(initialSub) <- cells

  urbanCells <- NULL
  if (!is.null(urbanHr)) {
    stopifnot("urbanHr lacks grid cells of startHr" = all(cells %in% getCells(urbanHr)))
    urbanCells <- matrix(as.array(urbanHr)[cells, , 1], length(cells), dimnames = list(cells, getYears(urbanHr)))
  }
  # main types in the constraints stand for their sub-types
  constraints <- lapply(constraints, function(x) {
    x[["types"]] <- unique(c(setdiff(x[["types"]], mainOf), names(mainOf)[mainOf %in% x[["types"]]]))
    return(x)
  })

  # compensation for the urban target: the rules of the land module, and no new inflow to
  # secondary forest (only harvest regrowth and restoration), youngsecdf and primary forest
  compensationAllowed <- allowed
  compensationAllowed[, c("secdforest", "youngsecdf", "primforest")] <- FALSE

  # cropland relocated for SNV in each grid cell before the year (updated at the first stage of a year)
  snvRelocated <- new.env()
  snvRelocated[["year"]] <- ""
  snvRelocated[["area"]] <- 0

  stageRules <- function(year, stage, state, previous, peat) {
    input <- list()
    if (snvRelocated[["year"]] != year) {
      snvRelocated[["year"]] <- year
      if (!is.null(previous[["snvRelocation"]])) {
        snvRelocated[["area"]] <- snvRelocated[["area"]] + rowSums(previous[["snvRelocation"]])
      }
    }

    # recovery on the othernat inflow of the previous year
    if (stage %in% c("recoveryForestry", "recovery") && !is.null(previous)) {
      input[["fromWeight"]] <- .recoveryWeight(stage, state,
                                               previous[c("snvRelocation", "peatRewetting", "optimisation")],
                                               mainOf)
    }

    # SNV relocation in proportion to the relocation target of the grid cells not relocated in earlier
    # years; the small weight 1e-6 lets grid cells without target supply the rest where the target grid
    # cells lack cropland
    if (stage == "snvRelocation" && !is.null(snvTargetHr)) {
      remaining <- pmax(snvTargetHr - snvRelocated[["area"]], 0)
      weight <- ifelse(state[, "crop"] > 0, remaining / state[, "crop"], 0) + 1e-6
      input[["fromWeight"]] <- matrix(weight, ncol = 1, dimnames = list(NULL, "crop"))
    }

    # land reduction on drained peat in proportion to the drained peat of the grid cells (peat of the
    # previous year)
    if (stage == "peatRewetting") {
      input[["fromWeight"]] <- .peatWeight(state, peat[["current"]], drained)
    }

    # urban land from the data (scaled by toolAllocateTransitions to the cluster urban land of the model)
    if (stage == "optimisation" && !is.null(urbanCells)) {
      input[["targetArea"]] <- list(urban = urbanCells[, year])
      input[["allowedTransitions"]] <- compensationAllowed
    }

    return(input)
  }

  # the cluster totals hold exactly (see .fitSourceArea), so the caps need not converge further
  # than 1e-6 Mha per cluster
  result <- luscale::toolAllocateTransitions(initialSub, mapping, statesLr, stages, stageRules, constraints,
                                             layer = peat, updateLayer = updatePeat, tol = 1e-6)

  peatOut <- NULL
  peatExcess <- NULL
  if (!is.null(peat)) {
    peatOut <- mbind(result[["layer"]][["out"]])
    peatExcess <- do.call(rbind, unname(result[["layer"]][["excess"]]))
  }

  return(list(landHr = .landOutput(result[["stateHr"]], mainOf, cells),
              transitionsHr = .transitionOutput(result[["transitionsHr"]], mainOf, cells),
              diagnostics = result[["diagnostics"]], excess = result[["excess"]],
              transitionChange = result[["transitionChange"]], deviationLr = result[["deviationLr"]],
              peatlandHr = peatOut, peatlandExcess = peatExcess))
}
