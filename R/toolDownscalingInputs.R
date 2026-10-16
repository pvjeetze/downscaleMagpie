#' Downscaling inputs of a MAgPIE run
#'
#' Reads the inputs of the MAgPIE land disaggregation from the output folder of a run: the
#' fulldata.gdx, the config and the 0.5 degree files. All settings of the run (cfg switches and gdx
#' values) are applied here, so that the downscaling functions take the result as data.
#'
#' @details The constraints of the land disaggregation follow the model: forest establishment
#' within the potential forest area (\code{q35_max_forest_establishment}, potential forest by
#' \code{c35_pot_forest_correction} and \code{c35_pot_forest_scenario}), cropland within the
#' available cropland (\code{c29_marginal_land}) reduced by the SNV share and its fader
#' (\code{q29_avl_cropland}), and protected natural vegetation and pasture (\code{vm_land.lo}) from
#' the conservation data selected by the switches of module 22 and fitted with
#' \code{\link{downscaleLandConservation}}. Where they cannot all hold within a cluster, the one with
#' the lowest priority is exceeded first.
#'
#' @param outputdir Output folder of a MAgPIE run, with fulldata.gdx, config.yml, the clustermap and
#' the 0.5 degree input files.
#' @param priorities Named vector with the priorities of the constraints \code{protection},
#' \code{potentialForest} and \code{crop} (see \code{\link[luscale]{toolAllocateTransitions}}).
#' @return Named list with
#' \item{outputdir, gdx}{the output folder and the path of its fulldata.gdx}
#' \item{mapping}{the clustermap of the run}
#' \item{landStartHr, landIniFile}{the land-use initialisation (1995) and its file}
#' \item{dataLr}{the gdx symbols of the land, as \code{dataLr} of \code{\link{downscaleLand}}}
#' \item{urbanHr, constraints, snvTargetHr, peatlandLr, peatlandStartHr}{inputs of
#' \code{\link{downscaleLand}} (\code{NULL} where the run has none)}
#' \item{peatlandCategories, peatlandFirstYear}{categories and first year of the peatland output}
#' \item{landConservationLr, conservationDataHr}{inputs of \code{\link{downscaleLandConservation}}
#' (\code{NULL} without WDPA data)}
#' \item{cropSplitLr, cropAreaLr}{cluster cropland by subcategory and crop area by crop and water,
#' for \code{\link{toolSplitGriddedCropland}}}
#' \item{biiLr, sideLayersHr}{BII coefficients of the clusters and the side layers, for
#' \code{\link{downscaleBII}}}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{downscaleLand}}, \code{\link{toolWriteDownscaledOutputs}}
#' @export
toolDownscalingInputs <- function(outputdir, # nolint: cyclocomp_linter
                                  priorities = c(protection = 1, potentialForest = 2, crop = 3)) {
  gdx <- file.path(outputdir, "fulldata.gdx")
  cfg <- gms::loadConfig(file.path(outputdir, "config.yml"))
  inFile <- list(
    landIni = file.path(outputdir, "avl_land_full_t_0.5.mz"),
    cropland = file.path(outputdir, "avl_cropland_0.5.mz"),
    snvTarget = file.path(outputdir, "SNVTargetCropland_0.5.mz"),
    urban = file.path(outputdir, "f34_urbanland_0.5.mz"),
    wdpa = file.path(outputdir, "wdpa_baseline_0.5.mz"),
    consvPrio = file.path(outputdir, "consv_prio_areas_0.5.mz"),
    sideLayers = file.path(outputdir, "luh2_side_layers_0.5.mz"),
    potentialForest = file.path(outputdir, "pot_forest_area_0.5.mz"),
    potentialForestUncorrected = file.path(outputdir, "pot_forest_area_uncorrected_0.5.mz"),
    peatland = file.path(outputdir, "f58_peatland_area_0.5.mz")
  )

  # Level of a variable of the run.
  .lvl <- function(symbol) {
    return(gdx2::readGDX(gdx, symbol, select = list(type = "level")))
  }

  # Potential forest area of the run: corrected for grassland ecoregions or not
  # (c35_pot_forest_correction), and its climate change scenario (c35_pot_forest_scenario) as in
  # modules/35_natveg/pot_forest_may24/input.gms.
  .potentialForest <- function() {
    if (identical(cfg[["gms"]][["c35_pot_forest_correction"]], "off")) {
      x <- read.magpie(inFile[["potentialForestUncorrected"]])
    } else {
      x <- read.magpie(inFile[["potentialForest"]])
    }
    years <- getYears(x, as.integer = TRUE)
    sourceYears <- years
    if (identical(cfg[["gms"]][["c35_pot_forest_scenario"]], "nocc")) {
      sourceYears[] <- 1995
    } else if (identical(cfg[["gms"]][["c35_pot_forest_scenario"]], "nocc_hist")) {
      fixYear <- as.numeric(gdx2::readGDX(gdx, "sm_fix_cc"))
      sourceYears[years > fixYear] <- fixYear
    }
    return(setYears(x[, paste0("y", sourceYears), ], getYears(x)))
  }

  # Land conservation of the model by cluster, and the 0.5 degree data selected by the switches of
  # module 22 (base protection from WDPA, future conservation from the conservation priority areas
  # with the fader, each for the policy countries and the others); NULL without WDPA data.
  .conservation <- function() {
    if (!file.exists(inFile[["wdpa"]])) {
      return(NULL)
    }
    wdpaHr <- read.magpie(inFile[["wdpa"]])
    landConservationLr <- gdx2::readGDX(gdx, "pm_land_conservation", react = "silent")
    conservationDataHr <- new.magpie(mapping[["cell"]], getYears(landConservationLr), getItems(wdpaHr, dim = 3.2),
                                     fill = 0, sets = c("x.y.iso", "year", "data"))
    consvIso <- gdx2::readGDX(gdx, "policy_countries22")
    consvIso <- consvIso[consvIso %in% getItems(wdpaHr, dim = 1.3)]
    if (length(consvIso) == 0) {
      warning("No countries selected in land conservation disaggregation. Results may be erroneous")
    }

    baseProtect <- cfg[["gms"]][["c22_base_protect"]]
    baseProtectNoselect <- cfg[["gms"]][["c22_base_protect_noselect"]]
    if (!all(c(baseProtect, baseProtectNoselect) %in% "none")) {
      if (baseProtectNoselect != "none") {
        conservationDataHr[, getYears(conservationDataHr), ] <-
          collapseDim(wdpaHr[, nyears(wdpaHr), baseProtectNoselect], dim = 3.1)
        conservationDataHr[, getYears(wdpaHr), ] <- collapseDim(wdpaHr[, , baseProtectNoselect], dim = 3.1)
      }
      if (baseProtect != "none") {
        conservationDataHr[consvIso, , ] <- collapseDim(wdpaHr[consvIso, nyears(wdpaHr), baseProtect], dim = 3.1)
      } else {
        conservationDataHr[consvIso, , ] <- 0
      }
    }

    consvPrioAll <- NULL
    if (file.exists(inFile[["consvPrio"]])) {
      consvPrioAll <- read.magpie(inFile[["consvPrio"]])
    } else {
      warning("Future land conservation used in MAgPIE run but high resolution ",
              "conservation priority data for disaggregation not found.")
    }
    protectScenario <- cfg[["gms"]][["c22_protect_scenario"]]
    protectScenarioNoselect <- cfg[["gms"]][["c22_protect_scenario_noselect"]]
    if (!all(c(protectScenario, protectScenarioNoselect) %in% "none") && !is.null(consvPrioAll)) {
      consvPrioHr <- new.magpie(cells_and_regions = mapping[["cell"]], names = getNames(consvPrioAll, dim = 2),
                                fill = 0, sets = c("x.y.iso", "year", "data"))
      if (protectScenarioNoselect != "none") {
        consvPrioHr <- collapseDim(consvPrioAll[, , protectScenarioNoselect], dim = 3.1)
      }
      if (protectScenario != "none") {
        consvPrioHr[consvIso, , ] <- collapseDim(consvPrioAll[consvIso, , protectScenario], dim = 3.1)
      } else {
        consvPrioHr[consvIso, , ] <- 0
      }
      # future conservation only pertains to natural vegetation
      consvPrioHr[, , c("crop", "past", "forestry", "urban")] <- 0
      consvFader <- gdx2::readGDX(gdx, "p22_conservation_fader", format = "first_found")
      conservationDataHr <- conservationDataHr + consvPrioHr * consvFader[, getYears(conservationDataHr), ]
    }
    return(list(landConservationLr = landConservationLr, conservationDataHr = conservationDataHr))
  }

  # clustermap of the run
  mapFile <- Sys.glob(file.path(outputdir, "clustermap_*.rds"))
  if (length(mapFile) == 0) stop("Could not find map file!")
  if (length(mapFile) > 1) {
    warning("More than one map file found. First occurrence will be used!")
    mapFile <- mapFile[1]
  }
  mapping <- readRDS(mapFile)

  # land-use initialisation (1995), available cropland, potential forest and urban land
  landIniHr <- read.magpie(inFile[["landIni"]])[, "y1995", ]
  croplandHr <- read.magpie(inFile[["cropland"]])[, , cfg[["gms"]][["c29_marginal_land"]]]
  potentialForestHr <- .potentialForest()
  urbanHr <- NULL
  if (cfg[["gms"]][["urban"]] == "exo_nov21") {
    urbanHr <- read.magpie(inFile[["urban"]])[, , cfg[["gms"]][["c34_urban_scenario"]]]
  }

  # land conservation, and the protected land of the model as the data fitted to the cluster values
  conservation <- .conservation()
  protectionHr <- NULL
  if (!is.null(conservation)) {
    protectionHr <- downscaleLandConservation(conservation[["landConservationLr"]],
                                              conservation[["conservationDataHr"]], mapping)
  }

  # semi-natural vegetation (SNV): country shares in cropland and their fader
  iso <- gdx2::readGDX(gdx, "iso")
  snvShareIso <- stats::setNames(rep(as.numeric(gdx2::readGDX(gdx, "s29_snv_shr_noselect")), length(iso)), iso)
  snvShareIso[gdx2::readGDX(gdx, "policy_countries29")] <- as.numeric(gdx2::readGDX(gdx, "s29_snv_shr"))
  snvShareHr <- unname(snvShareIso[as.character(mapping[["country"]])])
  snvFader <- gdx2::readGDX(gdx, "i29_snv_scenario_fader")
  snvFader <- stats::setNames(as.vector(snvFader), getYears(snvFader))

  # SNV cropland relocation of the model (NULL in other realisations of module 29) and the cropland
  # that requires relocation (Copernicus land cover), interpolated between the 20 % and 50 % targets
  # for the SNV share of the run as i29_snv_relocation_target
  snvRelocationLr <- gdx2::readGDX(gdx, "p29_snv_relocation", react = "silent")
  snvTargetHr <- NULL
  if (!is.null(snvRelocationLr) && sum(snvRelocationLr) > 0) {
    if (!file.exists(inFile[["snvTarget"]])) {
      stop("SNV relocation in the run, but ", inFile[["snvTarget"]], " not found")
    }
    snvTargetData <- read.magpie(inFile[["snvTarget"]])
    stopifnot("cell order of the 0.5 degree inputs differs from the cluster map" =
                identical(getCells(snvTargetData), as.character(mapping[["cell"]])))
    target20 <- as.vector(snvTargetData[, , "SNV20TargetCropland"])
    target50 <- as.vector(snvTargetData[, , "SNV50TargetCropland"])
    x1 <- as.numeric(gdx2::readGDX(gdx, "s29_snv_relocation_data_x1"))
    x2 <- as.numeric(gdx2::readGDX(gdx, "s29_snv_relocation_data_x2"))
    snvTargetHr <- pmin(snvShareHr / x1, 1) * target20 + pmax(snvShareHr - x1, 0) / (x2 - x1) * (target50 - target20)
  }

  # peatland: cluster values of the run and initial grid cell values
  peatLr <- NULL
  peatIniHr <- NULL
  peatCategories <- NULL
  if (cfg[["gms"]][["peatland"]] == "v2") {
    peatLr <- magpie4::PeatlandArea(gdx, level = "cell", sum = FALSE)
    peatIniHr <- read.magpie(inFile[["peatland"]])
    peatCategories <- c(getNames(peatIniHr), "rewetted", "unused")
  }

  # all 0.5 degree inputs in the grid cell order of the cluster map
  for (x in list(landIniHr, urbanHr, croplandHr, potentialForestHr, protectionHr, peatIniHr)) {
    if (!is.null(x) && !identical(getCells(x), as.character(mapping[["cell"]]))) {
      stop("cell order of the 0.5 degree inputs differs from the cluster map")
    }
  }

  # constraints: forest establishment within the potential forest area, cropland within the
  # available cropland reduced by the SNV share, and protected natural vegetation (as a whole,
  # because the land types of the protected-area data need not match those of the initial map) and
  # pasture
  cropCapHr <- as.vector(croplandHr) * (1 - outer(snvShareHr, snvFader))
  constraints <- list(
    potentialForest = list(types = c("forestry", "primforest", "secdforest", "youngsecdf"),
                           upper = potentialForestHr, priority = priorities[["potentialForest"]]),
    crop = list(types = "crop", upper = cropCapHr, priority = priorities[["crop"]])
  )
  if (!is.null(protectionHr)) {
    constraints[["natveg"]] <- list(types = c("primforest", "secdforest", "other"),
                                    priority = priorities[["protection"]],
                                    lower = dimSums(protectionHr[, , c("primforest", "secdforest", "other")], dim = 3))
    constraints[["past"]] <- list(types = "past", lower = protectionHr[, , "past"],
                                  priority = priorities[["protection"]])
  }

  # land variables and presolve shifts of the run
  dataLr <- list(
    transitions = .lvl("ov_lu_transitions"),
    transitionsUpper = gdx2::readGDX(gdx, "ov_lu_transitions", select = list(type = "upper")),
    land = .lvl("ov_land"),
    forestry = .lvl("ov_land_forestry"),
    other = .lvl("ov_land_other"),
    disturbance = gdx2::readGDX(gdx, "p35_disturbance_loss_primf"),
    recovery = gdx2::readGDX(gdx, "p35_forest_recovery_area"),
    maturation = gdx2::readGDX(gdx, "p35_maturesecdf"),
    forestryExpansion = .lvl("ov_landexpansion_forestry"),
    forestryReduction = .lvl("ov_landreduction_forestry"),
    otherReduction = .lvl("ov35_other_reduction"),
    otherHarvest = .lvl("ov35_hvarea_other"),
    snvRelocation = snvRelocationLr
  )

  # cropland by subcategory, for the split of the gridded cropland
  cropSplitLr <- magpie4::land(gdx, level = "cell", subcategories = "crop")
  cropSplitLr <- cropSplitLr[, , c("crop_area", "crop_fallow", "crop_treecover")]

  return(list(
    outputdir = outputdir, gdx = gdx, mapping = mapping, landStartHr = landIniHr, landIniFile = inFile[["landIni"]],
    dataLr = dataLr, urbanHr = urbanHr, constraints = constraints, snvTargetHr = snvTargetHr,
    peatlandLr = peatLr, peatlandStartHr = peatIniHr, peatlandCategories = peatCategories,
    peatlandFirstYear = cfg[["gms"]][["s58_fix_peatland"]],
    landConservationLr = conservation[["landConservationLr"]],
    conservationDataHr = conservation[["conservationDataHr"]],
    cropSplitLr = cropSplitLr,
    cropAreaLr = magpie4::croparea(gdx, level = "cell", products = "kcr", product_aggr = FALSE, water_aggr = FALSE),
    biiLr = magpie4::BII(gdx, file = NULL, level = "cell", mode = "auto", landClass = "all",
                         bii_coeff = NULL, side_layers = NULL),
    sideLayersHr = read.magpie(inFile[["sideLayers"]])
  ))
}
