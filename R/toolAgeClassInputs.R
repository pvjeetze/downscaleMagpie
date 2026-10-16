#' Age-class inputs of a MAgPIE run
#'
#' Reads the inputs of the age-class disaggregation from the output folder of a run: the age-class
#' symbols of modules 32 and 35 from the fulldata.gdx, the clustermap, and the gridded land split and
#' transitions written by the land disaggregation (\code{\link{toolWriteDownscaledOutputs}}), which has to
#' run first.
#'
#' @inheritParams toolDownscalingInputs
#' @return Named list with \code{outputdir}, \code{mapping}, \code{dataLr} (the gdx symbols, as
#' \code{dataLr} of \code{\link{downscaleAgeClasses}}), \code{landSplitHr} and \code{transitionsHr}
#' (in the years of the age classes).
#' @author Patrick v. Jeetze
#' @seealso \code{\link{downscaleAgeClasses}}
#' @export
toolAgeClassInputs <- function(outputdir) {
  gdx <- file.path(outputdir, "fulldata.gdx")
  inFile <- list(landSplit = file.path(outputdir, "cell.land_split_0.5.mz"),
                 transitions = file.path(outputdir, "cell.land_transitions_0.5.mz"))
  for (f in inFile) {
    if (!file.exists(f)) stop(f, " not found; run scripts/output/extra/disaggregation.R first")
  }

  # Level of a variable of the run.
  .lvl <- function(symbol) {
    return(gdx2::readGDX(gdx, symbol, select = list(type = "level")))
  }

  # clustermap of the run
  mapFile <- Sys.glob(file.path(outputdir, "clustermap_*.rds"))
  if (length(mapFile) == 0) stop("Could not find map file!")
  if (length(mapFile) > 1) {
    warning("More than one map file found. First occurrence will be used!")
    mapFile <- mapFile[1]
  }

  dataLr <- list(
    secdforest = .lvl("ov35_secdforest"),
    forestry = .lvl("ov32_land"),
    other = .lvl("ov_land_other"),
    secdforestDisturbance = gdx2::readGDX(gdx, "p35_disturbance_loss_secdf"),
    forestryDisturbance = gdx2::readGDX(gdx, "p32_disturbance_loss_ftype32", react = "silent"),
    recovery = gdx2::readGDX(gdx, "p35_forest_recovery_area"),
    maturation = gdx2::readGDX(gdx, "p35_maturesecdf"),
    secdforestReduction = .lvl("ov35_secdforest_reduction"),
    secdforestHarvest = .lvl("ov35_hvarea_secdforest"),
    otherReduction = .lvl("ov35_other_reduction"),
    otherHarvest = .lvl("ov35_hvarea_other"),
    forestryReduction = .lvl("ov32_land_reduction"),
    forestryReplant = .lvl("ov32_land_replant")
  )
  years <- getYears(dataLr[["secdforest"]])

  return(list(outputdir = outputdir, mapping = readRDS(mapFile[1]), dataLr = dataLr,
              landSplitHr = read.magpie(inFile[["landSplit"]])[, years, ],
              transitionsHr = read.magpie(inFile[["transitions"]])[, years, ]))
}
