# helpers ------------------------------------------------------------------------------------

# A MAgPIE run with one cluster of three cells and two years, consistent with its land accounting.
# y2000: disturbance 0.1, recovery 0.1, maturation 0.05; past to crop 0.3, othernat to crop 0.2,
# crop to plant 0.1, aff to othernat 0.05, harvest of youngsecdf 0.02 (regrowing as othernat).
# y2005: recovery 0.04 (all from the forestry abandonment of y2000); crop to past 0.1.
magpieCase <- function() {
  cells <- c("1.1.A", "2.2.A", "3.3.A")
  years <- c("y2000", "y2005")
  main <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
  .new <- function(names, sets = NULL) magclass::new.magpie("A.1", years, names, fill = 0, sets = sets)

  startHr <- magclass::new.magpie(cells, "y1995", c("crop", "past", "range", "forestry", "primforest",
                                                    "secdforest", "urban", "primother", "secdother"),
                                  fill = c(1, 0.2, 0.5, 0.5, 1, 0.5, 0, 0.5, 0, 0.2, 0, 0.3,
                                           1, 0.5, 0, 0.5, 1, 1, 0.1, 0, 0.2, 0.5, 0.3, 1, 0.2, 0.5, 0.5))

  dataLr <- list(transitions = .new(paste(rep(main, each = 7), rep(main, 7), sep = ".")),
                 land = .new(main),
                 forestry = .new(c("aff", "plant"), c("i", "j", "t", "type32")),
                 other = .new(c("othernat.ac0", "youngsecdf.ac0"), c("i", "j", "t", "othertype35", "ac")),
                 disturbance = .new("ac0"), recovery = .new("ac0"), maturation = .new("ac0"),
                 forestryExpansion = .new(c("aff", "plant")), forestryReduction = .new(c("aff", "plant")),
                 otherReduction = .new(c("othernat.ac0", "youngsecdf.ac0"), c("i", "j", "t", "othertype35", "ac")),
                 otherHarvest = .new(c("othernat.ac0", "youngsecdf.ac0"), c("i", "j", "t", "othertype35", "ac")))
  # the transitions that the land module fixes to 0
  dataLr[["transitionsUpper"]] <- .new(paste(rep(main, each = 7), rep(main, 7), sep = "."))
  dataLr[["transitionsUpper"]][, , ] <- Inf
  dataLr[["transitionsUpper"]][, , c(paste(setdiff(main, "primforest"), "primforest", sep = "."),
                                     "primforest.forestry", "primforest.other", "secdforest.other")] <- 0
  dataLr[["transitions"]][, "y2000", c("past.crop", "other.crop", "crop.forestry", "forestry.other")] <-
    c(0.3, 0.2, 0.1, 0.05)
  dataLr[["transitions"]][, "y2005", "crop.past"] <- 0.1
  dataLr[["land"]][, "y2000", ] <- c(2.1, 2.2, 0.55, 1.4, 2.65, 0.3, 2.8)
  dataLr[["land"]][, "y2005", ] <- c(2.0, 2.3, 0.55, 1.4, 2.65, 0.3, 2.8)
  dataLr[["forestry"]][, , ] <- c(0.15, 0.15, 0.4, 0.4)
  dataLr[["other"]][, , ] <- c(2.57, 2.53, 0.23, 0.27)
  dataLr[["disturbance"]][, "y2000", ] <- 0.1
  dataLr[["recovery"]][, , ] <- c(0.1, 0.04)
  dataLr[["maturation"]][, "y2000", ] <- 0.05
  dataLr[["forestryExpansion"]][, "y2000", "plant"] <- 0.1
  dataLr[["forestryReduction"]][, "y2000", "aff"] <- 0.05
  dataLr[["otherReduction"]][, "y2000", ] <- c(0.2, 0.02)
  dataLr[["otherHarvest"]][, "y2000", "youngsecdf"] <- 0.02

  return(list(dataLr = dataLr, startHr = startHr, mapping = data.frame(cell = cells, cluster = "A.1")))
}

# Cluster sum of a cell transition in a year.
clusterFlow <- function(transitionsHr, year, item) sum(transitionsHr[, year, item])

# tests --------------------------------------------------------------------------------------

test_that("cluster states, transitions by stage and sub-type, and output names follow the run", {
  case <- magpieCase()
  out <- downscaleLand(case[["dataLr"]], case[["startHr"]], case[["mapping"]])
  expect_lt(out[["deviationLr"]], 1e-6)
  expect_equal(magclass::getYears(out[["landHr"]]), c("y1985", "y2000", "y2005"))
  expect_equal(magclass::getNames(out[["landHr"]]),
               c("crop.total", "past.total", "forestry.aff", "forestry.plant", "primforest.total",
                 "secdforest.total", "urban.total", "other.othernat", "other.youngsecdf"))
  # cell area stays the same
  area <- magclass::dimSums(out[["landHr"]], dim = 3)
  expect_equal(as.vector(area[, "y2005", ]), as.vector(area[, "y1985", ]))

  tr <- out[["transitionsHr"]]
  expect_equal(magclass::getSets(tr, fulldim = FALSE)[3], "stage.from.fromType.to.toType")
  expect_equal(clusterFlow(tr, "y2000", "disturbance.primforest.total.secdforest.total"), 0.1)
  expect_equal(clusterFlow(tr, "y2000", "maturation.other.youngsecdf.secdforest.total"), 0.05)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.past.total.crop.total"), 0.3)
  # sub-types from the module flows: othernat to crop, crop to plant, aff to othernat, harvest of youngsecdf
  expect_equal(clusterFlow(tr, "y2000", "optimisation.other.othernat.crop.total"), 0.2, tolerance = 1e-6)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.crop.total.forestry.plant"), 0.1, tolerance = 1e-6)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.forestry.aff.other.othernat"), 0.05, tolerance = 1e-6)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.other.youngsecdf.other.othernat"), 0.02, tolerance = 1e-6)
  # the main transitions are the sum over stages and sub-types
  mainTr <- magclass::dimSums(tr, dim = c("stage", "fromType", "toType"))
  expect_equal(sum(mainTr[, "y2005", "crop.past"]), 0.1)
})

test_that("the initial state is split into sub-types by the cluster shares", {
  case <- magpieCase()
  land <- downscaleLand(case[["dataLr"]], case[["startHr"]], case[["mapping"]])[["landHr"]]
  # cluster start: aff 0.2 and plant 0.3 of forestry 0.5; youngsecdf 0.2 of other 3
  expect_equal(as.vector(land[, "y1985", "forestry.aff"]), c(0.2, 0, 0.3) * 0.4, tolerance = 1e-6)
  expect_equal(as.vector(land[, "y1985", "other.youngsecdf"]), c(0.7, 0.8, 1.5) * 0.2 / 3, tolerance = 1e-6)
})

test_that("recovery is placed on the othernat inflow from forestry of the previous year", {
  case <- magpieCase()
  tr <- downscaleLand(case[["dataLr"]], case[["startHr"]], case[["mapping"]])[["transitionsHr"]]
  recovery <- as.vector(tr[, "y2005", "recoveryForestry.other.othernat.other.youngsecdf"])
  expect_equal(sum(recovery), 0.04)
  # cell 2 has no forestry, so it had no abandonment
  expect_lt(recovery[2], 1e-4)
  inflow <- as.vector(tr[, "y2000", "optimisation.forestry.aff.other.othernat"])
  expect_true(all(inflow[c(1, 3)] > 0))
})

test_that("constraints on main types apply to their sub-types and take magclass bounds", {
  case <- magpieCase()
  upper <- magclass::new.magpie(case[["mapping"]][["cell"]], "y2000", "x", fill = c(0.2, 0, 1))
  out <- downscaleLand(case[["dataLr"]], case[["startHr"]], case[["mapping"]],
                       constraints = list(forest = list(types = "forestry", upper = upper)))
  # the caps hold to the convergence tolerance of 1e-6 Mha per cluster
  forestry <- magclass::dimSums(out[["landHr"]][, , "forestry"], dim = 3)
  expect_true(all(as.vector(forestry[, c("y2000", "y2005"), ]) <= c(0.2, 0, 1) + 1e-5))
  expect_lt(as.vector(forestry[2, "y2005", ]), 1e-5)
  expect_lt(out[["deviationLr"]], 1e-6)
})

test_that("land reduction on drained peat follows the drained peat, and peat equals downscalePeatland", {
  case <- magpieCase()
  cells <- case[["mapping"]][["cell"]]
  # the drained peat under cropland of the cluster declines by 0.1 in y2005
  peatlandLr <- magclass::new.magpie("A.1", c("y2000", "y2005"), c("crop", "intact", "rewetted"),
                                     fill = c(0.4, 0.3, 0.5, 0.5, 0, 0.1))
  peatlandStartHr <- magclass::new.magpie(cells, "y1995", c("crop", "intact"), fill = c(0.4, 0, 0, 0.5, 0, 0))
  out <- downscaleLand(case[["dataLr"]], case[["startHr"]], case[["mapping"]],
                       peatlandLr = peatlandLr, peatlandStartHr = peatlandStartHr)
  rewetting <- as.vector(out[["transitionsHr"]][, "y2005", "peatRewetting.crop.total.past.total"])
  expect_equal(rewetting, c(0.1, 0, 0), tolerance = 1e-5)
  expect_equal(as.vector(out[["peatlandHr"]][1, "y2005", ]), c(0.3, 0.5, 0.1), tolerance = 1e-5)
  expect_lt(out[["deviationLr"]], 1e-6)

  rewettingHr <- magclass::dimSums(out[["transitionsHr"]][, , "peatRewetting"],
                                   dim = c("stage", "fromType", "to", "toType"))
  standalone <- downscalePeatland(out[["landHr"]][, c("y2000", "y2005"), ], peatlandLr, peatlandStartHr,
                                  case[["mapping"]], rewettingHr)
  expect_equal(as.vector(standalone[["peatlandHr"]]), as.vector(out[["peatlandHr"]]))
})
