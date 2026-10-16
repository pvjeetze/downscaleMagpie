.transitionsMagpie <- function(flows, cells) {
  # flows: list by year of named lists by stage of matrices (cells x "from.to")
  .names <- function(x) unlist(lapply(names(x), function(s) paste(s, colnames(x[[s]]), sep = ".")))
  names <- unique(unlist(lapply(flows, .names)))
  arr <- array(0, c(length(cells), length(flows), length(names)), list(cells, names(flows), names))
  for (t in names(flows)) {
    for (s in names(flows[[t]])) arr[, t, paste(s, colnames(flows[[t]][[s]]), sep = ".")] <- flows[[t]][[s]]
  }
  x <- magclass::as.magpie(arr, spatial = 1, temporal = 2)
  magclass::getSets(x, fulldim = FALSE)[3] <- "stage.from.to"
  return(x)
}

test_that("outflow takes the newest area first or all area in proportion", {
  initial <- magclass::new.magpie("A.1", "y2000", c("crop", "other"), c(2, 1))
  flows <- list(y2010 = list(change = cbind(crop.other = 0.5)), y2020 = list(change = cbind(other.crop = 0.7)))
  transitions <- .transitionsMagpie(flows, "A.1")

  newest <- griddedLandOrigin(transitions, initial, "other")
  expect_equal(magclass::getNames(newest), c("other.other", "other.crop"))
  expect_equal(as.vector(newest[, , "other.crop"]), c(0, 0.5, 0))
  expect_equal(as.vector(newest[, , "other.other"]), c(1, 1, 0.8))

  proportional <- griddedLandOrigin(transitions, initial, "other", order = "proportional")
  expect_equal(as.vector(proportional[, "y2020", "other.crop"]), 0.5 * (1 - 0.7 / 1.5))
  expect_equal(as.vector(proportional[, "y2020", "other.other"]), 1 - 0.7 / 1.5)
})

test_that("area keeps its origin between tracked types, with the stages in order", {
  initial <- magclass::new.magpie("A.1", "y2000", c("crop", "othernat", "youngsecdf"), c(1, 1, 0))
  .year <- function(recovery, maturation, optimisation) {
    list(recovery = cbind(othernat.youngsecdf = recovery), maturation = cbind(youngsecdf.secdforest = maturation),
         optimisation = cbind(crop.othernat = optimisation))
  }
  flows <- list(y2010 = .year(0, 0, 0.4), y2020 = .year(0.3, 0.1, 0))
  out <- griddedLandOrigin(.transitionsMagpie(flows, "A.1"), initial, c("othernat", "youngsecdf"))
  expect_equal(as.vector(out[, "y2020", c("othernat.othernat", "othernat.crop")]), c(1, 0.1))
  expect_equal(as.vector(out[, "y2020", c("youngsecdf.othernat", "youngsecdf.crop")]), c(0, 0.2))
})

test_that("outflow beyond the start area passes through from the inflow of the stage", {
  initial <- magclass::new.magpie("A.1", "y2000", c("urban", "other", "past"), c(1, 0.1, 1))
  flows <- list(y2010 = list(change = cbind(urban.other = 0.3, other.past = 0.3)))
  out <- expect_silent(griddedLandOrigin(.transitionsMagpie(flows, "A.1"), initial, "other"))
  expect_equal(as.vector(out[, "y2010", c("other.other", "other.urban")]), c(0, 0.1))
})

test_that("the origins sum to the states of the tracked types", {
  set.seed(4)
  types <- c("crop", "past", "other")
  cells <- c("A.1", "A.2", "B.3")
  years <- c("y2000", "y2010", "y2020", "y2030")
  states <- array(runif(3 * 3 * 4), c(3, 3, 4), list(cells, types, years))
  for (i in 2:4) states[, , i] <- states[, , i] / rowSums(states[, , i]) * rowSums(states[, , 1])

  transitions <- .transitionsMagpie(luscale::toolMinimalTransitions(states), cells)
  initial <- array(states[, , 1], c(3, 1, 3), list(cells, "y2000", types))
  initial <- magclass::as.magpie(initial, spatial = 1, temporal = 2)
  for (order in c("newest", "proportional")) {
    out <- griddedLandOrigin(transitions, initial, c("past", "other"), order = order)
    total <- magclass::dimSums(out, dim = "origin")
    for (k in c("past", "other")) {
      expect_equal(unname(as.array(total[, , k])[, , 1]), unname(states[, k, ]), tolerance = 1e-12)
    }
  }
})
