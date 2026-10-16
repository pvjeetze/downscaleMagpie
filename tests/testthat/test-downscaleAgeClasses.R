ac <- c("ac0", "ac5", "ac10", "acx")

# cluster data of one cluster A with values by age class (vectors by year, in the order of ac)
.clusterAc <- function(values, items = ac) {
  x <- magclass::new.magpie("A", c("y2000", "y2005"), items, fill = 0)
  for (year in names(values)) x[, year, ] <- values[[year]]
  return(x)
}

test_that("age classes follow the cell land and the age classes the model removes", {
  zero <- .clusterAc(list())
  dataLr <- list(
    secdforest = .clusterAc(list(y2000 = c(0, 1, 1, 0), y2005 = c(0, 0, 1, 0.5))),
    forestry = .clusterAc(list(), paste0("plant.", ac)),
    other = .clusterAc(list(y2005 = c(0.2, 0, 0, 0, 0, 0, 0, 0)),
                       c(paste0("othernat.", ac), paste0("youngsecdf.", ac))),
    secdforestDisturbance = zero, recovery = zero, maturation = zero,
    # the model deforests 0.5 Mha of the oldest secondary forest
    secdforestReduction = .clusterAc(list(y2005 = c(0, 0, 0, 0.5))), secdforestHarvest = zero,
    otherReduction = .clusterAc(list(), c(paste0("othernat.", ac), paste0("youngsecdf.", ac))),
    otherHarvest = .clusterAc(list(), c(paste0("othernat.", ac), paste0("youngsecdf.", ac))),
    forestryReduction = .clusterAc(list(), paste0("plant.", ac)),
    forestryReplant = magclass::new.magpie("A", c("y2000", "y2005"), "plant", fill = 0)
  )
  landItems <- c("crop.total", "forestry.plant", "secdforest.total", "other.othernat", "other.youngsecdf")
  land <- magclass::new.magpie(c("c1", "c2"), c("y2000", "y2005"), landItems, fill = 0)
  land[, , "secdforest.total"] <- c(1, 1, 0.5, 1)
  land[, "y2005", "other.othernat"] <- c(0, 0.2)
  transitions <- magclass::new.magpie(c("c1", "c2"), c("y2000", "y2005"),
                                      c("optimisation.secdforest.total.crop.total",
                                        "optimisation.crop.total.other.othernat"), fill = 0)
  transitions["c1", "y2005", "optimisation.secdforest.total.crop.total"] <- 0.5
  transitions["c2", "y2005", "optimisation.crop.total.other.othernat"] <- 0.2

  out <- downscaleAgeClasses(dataLr, land, transitions, data.frame(cell = c("c1", "c2"), cluster = "A"))
  # cell 1 loses the oldest age class (acx after the shift), cell 2 keeps all
  expect_equal(as.vector(out[["secdforest"]]["c1", "y2005", ]), c(0, 0, 0.5, 0), tolerance = 1e-8)
  expect_equal(as.vector(out[["secdforest"]]["c2", "y2005", ]), c(0, 0, 0.5, 0.5), tolerance = 1e-8)
  # new other land enters the youngest age class
  expect_equal(as.vector(out[["other"]]["c2", "y2005", ]), c(0.2, 0, 0, 0), tolerance = 1e-8)
  expect_lt(max(out[["diagnostics"]][["gridCellDeviation"]]), 1e-10)
  expect_lt(max(out[["diagnostics"]][["clusterDeviation"]]), 1e-8)
})
