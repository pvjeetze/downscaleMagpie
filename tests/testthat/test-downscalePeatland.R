test_that("peat follows the cluster values and drained peat stays within the land of its type", {
  peatland <- magclass::new.magpie(c("A.1"), "y2000", c("crop", "intact", "rewetted"), fill = c(0.4, 1, 0.6))
  initial <- magclass::new.magpie(c("c1", "c2"), "y1995", c("crop", "intact"), fill = c(0.5, 0.5, 0.5, 0.5))
  # cell 1 keeps only 0.1 cropland, so its drained peat under cropland cannot exceed 0.1
  land <- magclass::new.magpie(c("c1", "c2"), "y2000", "crop", fill = c(0.1, 1))
  out <- downscalePeatland(land, peatland, initial, data.frame(cell = c("c1", "c2"), cluster = "A.1"))
  peat <- out[["peatlandHr"]]
  expect_equal(as.vector(magclass::dimSums(peat, dim = 1)), c(0.4, 1, 0.6), tolerance = 1e-6)
  expect_equal(as.vector(magclass::dimSums(peat, dim = 3)), c(1, 1), tolerance = 1e-5)
  expect_lte(as.vector(peat["c1", , "crop"]), 0.1 + 1e-5)
  expect_equal(out[["excess"]][["category"]], "crop")
})

test_that("downscalePeatland follows the land reduction on drained peat", {
  peatland <- magclass::new.magpie("A.1", c("y2000", "y2010"), c("crop", "rewetted"),
                                   fill = c(1, 0.5, 0, 0.5))
  initial <- magclass::new.magpie(c("c1", "c2"), "y1995", c("crop", "rewetted"), fill = c(0.5, 0.5, 0, 0))
  land <- magclass::new.magpie(c("c1", "c2"), c("y2000", "y2010"), c("crop.total", "forestry.aff"), fill = 1)
  rewetting <- magclass::new.magpie(c("c1", "c2"), c("y2000", "y2010"), "crop", fill = 0)
  # the cropland on peat of cell 2 is reduced, so its drained peat is rewetted
  rewetting["c2", "y2010", "crop"] <- 0.5
  mapping <- data.frame(cell = c("c2", "c1"), cluster = "A.1")
  out <- downscalePeatland(land, peatland, initial, mapping, rewetting)
  expect_equal(as.vector(out[["peatlandHr"]][, "y2010", "crop"]), c(0.5, 0), tolerance = 1e-6)
  expect_equal(as.vector(out[["peatlandHr"]][, "y2010", "rewetted"]), c(0, 0.5), tolerance = 1e-6)
  # without the land reduction both cells lose half of their drained peat
  out <- downscalePeatland(land, peatland, initial, mapping)
  expect_equal(as.vector(out[["peatlandHr"]][, "y2010", "crop"]), c(0.25, 0.25), tolerance = 1e-6)
})
