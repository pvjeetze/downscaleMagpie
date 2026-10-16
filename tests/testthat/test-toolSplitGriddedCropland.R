# Two clusters: AAA with two grid cells and cropland, BBB with one grid cell whose cluster has no
# cropland in the model (it takes the shares of all clusters).
splitMapping <- data.frame(cell = c("1.1.AAA", "2.2.AAA", "3.3.BBB"), cluster = c("AAA.1", "AAA.1", "BBB.2"))

splitCase <- function() {
  landHr <- magclass::new.magpie(splitMapping[["cell"]], c("y1985", "y2000"),
                                 c("crop.total", "past.total", "forestry.aff", "other.othernat"), fill = 1)
  landHr[, , "crop.total"] <- rep(c(1, 2, 1), 2)
  cropSplitLr <- magclass::new.magpie(c("AAA.1", "BBB.2"), "y2000", c("crop_area", "crop_fallow", "crop_treecover"),
                                      fill = 0)
  cropSplitLr["AAA.1", , ] <- c(1.8, 0.9, 0.3)
  cropAreaLr <- magclass::new.magpie(c("AAA.1", "BBB.2"), "y2000",
                                     c("tece.rainfed", "tece.irrigated", "betr.rainfed", "betr.irrigated",
                                       "begr.rainfed", "begr.irrigated"), fill = 0)
  cropAreaLr["AAA.1", , ] <- c(1.2, 0.3, 0.2, 0, 0.1, 0)
  return(list(landHr = landHr, cropSplitLr = cropSplitLr, cropAreaLr = cropAreaLr))
}

test_that("gridded cropland is split by the cluster shares and the land split covers the land", {
  x <- splitCase()
  out <- toolSplitGriddedCropland(x[["landHr"]], x[["cropSplitLr"]], x[["cropAreaLr"]], splitMapping)
  split <- out[["landSplitHr"]]
  expect_equal(magclass::getNames(split),
               c("crop.kfo.rainfed", "crop.kfo.irrigated", "crop.kbe.rainfed", "crop.kbe.irrigated",
                 "crop.fallow.total", "crop.treecover.total", "past.total.total", "forestry.aff.total",
                 "other.othernat.total"))
  # grid cell 2: cropland 2, crop area 60 %, fallow 30 %, tree cover 10 %; crop area 2/3 rainfed
  # food crops, 1/6 irrigated food crops, 1/6 rainfed bioenergy crops
  expect_equal(as.vector(split["2.2.AAA", "y2000", ]), c(0.8, 0.2, 0.2, 0, 0.6, 0.2, 1, 1, 1))
  # the start takes the shares of the first year, the cluster without cropland the shares of all
  expect_equal(as.vector(split["2.2.AAA", "y1985", ]), as.vector(split["2.2.AAA", "y2000", ]))
  expect_equal(as.vector(split["3.3.BBB", "y2000", 1:6]), as.vector(split["2.2.AAA", "y2000", 1:6]) / 2)
  expect_equal(as.vector(out[["landHr"]][, , "crop"]), rep(c(1, 2, 1), 2))
  # crop area as share of the grid cell area, in the years of the model
  expect_equal(magclass::getYears(out[["cropAreaShareHr"]]), "y2000")
  expect_equal(as.vector(out[["cropAreaShareHr"]]["2.2.AAA", , "tece.rainfed"]), 0.8 / 5)
})
