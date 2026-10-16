# Two clusters of two grid cells each, one year; land conservation of the model by cluster and of
# the data by grid cell.
landTypes <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
consvMapping <- data.frame(cell = c("1.1.AAA", "2.2.AAA", "3.3.BBB", "4.4.BBB"),
                           cluster = c("AAA.1", "AAA.1", "BBB.2", "BBB.2"))

consvCase <- function() {
  model <- magclass::new.magpie(c("AAA.1", "BBB.2"), "y2000",
                                paste(rep(landTypes, each = 2), c("protect", "restore"), sep = "."), fill = 0)
  data <- magclass::new.magpie(consvMapping[["cell"]], "y2000", landTypes, fill = 0,
                               sets = c("x.y.iso", "year", "data"))
  return(list(model = model, data = data))
}

test_that("without land, the data are scaled to the cluster values of the model", {
  x <- consvCase()
  # cluster AAA: the model protects 2 Mha natural vegetation (1 by restoration) and 0.4 Mha pasture
  x[["model"]]["AAA.1", , c("primforest.protect", "other.protect", "other.restore")] <- c(0.5, 0.5, 1)
  x[["model"]]["AAA.1", , "past.protect"] <- 0.4
  x[["data"]]["1.1.AAA", , c("primforest", "past")] <- c(0.5, 0.2)
  x[["data"]]["2.2.AAA", , "other"] <- 0.5
  # cluster BBB: data without protection in the model
  x[["data"]]["3.3.BBB", , "other"] <- 1
  out <- downscaleLandConservation(x[["model"]], x[["data"]], consvMapping)
  expect_equal(as.vector(out["1.1.AAA", , c("primforest", "past")]), c(1, 0.4))
  expect_equal(as.vector(out["2.2.AAA", , "other"]), 1)
  expect_equal(sum(out[c("3.3.BBB", "4.4.BBB"), , ]), 0)
})

test_that("with land, protection stays within the land of its type and meets the model where possible", {
  x <- consvCase()
  x[["model"]]["AAA.1", , "primforest.protect"] <- 1.5
  x[["model"]]["AAA.1", , "past.protect"] <- 0.05
  # grid cell 1 has less natural vegetation (1 Mha) than the data protect (1.5 Mha)
  x[["data"]]["1.1.AAA", , c("primforest", "past")] <- c(1.5, 0.1)
  x[["data"]]["2.2.AAA", , "other"] <- 0.2
  land <- magclass::new.magpie(consvMapping[["cell"]], "y2000", landTypes, fill = 0)
  land["1.1.AAA", , c("primforest", "secdforest", "past")] <- c(0.6, 0.4, 0.1)
  land["2.2.AAA", , c("other", "past")] <- c(2, 1)
  out <- downscaleLandConservation(x[["model"]], x[["data"]], consvMapping, landHr = land)
  # natural vegetation: 1 Mha in grid cell 1 (all of it, split as its land), the rest of 1.5 in
  # grid cell 2 (0.2 + 0.3 on its unprotected land); pasture: the model protects half of the data
  expect_equal(as.vector(out["1.1.AAA", , c("primforest", "secdforest")]), c(0.6, 0.4), tolerance = 1e-5)
  expect_equal(as.vector(out["2.2.AAA", , "other"]), 0.5, tolerance = 1e-5)
  expect_equal(as.vector(out["1.1.AAA", , "past"]), 0.05, tolerance = 1e-5)
  expect_equal(sum(out[c("3.3.BBB", "4.4.BBB"), , ]), 0)
})

test_that("with land, protection beyond the data fills grid cells with protection in the data first", {
  x <- consvCase()
  land <- magclass::new.magpie(consvMapping[["cell"]], "y2000", landTypes, fill = 0)
  # cluster BBB: other land of 1 Mha in both grid cells, data protection of 0.2 Mha in grid cell 3 only
  land[c("3.3.BBB", "4.4.BBB"), , "other"] <- 1
  x[["data"]]["3.3.BBB", , "other"] <- 0.2
  # the model protects 0.6 Mha: all in grid cell 3, which still has unprotected land
  x[["model"]]["BBB.2", , "other.protect"] <- 0.6
  out <- downscaleLandConservation(x[["model"]], x[["data"]], consvMapping, landHr = land)
  expect_equal(as.vector(out[c("3.3.BBB", "4.4.BBB"), , "other"]), c(0.6, 0), tolerance = 1e-4)
  # the model protects 1.5 Mha: grid cell 3 is full, the rest goes to grid cell 4
  x[["model"]]["BBB.2", , "other.protect"] <- 1.5
  out <- downscaleLandConservation(x[["model"]], x[["data"]], consvMapping, landHr = land)
  expect_equal(as.vector(out[c("3.3.BBB", "4.4.BBB"), , "other"]), c(1, 0.5), tolerance = 1e-4)
  # natural vegetation as one group: the protection of grid cell 4 follows its natural vegetation
  land["4.4.BBB", , c("secdforest", "other")] <- c(0.5, 0.5)
  out <- downscaleLandConservation(x[["model"]], x[["data"]], consvMapping, landHr = land)
  expect_equal(as.vector(out["4.4.BBB", , c("secdforest", "other")]), c(0.25, 0.25), tolerance = 1e-4)
})

test_that("inputs are matched by grid cell name, and mapping can be a file", {
  x <- consvCase()
  x[["model"]]["AAA.1", , "primforest.protect"] <- 1
  x[["data"]]["1.1.AAA", , "primforest"] <- 0.5
  x[["data"]]["2.2.AAA", , "primforest"] <- 1.5
  expected <- downscaleLandConservation(x[["model"]], x[["data"]], consvMapping)
  # the data in another grid cell order
  out <- downscaleLandConservation(x[["model"]], x[["data"]][rev(consvMapping[["cell"]]), , ], consvMapping)
  expect_equal(magclass::getCells(out), consvMapping[["cell"]])
  expect_equal(as.vector(out), as.vector(expected))
  mapFile <- tempfile(fileext = ".csv")
  utils::write.table(consvMapping, mapFile, sep = ";", row.names = FALSE, quote = FALSE)
  expect_equal(as.vector(downscaleLandConservation(x[["model"]], x[["data"]], mapFile)), as.vector(expected))
  expect_error(downscaleLandConservation(x[["model"]], x[["data"]][1:3, , ], consvMapping), "lacks grid cells")
})
