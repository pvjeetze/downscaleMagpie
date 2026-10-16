test_that("products are written with the file names of the MAgPIE outputs", {
  outputdir <- withr::local_tempdir()
  cells <- c("1.1.AAA", "2.2.AAA")
  landSplitHr <- magclass::new.magpie(cells, "y2000", c("crop.kfo.rainfed", "past.total.total"), fill = c(1, 2, 3, 4))
  cropland <- list(landSplitHr = landSplitHr, landHr = magclass::dimSums(landSplitHr, dim = c(3.2, 3.3)),
                   cropAreaShareHr = magclass::new.magpie(cells, "y2000", "tece.rainfed", fill = 0.1))
  conservation <- magclass::new.magpie(cells, "y2000", "primforest", fill = 0.5)
  stock <- magclass::new.magpie(cells, "y2000", c("ac0", "ac5"), fill = c(1, 0, 1, 1))
  ageClasses <- list(secdforest = stock, forestry = stock, other = stock,
                     diagnostics = data.frame(gridCellDeviation = 0, clusterDeviation = 0))
  expect_message(toolWriteDownscaledOutputs(list(outputdir = outputdir), cropland = cropland,
                                            conservation = conservation,
                                            ageClasses = ageClasses), "cell.land_split")
  files <- c("cell.croparea_0.5_share.mz", "cell.land_split_0.5.mz", "cell.land_split_0.5_share.mz",
             "cell.land_0.5.mz", "cell.land_0.5_share.mz", "cell.conservation_land_0.5.mz",
             "cell.ageclass_secdforest_0.5_share.mz", "cell.ageclass_diagnostics.csv")
  expect_true(all(file.exists(file.path(outputdir, files))))
  land <- magclass::read.magpie(file.path(outputdir, "cell.land_0.5.mz"))
  expect_equal(magclass::getSets(land, fulldim = FALSE)[1], "x.y.iso")
  expect_equal(as.vector(land), c(1, 2, 3, 4))
  share <- magclass::read.magpie(file.path(outputdir, "cell.ageclass_secdforest_0.5_share.mz"))
  # grid cell 1 holds both age classes equally, grid cell 2 only ac5
  expect_equal(as.vector(share), c(0.5, 0, 0.5, 1))
})
