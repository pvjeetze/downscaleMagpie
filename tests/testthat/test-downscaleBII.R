test_that("BII follows the cluster coefficients on the land shares, with BII 1 on primary other land", {
  skip_if_not_installed("magpie4")
  cells <- c("1.1.AAA", "2.2.AAA")
  types <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
  landHr <- magclass::new.magpie(cells, c("y1985", "y2000", "y2005"), types, fill = 0.1)
  landHr[, , "other"] <- 0.5
  # other land is primary in grid cell 1 and secondary in grid cell 2
  landIni <- magclass::new.magpie(cells, "y1995", c("primother", "secdother"), fill = c(0.5, 0, 0, 0.5))
  sideLayers <- magclass::new.magpie(cells, NULL, c("manpast", "rangeland", "forested", "nonforested"),
                                     fill = c(0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.7, 0.7))
  classes <- c("crop", "manpast", "rangeland", "forestry", "primforest", "secdforest", "urban", "other")
  biiLr <- magclass::new.magpie("AAA.1", c("y2000", "y2005"),
                                paste(rep(classes, each = 2), c("forested", "nonforested"), sep = "."), fill = 0.5)
  out <- downscaleBII(biiLr, landHr, sideLayers, landIni, data.frame(cell = cells, cluster = "AAA.1"))
  expect_equal(magclass::getYears(out), c("y2000", "y2005"))
  # grid cell 1: 0.5 of 1.1 Mha primary other land (BII 1), the rest BII 0.5
  expect_equal(as.vector(out["1.1.AAA", , ]), rep(0.5 / 1.1 + 0.5 * 0.6 / 1.1, 2), tolerance = 1e-8)
  expect_equal(as.vector(out["2.2.AAA", , ]), c(0.5, 0.5), tolerance = 1e-8)
})
