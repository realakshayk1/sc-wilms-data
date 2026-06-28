library(testthat)

test_that("repo root resolves", {
  source(file.path("phase1_mechanotypes", "utils.R"))
  root <- REPO_ROOT()
  expect_true(file.exists(file.path(root, "config", "paths.yaml")))
})

test_that("feature config defines 1-D score features", {
  source(file.path("phase1_mechanotypes", "utils.R"))
  cfg <- load_config()
  expect_true(length(cfg$features$features) >= 1)
  ids <- vapply(cfg$features$features, `[[`, "", "id")
  expect_true(all(nzchar(ids)))
})

test_that("demo seurat has cell_state and histology", {
  source(file.path("phase1_mechanotypes", "utils.R"))
  cfg <- load_config()
  demo <- make_demo_seurat(cfg)
  expect_true("cell_state" %in% colnames(demo$meta))
  expect_true("histology" %in% colnames(demo$meta))
  expect_equal(nrow(demo$counts), length(rownames(demo$counts)))
})
