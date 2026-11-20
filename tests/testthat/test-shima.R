test_that("the shima namespace loads", {
  expect_true("shima" %in% loadedNamespaces())
})

test_that("run_app remains exported with the original interface", {
  exports <- getNamespaceExports("shima")

  expect_true("run_app" %in% exports)
  expect_true(is.function(getExportedValue("shima", "run_app")))

  run_app <- getExportedValue("shima", "run_app")
  expect_identical(
    names(formals(run_app)),
    c("onStart", "options", "enableBookmarking", "uiPattern", "...")
  )
  expect_identical(formals(run_app)$onStart, quote(NULL))
  expect_identical(formals(run_app)$options, quote(list()))
  expect_identical(formals(run_app)$enableBookmarking, quote(NULL))
  expect_identical(formals(run_app)$uiPattern, quote("/"))
})

test_that("bundled application resources are discoverable", {
  root <- system.file(package = "shima")
  config <- system.file("golem-config.yml", package = "shima")
  www <- system.file("app/www", package = "shima")

  expect_true(nzchar(root))
  expect_true(nzchar(config))
  expect_true(file.exists(config))
  expect_true(nzchar(www))
  expect_true(dir.exists(www))
})

test_that("the golem configuration loads", {
  get_golem_config <- getFromNamespace("get_golem_config", "shima")

  expect_identical(get_golem_config("golem_name"), "Shima")
  expect_identical(get_golem_config("golem_version"), "0.0.0.9000")
})

test_that("basic Shiny application objects can be constructed", {
  app_ui <- getFromNamespace("app_ui", "shima")
  app_server <- getFromNamespace("app_server", "shima")

  ui_object <- app_ui(list())
  app_object <- shiny::shinyApp(ui = app_ui, server = app_server)

  expect_false(is.null(ui_object))
  expect_s3_class(app_object, "shiny.appobj")
})
