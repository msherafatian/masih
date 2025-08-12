# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
usethis::use_package("shiny")
usethis::use_package("shinydashboard")
usethis::use_package("Seurat")
usethis::use_package("ggplot2")
usethis::use_package("dplyr")
usethis::use_package("plotly")
usethis::use_package("DT")
usethis::use_package("viridis")
usethis::use_package("tidyr")
usethis::use_package("cancersea")
usethis::use_package("openxlsx")
usethis::use_package("corrplot")
usethis::use_package("slingshot")
usethis::use_package("SingleCellExperiment")
usethis::use_package("RColorBrewer")
usethis::use_package("shinyjs")

## Add modules ----
## Create a module infrastructure in R/
golem::add_module(name = "upload")
golem::add_module(name = "cluster_analysis")
golem::add_module(name = "cancersea")
golem::add_module(name = "cell_cycle")
golem::add_module(name = "explorer")
golem::add_module(name = "markers")
golem::add_module(name = "trajectory")
golem::add_module(name = "compare")
golem::add_module(name = "export")

## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct("seurat", module = "utils")
golem::add_fct("plotting", module = "utils")
golem::add_fct("cancersea")
golem::add_fct("markers")
golem::add_fct("trajectory")
golem::add_fct("comparative")
golem::add_fct("export")

## External resources
## Creates .js and .css files at inst/app/www
golem::add_css_file("custom")

## Add internal datasets ----
## If you have data in your package
# usethis::use_data_raw(name = "my_dataset", open = FALSE)

## Tests ----
## Add one line by test you want to create
usethis::use_test("app")

# Documentation

## Vignette ----
usethis::use_vignette("masih")
devtools::build_vignettes()

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
## 
## (You'll need GitHub there)
usethis::use_github()

# GitHub Actions
usethis::use_github_action()
# Chose one of the three
# See https://usethis.r-lib.org/reference/use_github_action.html
usethis::use_github_action_check_release()
usethis::use_github_action_check_standard()
usethis::use_github_action_check_full()
# Add action for PR
usethis::use_github_action_pr_commands()

# Travis CI
usethis::use_travis()
usethis::use_travis_badge()

# AppVeyor
usethis::use_appveyor()
usethis::use_appveyor_badge()

# Circle CI
usethis::use_circleci()
usethis::use_circleci_badge()

# Jenkins
usethis::use_jenkins()

# GitLab CI
usethis::use_gitlab_ci()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")