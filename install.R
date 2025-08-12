install.packages("golem")

# In R console - Initial setup
install.packages("golem")
library(golem)

# Create the golem project
golem::create_golem(path = "/Users/user/Documents/scripts/R/Masih/", overwrite = TRUE)  # Using the name from your app

# Navigate to the project
setwd("/Users/user/Documents/scripts/R/Masih/")


# Set up the basic golem infrastructure
golem::set_golem_options()

# Add your package dependencies
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
