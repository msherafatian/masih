#' Access files in the current app
#'
#' @param ... character vectors, specifying subdirectory and file(s)
#' within your package. The default, none, returns the root of the app.
#'
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "masih")
}

#' Read App Config
#'
#' @param value Value to retrieve from the config file
#' @param config GOLEM_CONFIG_ACTIVE value. If unset, R_CONFIG_ACTIVE.
#' If unset, "default".
#' @param use_parent Logical, scan the parent directory for config file.
#'
#' @noRd
get_golem_config <- function(
    value,
    config = Sys.getenv(
      "GOLEM_CONFIG_ACTIVE",
      Sys.getenv(
        "R_CONFIG_ACTIVE",
        "default"
      )
    ),
    use_parent = TRUE
) {
  config::get(
    value = value,
    config = config,
    file = app_sys("golem-config.yml"),
    use_parent = use_parent
  )
}

#' App configuration parameters
#' @noRd
app_config <- list(
  # Upload settings
  max_upload_size = 1000 * 1024^2,  # 1000 MB
  
  # Processing defaults
  default_nfeatures = 3000,
  default_dims = 30,
  default_resolution = 0.8,
  
  # QC defaults
  default_max_mt = 5,
  default_min_ncount = 200,
  default_max_ncount = 2500,
  default_min_features = 200,
  
  # Marker analysis defaults
  default_marker_logfc = 0.25,
  default_marker_minpct = 0.1,
  default_marker_topn = 10,
  default_marker_display_genes = 30,
  
  # Trajectory defaults
  default_trajectory_dims = 10,
  
  # Export defaults
  default_export_width = 7,
  default_export_height = 7,
  default_export_dpi = 300
)