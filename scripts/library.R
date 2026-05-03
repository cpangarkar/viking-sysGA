# =============================================================================
# ATLANTIS API LIBRARY
# =============================================================================
# Servers:
#   atlantis.internal.vital.company  (default)
#   atlantis.develop
#   atlantis.staging
#   atlantis.research
#   atlantis.production
# =============================================================================

library(httr2)
library(jsonlite)
library(dplyr)
library(purrr)


# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------

ATLANTIS_SERVERS = list(
  internal    = "https://atlantis.internal.vital.company",
  develop     = "https://atlantis.develop.vital.company",
  staging     = "https://atlantis.staging.vital.company",
  research    = "https://atlantis.research.vital.company",
  production  = "https://atlantis.production.vital.company"
)

atlantis_url = function(endpoint) {
  paste0(ATLANTIS_SERVER, endpoint)
}


# -----------------------------------------------------------------------------
# HELPERS (private)
# -----------------------------------------------------------------------------

# Base POST request to Atlantis API
.atlantis_post = function(url, body, filepath, verbosity = 0) {
  url |>
    request() |>
    req_headers(MyHeader = "Content-Type: application/json") |>
    req_body_json(data = body) |>
    req_perform(verbosity = verbosity, path = filepath)
}

# Handle paginated downloads generically
# .fetch_fn must accept (url, ...) and return list(num_pages, df)
.download_paginated = function(base_url, fetch_fn, ...) {
  args = list(...)
  ans  = do.call(fetch_fn, c(list(url = base_url), args))
  n    = ans$num_pages

  if (is.na(n) || n < 1) return(ans$df)

  message(sprintf("Download has %d page(s). Fetching all...", n))

  pages_df =
    seq(2, n) |>
    map(~ {
      url_iter = paste0(base_url, "?page=", .x)
      message("  page ", .x, " / ", n)
      do.call(fetch_fn, c(list(url = url_iter), args))$df
    }) |>
    list_rbind()

  bind_rows(ans$df, pages_df)
}


# -----------------------------------------------------------------------------
# STUDY
# -----------------------------------------------------------------------------

#' Get all run IDs belonging to a study
get_runs_in_study = function(string_study_id, server = ATLANTIS_SERVER) {
  filepath = "./data/output_runs_in_study.json"
  url      = paste0(server, "/api/run/filter/")

  .atlantis_post(url, list(study_id = string_study_id), filepath, verbosity = 3)

  data = fromJSON(filepath)
  data$data$results$id |> as_tibble() |> rename(run_id = value)
}


# -----------------------------------------------------------------------------
# RUN RESULTS — Relative OD
# -----------------------------------------------------------------------------

#' Internal single-page fetch for relative OD
.fetch_relativeOD = function(url, string_run_id, filepath) {
  .atlantis_post(
    url,
    list(run_id               = string_run_id,
         reading_result_type  = "READING_RESULT_TYPE_RELATIVE_OD"),
    filepath
  )
  data = fromJSON(filepath)
  list(num_pages = data$data$num_pages,
       df        = data$data$results |> as_tibble())
}

#' Download all pages of relative OD data for a run
download_relativeOD = function(string_run_id, server = ATLANTIS_SERVER) {
  filepath = "./data/output_relative_od.json"
  base_url = paste0(server, "/api/readingresult/filter/")

  .download_paginated(base_url, .fetch_relativeOD,
                      string_run_id = string_run_id,
                      filepath      = filepath)
}


# -----------------------------------------------------------------------------
# RUN RESULTS — Temperatures
# -----------------------------------------------------------------------------

#' Internal single-page fetch for temperatures
.fetch_temperatures = function(url, string_run_id, filepath) {
  result = try(
    .atlantis_post(url, list(run_id = string_run_id), filepath),
    silent = TRUE
  )

  if (inherits(result, "try-error")) {
    return(list(num_pages = NA, df = data.frame()))
  }

  data = fromJSON(filepath)
  list(num_pages = data$data$num_pages,
       df        = data$data$results |> as_tibble())
}

#' Download all pages of temperature data for a run
download_temperatures = function(string_run_id, server = ATLANTIS_SERVER) {
  filepath = "./data/output_temps.json"
  base_url = paste0(server, "/api/temperature/filter/")

  ans = .fetch_temperatures(base_url, string_run_id, filepath)

  if (nrow(ans$df) == 0) {
    message("No temperature data found for run: ", string_run_id)
    return(ans$df)
  }

  .download_paginated(base_url, .fetch_temperatures,
                      string_run_id = string_run_id,
                      filepath      = filepath)
}


# -----------------------------------------------------------------------------
# FEATURES
# -----------------------------------------------------------------------------

#' Trigger server-side feature generation for a run
generate_features = function(string_run_id,
                             output_filepath = "./data/dummy.json",
                             server          = ATLANTIS_SERVER) {
  url = paste0(server, "/custom/postprocessing/generate_features")
  .atlantis_post(url, list(id = string_run_id), output_filepath, verbosity = 3)
}

#' Download and parse features for a run, extracting relative OD signal
download_features = function(string_run_id,
                             output_filepath = "./data/output_features.json",
                             server          = ATLANTIS_SERVER) {
  url = paste0(server, "/api/feature/filter/")
  .atlantis_post(url, list(run_id = string_run_id), output_filepath, verbosity = 3)

  data = fromJSON(output_filepath)
  df   = data$data$results$metadata |> as_tibble()

  # Extract relative OD signal — one value per row across spread columns
  ddf    = df$input_cols |> select(starts_with("rela"))
  signal = map_dbl(seq(nrow(ddf)), function(i) {
    vals = ddf[i, ]
    if (any(!is.na(vals))) vals[[which(!is.na(vals))[1]]] else NA_real_
  })

  bind_cols(
    assay  = df$context$assay_code,
    well   = df$context$well_master_code,
    signal = signal
  ) |>
    mutate(run_id = string_run_id)
}

# -----------------------------------------------------------------------------
# Artifacts
# -----------------------------------------------------------------------------

download_relOD_csv = function(object_id,
  artifact_type = "ARTIFACT_TYPE_RELATIVE_OD",
  server        = ATLANTIS_SERVER) {
  filepath = "./data/output_artifacts.json"

  tryCatch({
    .atlantis_post(
    paste0(server, "/api/artifact/filter/"),
    list(object_id     = object_id,
    artifact_type = artifact_type),
    filepath
    )

    data = fromJSON(filepath)
    data$data$results |> as_tibble()

    }, error = function(e) {
    message("Skipping ", object_id, ": ", conditionMessage(e))
    return(NULL)
  })
}

download_run_results = function(run_id, server = ATLANTIS_SERVER) {
  filepath = "./data/output_run_results.json"

  .atlantis_post(
    paste0(server, "/api/runresult/filter/"),
    list(run_id = run_id),
    filepath
  )

  data = fromJSON(filepath)
  data$data$results |> as_tibble()
}