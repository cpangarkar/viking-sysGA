rm(list=ls())
graphics.off()

library(tidyverse)
library(ggh4x)
library(gt)
library(pdftools)
library(mcr)
library(patchwork)
library(VCA)
library(EnvStats)
library(ggpmisc)
library(readxl)
library(magrittr)
library(sysfonts)
library(showtext)
library(jsonlite)
library(tidyverse)
library(httr2)
library(gt)
library(patchwork)
library(tidyxl)
library(qs2)
library(ggrepel)


path = "C:/Users/cpang/analysis/viking-sysGA/"
setwd(path)
source("./scripts/config.R")
source("./scripts/library_functions.R")
source("./scripts/library.R")

## Configs
  # config for MAR - download features table one time
  # fdf = run_id |> download_features()
  mar = read_csv("./config/mar.csv")

  # read operator log and extract run_id URLS
  operator_log = read_csv("./config/operator-log-index.csv")

  # This creates a tibble named ares with run_name to run_id mapping

  # read primary wavelengths <-> assay mapping
  assay_lambdas = read_csv("./config/assay_wavelength_mapping.csv")

  # Set active server here
  ATLANTIS_SERVER = ATLANTIS_SERVERS$develop

## run and assay of interest
run_input = qs_read("./data/all_run_inputs.qs2")
assayname = "Ca"

# download CSVs for each run_id of interest: don't do it every time
df = 
  run_input |> 
  map(~{
    result = tryCatch({
      operator_log |> 
        filter(run_name == .x) |> 
        pull(run_id) |> 
        download_relOD_csv() |> 
        pull(ro_presigned_url) |> 
        read_csv(show_col_types = FALSE) |>
        mutate(run_name = .x)
    }, error = function(e) {
      message("Skipping ", .x, ": ", conditionMessage(e))
      NULL
    })
    return(result)
  }) |> 
  compact() |> 
  list_rbind()
  
qs_save(df, "./data/df.qs2")

df = 
  df |> 
  group_by(run_name) |> 
  group_modify(~{
    .x |> 
    mutate(tassay = (timestamp - first(timestamp)) |> as.numeric()) |> 
    select(-metadata, -timestamp) |> 
    pivot_longer(col = -c(tassay, read_path, read_id)) |> 
    magrittr::set_colnames(c("read_path", "read_id", "tassay", "well", "relOD")) |> 
    left_join(mar, by="well") |> 
    mutate(read_path = read_path |> str_split_i("_",4) |> tolower()) -> .x2
    return(.x2)
  })

df =
  df |>
  group_by(run_name, well, read_path, read_id) |>
  group_map(~ {
    .x = .x |> distinct()
    if (.y$read_id == "endpoint") {
      if (nrow(.x) != 2) {
        .x = .x |> mutate(read_id = c("2r-endpoint"))  
      }
      .x = .x |> mutate(read_id = c("1r-endpoint", "2r-endpoint"))
    } else {
      .x = .x |> mutate(read_id = .y$read_id)
    }
    return(.x)
  }, .keep=T) |>
  list_rbind()

qs_save(df, "./data/df_long.qs2")
