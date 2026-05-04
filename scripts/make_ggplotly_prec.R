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
library(plotly)
library(htmlwidgets)

path = "C:/Users/cpang/analysis/viking-sysGA/"
setwd(path)
source("./scripts/config.R")
source("./scripts/library_functions.R")
source("./scripts/library.R")



## Config---- 
  ## run and assay of interest - copy a column of run_names from googlesheet (operator log)
  run_input = qs_read("./data/all_run_inputs.qs2")


  # read primary wavelengths <-> assay mapping
  assay_lambdas = read_csv("./config/assay_wavelength_mapping.csv")

## Read previously saved data
assay_df =
  "./data/precision_assay_df.qs" |> 
  qs_read() |> 
  rename(run_name = disc) |> 
  filter(sample |> tolower()|> str_detect("control"))

df = qs_read("./data/df_long.qs2")

df = 
  df |> 
  left_join(assay_df |> select(run_name, sample), by="run_name")

make_ggplotly_method_comparison = function(.ddf, label_df, ylims, assayname, sampless) {
  p = 
    .ddf |>
    arrange(run_name, tassay) |>
    ggplot(aes(x = tassay, y = relOD, group = run_name,
      text = paste0("Run: ", run_name, "<br>Recovery: ", recovery))) +
    geom_line(linewidth = 0.5, col = "grey80") +
    geom_text(
      data  = label_df,
      aes(label = run_name),
      hjust = -0.1,
      size  = 2.5,
      col   = "grey60"
    ) +
    # geom_line(
    #   data = .ddf |> filter(run_name == highlight_run) |> arrange(read_id, tassay),
    #   aes(col = read_id, group = read_id),
    #   linewidth = 0.75
    # ) +
    scale_y_continuous(limits = ylims) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      title = paste("OD-time traces for Precision - ", assayname)
    )
  
  if (.ddf |> pull(read_path) |> unique() |> length() > 1) {
    p = p + facet_grid(row = vars(read_path))
  }

  
  pl = ggplotly(p, tooltip = "text") |>
    highlight(
      on        = "plotly_click",
      off       = "plotly_doubleclick",
      color     = "red",
      opacityDim = 0.1        # dims non-selected traces to 10% opacity
    ) |> 
    layout(
      annotations = list(
        list(
          text      = paste0("Hover over a trace to see run_name and recovery"),
          x         = 0,
          y         = 1.02,
          xref      = "paper",
          yref      = "paper",
          xanchor   = "left",
          showarrow = FALSE,
          font      = list(size = 11, color = "grey40")
        )
      )
    )
  saveWidget(pl |> partial_bundle(), paste0("./results/explorer_charts/OD_traces_precision_", sampless, "_", assayname, ".html"), selfcontained = TRUE)
  return()
}



# plot method comparison runs
assay_list = c("ALT","ALB","ALP","AST","TBIL","BUN","Ca","CREA","GLU","CHOL","TP","TRIG","AMY","CK","LAC","PHOS","Mg")

cross2(assay_df |> pull(sample) |> unique(), assay_list) |> 
  map(~{
    ## Analysis
    sampless = .x[[1]]
    assayname = .x[[2]]
    print(assayname)
    print(sampless)
    lambdas = assay_lambdas |> filter(assay == assayname) |> select(-assay) |> as.character() |> na.omit() |> as.character()

    .ddf = 
      df |> 
      filter(sample==sampless) |> 
      filter(assay == assayname) |> 
      filter(read_path %in% lambdas) |> 
      filter(read_id |> str_detect("^[12]r-")) |> 
      left_join(assay_df |> select(run_name, assay, vital, predicate), by=join_by("run_name", "assay")) |> 
      mutate(recovery = round(vital/predicate,2))

    ylims = .ddf |> 
      filter(str_detect(read_id, "^[12]r-")) |> 
      filter(is.finite(relOD)) |> 
      pull(relOD) |> 
      range(na.rm = TRUE) * c(0.95, 1.05)
    print(ylims)
    # create labels
    label_df = .ddf |> 
      group_by(run_name) |> 
      slice_max(tassay, n = 1) |> 
      ungroup()

    make_ggplotly_method_comparison(.ddf, label_df, ylims, assayname, sampless)
    return(NULL)
})
