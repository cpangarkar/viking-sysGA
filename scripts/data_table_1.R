
make_composite = function(df) {
  df = 
    df %>% 
    mutate(qc1 = ifelse(uF_QC_1 == "PASS", T, F)) %>% 
    mutate(qc2 = ifelse(uF_QC_2 == "PASS", T, F)) %>% 
    mutate(composite = ifelse(qc1 & is.na(qc2), vital_1,
                       ifelse(qc1 & qc2, vital_1, 
                       ifelse(!qc1 & is.na(qc2), NA, 
                       ifelse(qc1 & !qc2, vital_1,
                       ifelse(!qc1 & qc2, vital_2,
                       ifelse(is.na(qc1) & qc2, vital_2,
                       ifelse(is.na(qc1) & !qc2, NA, 
                       ifelse(!qc1 & !qc2, NA, NA)))))))) )

  
    return(df)
}
    
make_dtable = function(sample_name) {
  df <- 
    mcdf %>% 
    select(-OD) %>% 
    filter(sampleID == sample_name) %>% 
    mutate(vital = round(vital, 2)) %>% 
    mutate(perErr = round(100*(vital-predicate)/predicate,1)) %>% 
    mutate(absErr = vital-predicate) %>% 
    left_join(select(ainfo, meas, rounding, units, tae_per, tae_abs) , by=join_by(assay_name == meas)) %>% 
    mutate(within_TEa = ifelse(abs(absErr) <= tae_abs | abs(perErr) <= tae_per, T, F)) %>% 
    select(-perErr, -absErr, -tae_per, -tae_abs) %>% 
    pivot_wider(names_from = run_num, values_from = c(vital, within_TEa, uF_QC)) %>% 
    make_composite() %>% 
    mutate(within_TEa_1 = ifelse(is.na(within_TEa_1), T, within_TEa_1)) %>% 
    mutate(within_TEa_2 = ifelse(is.na(within_TEa_2), T, within_TEa_2)) %>% 
    mutate(predicate = round(predicate, rounding)) %>% 
    mutate(vital_1 = round(vital_1, rounding)) %>% 
    mutate(vital_2 = round(vital_2, rounding)) %>% 
    mutate(composite = round(composite, rounding))
  
  
    
    
  dtable = 
    df |>
    select(assay_name, predicate, vital_1, uF_QC_1, vital_2, uF_QC_2, composite, units, within_TEa_1, within_TEa_2, rounding) |>
    gt(
    )|>
    fmt_number(
      decimals= from_column(column = "rounding")
    )|>
    tab_style(
      style = list(
        cell_text(
          color = vpal["vgreen"],
          weight = "bold",
          size = px(22)
        ),
        cell_fill(
          color = "grey90"
        )
      ),
      locations = list(
        cells_title()
      )
    )|>
    tab_style(
      style = cell_text(size=px(14)),
      locations = cells_body(columns = everything())
    )|>
    cols_hide(columns = c(within_TEa_1, within_TEa_2, rounding)) |> 
    cols_align(
        align = "center"
    )|>
    data_color(
      columns = c("uF_QC_1", "uF_QC_2"), 
      target_columns = c("uF_QC_1", "uF_QC_2"),
      method = "factor",
      palette = c(vpal["vpoppy"],"gray60"),
      apply_to = "text"
    )|> 
    tab_header(
      title = paste0("Results for ", sample_name),
      subtitle = "                            "
    )|>
    sub_values(
      fn = function(x) is.na(x), 
      replacement=""
    )|>
    cols_label(
      assay_name = "",
      predicate = "",
      vital_1 = "Value",
      vital_2 = "Value",
      uF_QC_1 = "QC flag",
      uF_QC_2 = "QC flag",
      units = "",
      composite=""
    )|>
    tab_spanner(
      label = "Vital Run 1",
      columns = c("vital_1","uF_QC_1"),
      id = "num_spanner_1"
    )|>
    tab_spanner(
      label = "Vital Run 2",
      columns = c("vital_2","uF_QC_2"),
      id = "num_spanner_2"
    )|>
    tab_spanner(
      label = "Analyte",
      columns = c("assay_name"),
      id = "num_spanner_3"
    )|>
    tab_spanner(
      label = "Predicate",
      columns = c("predicate"),
      id = "num_spanner_4"
    )|>
    tab_spanner(
      label = "Unit",
      columns = c("units"),
      id = "num_spanner_5"
    )|>
    tab_spanner(
      label = "Composite",
      columns = c("composite"),
      id = "num_spanner_6"
    )|>
    tab_options(
      table.width = px(600),
      table.font.name = "Vital",
      container.height = px(600),
      container.padding.y = px(50),
      table.border.top.width = px(1),
      table.border.right.width = px(1),
      table.border.bottom.width = px(1),
      table.border.left.width = px(1)
    )
    
  
  dtable %>% gtsave(paste0("./results/table_", sample_name, ".pdf"))
 

}


# data_color(
#   columns = c("within_TEa_1", "within_TEa_2", "uF_QC_1", "uF_QC_2"), 
#   target_columns = c("vital_1", "vital_2", "uF_QC_1", "uF_QC_2"),
#   method = "factor",
#   palette = c(vpal["vpoppy"],"gray90"),
#   apply_to = "text"
# )|> 
