make_dtable_single_sample = function(sample_name, rdf, plot_table = T) {
  df <- 
    rdf %>% 
    filter(sampleID == sample_name) %>% 
    filter(assay_name %in% c(assay_green, assay_yellow))
  
  if (df %>% select(run_num) %>% distinct() %>% pull() %>% length() == 1) {
    df = df %>% bind_rows(df[1,])
    df$run_num[nrow(df)] = 2
    df$vital[nrow(df)] = NA
    df$predicate[nrow(df)] = NA
    df$assay_name[nrow(df)] = "dummy"
    df$uF_QC[nrow(df)] = NA
    
  }
  df = df %>% 
    select(-run_id, -run_name, -module, -structure_label) %>% 
    mutate(vital = round(vital, 2)) %>% 
    mutate(perErr = round(100*(vital-predicate)/predicate,1)) %>% 
    mutate(absErr = vital-predicate) %>% 
    left_join(select(ainfo, meas, rounding, units, tae_per, tae_abs) , by=join_by(assay_name == meas)) %>% 
    mutate(within_TEa = ifelse(abs(absErr) <= tae_abs | abs(perErr) <= tae_per, T, F)) %>% 
    select(-perErr, -absErr, -tae_per, -tae_abs) %>% 
    distinct() %>% 
    pivot_wider(names_from = run_num, values_from = c(vital, within_TEa, uF_QC)) %>% 
    make_composite() %>% 
    mutate(within_TEa_1 = ifelse(is.na(within_TEa_1), T, within_TEa_1)) %>% 
    mutate(within_TEa_2 = ifelse(is.na(within_TEa_2), T, within_TEa_2)) %>% 
    mutate(predicate = round(predicate, rounding)) %>% 
    mutate(vital_1 = round(vital_1, rounding)) %>% 
    mutate(vital_2 = round(vital_2, rounding)) %>% 
    mutate(composite = round(composite, rounding)) %>% 
    filter(!is.na(predicate))
  
  
  dtable = 
    df |>
    select(assay_name, predicate, vital_1, uF_QC_1, vital_2, uF_QC_2, composite, units, within_TEa_1, within_TEa_2, rounding) |>
    gt(
    )|>
    cols_move(
      predicate, composite
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
        cells_title(group="title")
      )
    )|>
    tab_style(
      style = list(
        cell_text(
          color = "grey70",
          size = px(16)
        ),
        cell_fill(
          color = "grey90"
        )
      ),
      locations = list(
        cells_title(group="subtitle")
      )
    )|>
    tab_style(
      style = cell_text(size=px(12)),
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
      palette = c("gray70",vpal["vpoppy"]),
      ordered = T,
      apply_to = "text"
    )|> 
    tab_header(
      title = paste0("Results for ", sample_name),
      subtitle = "Archived serum sample ID - xx-xxxx"
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
      container.height = px(800),
      container.padding.y = px(0),
      table.border.top.width = px(1),
      table.border.right.width = px(1),
      table.border.bottom.width = px(1),
      table.border.left.width = px(1)
    )
  
  
  if (plot_table) {dtable %>% gtsave(paste0("./results/table_", sample_name, ".pdf"))}
  
   return(df |> select(assay_name, predicate, composite, units, rounding))
   
}