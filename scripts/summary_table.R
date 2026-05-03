



dtable = 
  sdf |>
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
  tab_style(
    style = cell_text(size=px(12), color=vpal["vgreen"]),
    locations = cells_column_labels(columns = everything())
  )|>
  tab_style(
    style = cell_text(size=px(10), color="grey45"),
    locations = cells_body(columns = "units")
  )|>
  cols_align(
    align = "center"
  )|>
  cols_hide(
    columns = "rounding"
  )|>
  tab_header(
    title = paste0("Summary Results")
  )|>
  sub_values(
    fn = function(x) is.na(x), 
    replacement=""
  )|>
  cols_label(
    assay_name = "",
    rounding = "",
    units = ""
  )|>
  cols_label(
    matches("composite") ~ "Vital",
    matches("predicate") ~ "Predicate"
  )|>
  tab_spanner(
    label = "Sample 1",
    columns = matches("L6_E1$"),
    id = "num_spanner_1"
  )|>
  tab_spanner(
    label = "Sample 2",
    columns = matches("L8_E3$"),
    id = "num_spanner_2"
  )|>
  tab_spanner(
    label = "Sample 3",
    columns = matches("L5_E2$"),
    id = "num_spanner_3"
  )|>
  tab_spanner(
    label = "Sample 4",
    columns = matches("L3_E1$"),
    id = "num_spanner_4"
  )|>
  tab_spanner(
    label = "Sample 5",
    columns = matches("L2_E2$"),
    id = "num_spanner_5"
  )|>
  tab_spanner(
    label = "Sample 6",
    columns = matches("L10_E7$"),
    id = "num_spanner_6"
  )|>
  tab_spanner(
    label = "Sample 7",
    columns = matches("Sample 7$"),
    id = "num_spanner_7"
  )|>
  tab_spanner(
    label = "Sample 8",
    columns = matches("Sample 8$"),
    id = "num_spanner_8"
  )|>
  tab_spanner(
    label = "Sample 9",
    columns = matches("Sample 9$"),
    id = "num_spanner_9"
  )|>
  tab_spanner(
    label = "Sample 10",
    columns = matches("Sample 10$"),
    id = "num_spanner_10"
  )|>
  tab_spanner(
    label = "Sample 11",
    columns = matches("Sample 11$"),
    id = "num_spanner_11"
  )|>
  data_color(
    columns = assay_name, 
    palette = vpal["vgreen"],
    apply_to = "text"
  )|>
  data_color(
    columns = assay_name, 
    palette = "#ccdcdb",
    apply_to = "fill"
  )|>
  tab_style(
    style = cell_borders(
      sides = c("left"),
      weight = px(1),
      color = "grey60"
      ),
    locations = cells_body(
      columns = matches("composite")
    )
  )|>
  tab_options(
    table.width = px(1100),
    table.font.name = "Vital",
    container.height = px(850),
    container.padding.y = px(00),
    table.border.top.width = px(1),
    table.border.right.width = px(1),
    table.border.bottom.width = px(1),
    table.border.left.width = px(1),
    page.orientation = "landscape"
  )

plot_table = T

if (plot_table) {dtable %>% gtsave(paste0("./results/table_summary.html"))}

