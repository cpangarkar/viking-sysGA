files = list.files("./results/explorer_charts/", pattern = "\\.html$", full.names = FALSE)

# extract unique assay names
assays = files |>
  str_extract("(?<=OD_traces_method_comp_).+(?=\\.html)") |>
  na.omit() |>
  unique() |>
  sort()

rows = assays |>
  map_chr(~ paste0(
    '<tr>',
    '<td>', .x, '</td>',
    '<td><a href="OD_traces_method_comp_', .x, '.html">Method Comparison</a></td>',
    '<td>',
      '<a href="OD_traces_precision_Control L_', .x, '.html">Control L</a> &nbsp;',
      '<a href="OD_traces_precision_Control M_', .x, '.html">Control M</a> &nbsp;',
      '<a href="OD_traces_precision_Control H_', .x, '.html">Control H</a>',
    '</td>',
    '</tr>'
  )) |>
  paste(collapse = "\n")

html = paste0('
<html>
<head>
  <style>
    body { font-family: Arial, sans-serif; padding: 2em; }
    h2   { color: #333; }
    table { border-collapse: collapse; width: 70%; }
    th   { background-color: #f2f2f2; text-align: left; padding: 10px; border-bottom: 2px solid #ccc; }
    td   { padding: 8px 10px; border-bottom: 1px solid #eee; vertical-align: middle; }
    tr:hover { background-color: #f9f9f9; }
    a    { color: #2a7ae2; text-decoration: none; margin-right: 4px; }
    a:hover { text-decoration: underline; }
  </style>
</head>
<body>
  <h2>Viking SysGA Plots for Data Exploration</h2>
  <table>
    <tr>
      <th>Assay</th>
      <th>Method Comparison</th>
      <th>Precision</th>
    </tr>
    ', rows, '
  </table>
</body>
</html>')
writeLines(html, "./results/explorer_charts/index.html")
