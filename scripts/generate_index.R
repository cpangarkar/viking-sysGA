files = list.files("./results/explorer_charts/", pattern = "\\.html$", full.names = FALSE)

# extract unique assay names
assays = files |>
  str_extract("(?<=OD_traces_method_comp_|OD_traces_precision_).+(?=\\.html)") |>
  na.omit() |>
  unique() |>
  sort()

rows = assays |>
  map_chr(~ paste0(
    '<tr>',
    '<td>', .x, '</td>',
    '<td><a href="OD_traces_method_comp_', .x, '.html">Method Comparison</a></td>',
    '<td><a href="OD_traces_precision_', .x, '.html">Precision</a></td>',
    '</tr>'
  )) |>
  paste(collapse = "\n")

html = paste0('
<html>
<head>
  <style>
    body { font-family: Arial, sans-serif; padding: 2em; }
    h2   { color: #333; }
    table { border-collapse: collapse; width: 60%; }
    th   { background-color: #f2f2f2; text-align: left; padding: 10px; border-bottom: 2px solid #ccc; }
    td   { padding: 8px 10px; border-bottom: 1px solid #eee; }
    tr:hover { background-color: #f9f9f9; }
    a    { color: #2a7ae2; text-decoration: none; }
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
