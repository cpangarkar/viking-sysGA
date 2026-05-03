library(xml2)
ares = c()

unzip("./config/Viking System GA Operator Log.xlsx", exdir = "./xlsx_extracted")

sheetnames = c("Erik", "Sigrid","Helga")

for (i in seq(length(sheetnames))) {
  # parse the relationships file
  xml = read_xml(paste0("./xlsx_extracted/xl/worksheets/_rels/sheet",i,".xml.rels"))

  # no namespace needed for rels files
  nodes = xml_find_all(xml, "//*[local-name()='Relationship']")
  links = tibble(
    r_id   = xml_attr(nodes, "Id"),
    target = xml_attr(nodes, "Target")
  ) |>
    filter(str_detect(target, "atlantis")) |>
    mutate(run_id = str_extract(target, "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"))

  sheet = read_xml(paste0("./xlsx_extracted/xl/worksheets/sheet",i,".xml"))

  # extract hyperlink nodes - cell ref to r:id mapping
  hl_nodes = xml_find_all(sheet, "//*[local-name()='hyperlink']")

  hyperlinks = tibble(
    cell = xml_attr(hl_nodes, "ref"),
    r_id = xml_attr(hl_nodes, "id")  # note: no namespace prefix needed here
  )

  cells = tidyxl::xlsx_cells("./config/Viking System GA Operator Log.xlsx", sheet=sheetnames[i])

  result = hyperlinks |>
    left_join(links, by = "r_id") |>
    left_join(
      cells |> select(address, character),
      by = c(cell = "address")
    )
  
  ares = rbind(ares, result)
}
  

ares = 
  ares |> 
  select(character, run_id) |> 
  rename(run_name = character) |> 
  filter(run_name |> str_detect("[EHS]-[0-9]{1,3}"))