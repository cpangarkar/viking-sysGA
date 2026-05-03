# 
# 
# mcdf = read_csv("./data/practice_data.csv")
# mcdf = mcdf %>% filter(grepl("demo",tolower(sampleID)))
# caldf = read_csv("./data/cal_data.csv")
# 
# ainfo = read_xlsx("./config/TAE_summary_Mystique.xlsx", sheet = "Sheet1", range = "A1:X25")
# 
# corr = read_csv("./data/corr_factors.csv")
# 
# assay_list = c( "GLU","Ca","TRIG","ALB","CHOL","TP","CREA","ALP","GGT","TBIL","AST")
# 
# # plist =
# #   assay_list %>%
# #   map(~{
# #     print(.x)
# #     recalibrate(.x, 1)
# #     print(1)
# #     recalibrate(.x, filter(corr, assay_name==.x)$corr)
# #     print(2)
# #     recalibrate(.x, with(filter(corr, assay_name==.x), c(corr_int, corr_slope)))
# #     return()
# #   }, .keep=T)
# # 
# # pdf_combine(dir("./results/", pattern = "MethodComparison_Plots",full.names = T), "output_plots.pdf")
# 
# sample_list =
#   mcdf %>%
#   select(sampleID) %>%
#   distinct() %>%
#   pull()
# 
# tdf =
#   sample_list %>%
#   map_dfr(~{
#     print(.x)
#     ddf = make_dtable_single_sample(.x, mcdf, plot_table = F)
#     ddf = ddf %>% mutate(sampleID = .x)
#     return(ddf)
#   }, .keep=T)
# #
# # flist = dir("./results/", pattern = "table", full.names = T)
# # flist = flist[sort(as.numeric(str_extract(flist,"[0-9]{1,2}")), index.return=T)$ix]
# # pdf_combine(flist, "output_tables.pdf")
# 

sdf =
  sample_list %>%
  map_dfr(~{
    print(.x)
    .y = make_dtable_single_sample(.x, rdf, F)
    .y = mutate(.y, sampleID = .x)
  })

sdf = 
  sdf %>% 
  pivot_wider(names_from = sampleID, values_from = c("predicate","composite"), names_vary="slowest")

aa = colnames(sdf)
aa = aa[-c(1:3)]
aa = aa[sort(as.numeric(str_extract(aa,"[0-9]{1,2}")), index.return=T)$ix]
aa = as_tibble(data.frame(aa))
aa = aa %>% 
      mutate(sid = as.numeric(str_extract(aa,"[0-9]{1,2}"))) %>% 
      group_by(sid) %>% 
      group_modify(~{arrange(.x, aa)}) %>% 
      pull(aa)
ii = as.numeric(sapply(aa, function(x) which(colnames(sdf)==x)))

sdf = sdf[,c(1,2,3,ii)]
  
#windows()
sdf %>% gt()
