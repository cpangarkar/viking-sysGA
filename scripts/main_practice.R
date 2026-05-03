rm(list=ls())
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


path = "C:/Users/cpang/Documents/demo_KH//"
setwd(path)
source("./scripts/config.R")
source("./scripts/library_functions.R")

##Inputs and Config##

mcdf = read_csv("./data/practice_data.csv")
mcdf = mcdf %>% filter(grepl("demo",tolower(sampleID)))
caldf = read_csv("./data/cal_data.csv")

ainfo = read_xlsx("./config/TAE_summary_Mystique.xlsx", sheet = "Sheet1", range = "A1:X46")

corr = read_csv("./data/corr_factors.csv")

assay_list = c( "GLU","Ca","TRIG","ALB","CHOL","TP","CREA","ALP","GGT","TBIL","AST")

cdf1 = read_xlsx("./data/Alpha plus data compilation.xlsx", sheet="bias correction", range = "BU169:CV198", col_names = F, col_types = "numeric")
ddf = read_xlsx("./data/Alpha plus data compilation.xlsx", sheet="bias correction", range = "A169:A198", col_names = F, col_types = "date")

cdf1 = bind_cols(ddf, cdf1)
cdf1 = cdf1 %>% 
      select(-c(...13, ...17, ...23,...28,...29)) %>% 
      magrittr::set_names(c("date","GLU","Ca","TRIG","ALB","CHOL","TP","CREAT","ALP","GGT","AST","TBIL","DDIM","hsCRP","TSH","HGB","WBC","LYMPH%","NEUT%","MONO%","EOS%","RBC","HCT","MCV")) %>% 
      mutate(srow_num = seq(169:198)) %>% 
      mutate(method = "vital") %>% 
      pivot_longer(c(-date, -srow_num, -method))

cdf2 = read_xlsx("./data/Alpha plus data compilation.xlsx", sheet="bias correction", range = "CW169:DX198", col_names = F, col_types = "numeric")

cdf2 = bind_cols(ddf, cdf2)
cdf2 = cdf2 %>% 
      select(-c(...3, ...16, ...17,...23,...28,...29)) %>% 
      magrittr::set_names(c("date","GLU","Ca","TRIG","ALB","CHOL","TP","CREAT","ALP","GGT","AST","TBIL","hsCRP","TSH","HGB","WBC","LYMPH%","NEUT%","MONO%","EOS%","RBC","HCT","MCV")) %>% 
      mutate(hsCRP = hsCRP*10) %>% 
      mutate(srow_num = seq(169:198)) %>% 
      mutate(method = "predicate") %>% 
      pivot_longer(c(-date, -srow_num, -method))


  
cdf = bind_rows(cdf1, cdf2) %>% 
      pivot_wider(names_from = method, values_from = value) %>% 
      filter(!is.na(predicate)) %>% 
      filter(!is.na(vital)) %>% 
      group_by(name) %>% 
      group_modify(~{detect_outliers(.x)}) %>% 
      filter(vital > 0)

#windows(18,18)
p = 
  cdf %>% 
  filter(date > as.Date("2024-05-06")) %>% 
  filter(srow_num <27) %>%  
  filter(!outlier) %>% 
  ggplot(aes(x=predicate, y=vital)) +  
  geom_abline(slope=1, intercept=0, col=vpal["vpoppy"]) + 
  stat_poly_line(formula = y ~ x, method="lm") + 
  stat_poly_eq(formula = y ~ x, use_label(c("eq","R2")), coef.digits=4) +
  geom_point(aes(col=factor(date), shape=outlier)) +
  facet_wrap(.~name, scales="free") +
  theme_bw() + 
  labs(
    title = "Libra May 9,10,11 runs concordance",
    subtitle = "poppy is y = x"
  )

#print(p)

fits = 
  cdf %>% 
  filter(date > as.Date("2024-05-06")) %>% 
  #filter(srow_num <27) %>% 
  filter(!outlier) %>% 
  group_by(name) %>% 
  group_modify(~{
    mod = with(.x, lm(vital ~ predicate))
    dfrow = data.frame(t(c(coef(mod), summary(mod)$r.squared)))
    colnames(dfrow) = c("intercept","slope","R2")
    return(dfrow)
  }) %>% 
  ungroup()

fits = read_xlsx("./bias_corr_analysis.xlsx", sheet="bias_corr", range="A2:D24")  

cdf = 
  cdf %>% 
  rename(assay_name = name) %>% 
  left_join(select(ainfo, meas, rounding, units, tae_per, tae_abs) , by=join_by(assay_name == meas)) 

cdf = 
  cdf %>% 
  group_by(assay_name) %>% 
  group_modify(~{
    coefs = fits %>% filter(name == .x$assay_name[1]) %>% select(intercept, slope)
    .x = .x %>% 
            mutate(vital_orig = vital) %>% 
            mutate(vital = (vital * coefs$slope[1]) + coefs$intercept[1])
    return(.x %>% select(-assay_name))
  }, .keep=T) 
  


  
# dum = 
#   cdf %>% 
#   group_by(assay_name) %>% 
#   group_map(~{
#     
#     print(.x)
#     assay = .x$assay_name[1]
#     pd    = compile_tae_limits(ainfo %>% filter(meas == assay))
#     lims  = pd$lims
#     tae   = pd$tae
#     mdl   = refr = pd$refr
#     
#     scdf = .x %>% 
#            filter(srow_num > 26) %>% 
#            mutate(uF_QC = "PASS")
#     
#     if (nrow(scdf)<2) {return()}
#     
#     mod  = with(filter(scdf), lm(vital ~ predicate))
#     
#     for (modifier in c(3)) {
#         print(modifier)
#         make_vitalPlot_Hdemo(scdf, assay, tae, lims, mdl, ainfo %>% filter(meas == assay), mod, modifier)
#     }
#     return()
#   }, .keep=T)

#pdf_combine(dir("./results/", pattern = "MethodComparison_Plots_*.*type-1",full.names = T), "output_plots_type-1.pdf")
#pdf_combine(dir("./results/", pattern = "MethodComparison_Plots_*.*type-2",full.names = T), "output_plots_type-2.pdf")
#pdf_combine(dir("./results/", pattern = "MethodComparison_Plots_*.*type-3",full.names = T), "output_plots_type-3.pdf")
#pdf_combine(dir("./results/", pattern = "MethodComparison_Plots_*.*type-4",full.names = T), "output_plots_type-4.pdf")

# pdf_file_list = dir("./results/", pattern = "MethodComparison_Plots_*.*type-3",full.names = T)
# pdf_combine(pdf_file_list, "output_plots_type-3.pdf")
# file.remove(pdf_file_list)



# plist =
#   assay_list %>%
#   map(~{
#     print(.x)
#     recalibrate(.x, 1)
#     recalibrate(.x, filter(corr, assay_name==.x)$corr)
#     return()
#   }, .keep=T)
# 
# pdf_combine(dir("./results/", pattern = "MethodComparison_Plots",full.names = T), "output_plots.pdf")

# sample_list = 
#   mcdf %>% 
#   select(sampleID) %>% 
#   distinct() %>% 
#   pull() 
# 
# plist = 
#   sample_list %>% 
#   map(~{
#     print(.x)
#     make_dtable(.x)
#     return()
#   }, .keep=T)
# 
# flist = dir("./results/", pattern = "table", full.names = T)
# flist = flist[sort(as.numeric(str_extract(flist,"[0-9]{1,2}")), index.return=T)$ix]
# pdf_combine(flist, "output_tables.pdf")

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
# plist =
#   assay_list %>%
#   map(~{
#     print(.x)
#     recalibrate(.x, 1)
#     print(1)
#     recalibrate(.x, filter(corr, assay_name==.x)$corr)
#     print(2)
#     recalibrate(.x, with(filter(corr, assay_name==.x), c(corr_int, corr_slope)))
#     return()
#   }, .keep=T)
# 
# pdf_combine(dir("./results/", pattern = "MethodComparison_Plots",full.names = T), "output_plots.pdf")

# sample_list = 
#   mcdf %>% 
#   select(sampleID) %>% 
#   distinct() %>% 
#   pull() 
# 
# df = 
#   sample_list %>% 
#   map_dfr(~{
#     print(.x)
#     ddf = make_dtable(.x, plot_table = F)
#     ddf = ddf %>% mutate(sampleID = .x)
#     return(ddf)
#   }, .keep=T)
# 
# flist = dir("./results/", pattern = "table", full.names = T)
# flist = flist[sort(as.numeric(str_extract(flist,"[0-9]{1,2}")), index.return=T)$ix]
# pdf_combine(flist, "output_tables.pdf")
