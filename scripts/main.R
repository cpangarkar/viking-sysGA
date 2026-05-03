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


path = "C:/Users/cpang/analysis/viking-sysGA/"
setwd(path)
source("./scripts/config.R")
source("./scripts/library_functions.R")

##Inputs and Config##

assay_list = c("ALT","ALB","ALP","AST","TBIL","BUN","Ca","CREA","GLU","CHOL","TP","TRIG","AMY","CK","LAC","PHOS","Mg")
aldf = as_tibble(data.frame(assay_name = assay_list))
ainfo = read_xlsx("./config/TAE_summary_Mystique.xlsx", sheet = "Sheet1", range = "A1:W24")
vis_outliers = 
  read_csv("./data/visual-outliers.csv") |> 
  pivot_longer(cols = everything(), names_to = "assay", values_to = "disc") |>
  filter(!is.na(disc)) |> 
  mutate(comstring = paste0(disc, assay, sep="-"))

samplelist = 
  read_csv("./config/sample-list.csv") |> 
  mutate(sample = sample |> str_replace("_","-")) |> 
  mutate(species = ifelse(grepl("^fe", tolower(species)), "feline","canine"))

samplehil = 
  read_csv("./config/sample-hil.csv") |> 
  magrittr::set_colnames(c("sample","hil","species")) |> 
  select(-species) |> 
  separate(hil, into=c("hem","ict","lip"), sep="\\s+") |> 
  mutate(hem = as.numeric(hem)) |> 
  mutate(sample = sample |> str_replace("_","-")) 

prec_levels =
  read_csv("./config/precision_levels.csv")

## Read main method comparison data

fL = "./data/OD_corrected 30apr am.csv"
#fL = "./data/data_for_luc 29apr am.csv"

df =
  read_csv(fL) |> 
  select(-contains("well")) |> 
  pivot_longer(cols=-c(disc, sample, timestamp, contains("well")),names_to = "assay", values_to = "value") |> 
  mutate(name=assay) |> 
  mutate(name = assay |> str_split_i("_",2)) |> 
  mutate(name = ifelse(is.na(name), "vital","predicate")) |> 
  rowwise() |> 
  mutate(assay = assay |> str_split_i("_",1) ) |> 
  pivot_wider(id_cols= c(disc, sample, timestamp, assay), names_from=name, values_from=value, values_fn = first) |> 
  mutate(uF_QC = "PASS") |> 
  mutate(outlier = F) |> 
  mutate(instr = disc |> str_split_i("-",1)) |> 
  rowwise() |> 
  mutate(outlier = ifelse(paste0(disc, assay, sep="-") %in% vis_outliers$comstring, T, F)) |> 
  left_join(samplelist, by="sample") |> 
  left_join(samplehil, by="sample")


df = 
    df |> 
    filter(sample |> str_detect("^\\d{7}-\\d{2}$"))  |> 
    group_by(assay)

## Plot method comp
cross2(c("canine","feline"), assay_list) %>%
  map(~{
    print(c(.x[[2]], .x[[1]]))
    make_plots_final(df %>% filter(assay==.x[[2]] & species == .x[[1]]))
    })

pdf_file_list = dir("./results/", pattern = "MethodComparison_Plots_*.*type-3",full.names = T)

analyte = pdf_file_list %>% str_split_i("/",3) %>% str_split_i("_",3) 
dddf =  as_tibble(data.frame(fL = pdf_file_list, assay = analyte))
        

pdf_combine(dddf$fL, "MC_output_plots_type - With Spectral Correction 30Apr AM.pdf")
file.remove(pdf_file_list)

## do initial method comp fits and set up bias correction
bias_table = 
  cross2(c("canine","feline"), assay_list) %>%
    map(~{
      print(c(.x[[2]], .x[[1]]))
      ainfo = ainfo |> filter(meas == .x[[2]])  
      pd = compile_tae_limits(ainfo)

      rdf = 
        df |> 
        filter(assay==.x[[2]] & species == .x[[1]]) |> 
        filter(!is.na(vital)) |> 
        filter(vital > 0) |> 
        filter(predicate > 0) |>
        filter(predicate >= ainfo$amr_low & predicate <= ainfo$amr_high)
      
      ## get non-regression stats


      rdf = 
        rdf |> 
        rowwise() |> 
        mutate(tea_max = predicate*(1 + 0.01*pd$tae$errlims[which.min(abs(pd$tae$xr - predicate))])) |> 
        mutate(tea_min = predicate*(1 - 0.01*pd$tae$errlims[which.min(abs(pd$tae$xr - predicate))])) |> 
        mutate(within_tea = ifelse(is.na(tea_max), T, ifelse(vital |> between(tea_min, tea_max), T, F))) |> 
        ungroup()


      withins =
        rdf |> 
        summarise(
          num_within = length(which(within_tea)),
          per_within = round(100*length(which(within_tea))/n(),1),
          tot_points = n()
        )


      if (nrow(rdf)>2) {
        mod = with(rdf |> filter(!outlier), mcreg(predicate, vital, error.ratio = 1, method.reg = "WDeming"))
        mdls = ainfo |> filter(meas==.x[[2]]) |> select(contains("ref"), contains("mdl")) |> as.numeric() |> na.omit() |> (\(x) x[x != 0])()
        mdls = mdls[(mdls >= ainfo$amr_low & mdls <= ainfo$amr_high)]

        bias = mod |> calcBias(x.levels = mdls, type="proportional") |> as_tibble() |> magrittr::set_colnames(c("level","bias_per","se","lci","uci"))
        rdf = 
          rdf |> 
          mutate(vitalc = (vital - mod@glob.coef[1])/mod@glob.coef[2]) |> 
          mutate(biascorr_slope = mod@glob.coef[2]) |> 
          mutate(biascorr_int = mod@glob.coef[1])

      } else {
        bias = bind_cols(level = c(NA), bias_per = c(NA), se = c(NA), lci = c(NA), uci = c(NA))
        rdf = 
          rdf |> 
          mutate(vitalc = (vital - 0)/1)
      }
      bias = 
        bias |> 
        mutate(assay = .x[[2]]) |> 
        mutate(species = .x[[1]]) |> 
        mutate(level = round(level, ainfo$rounding)) |> 
        mutate(bias_per = round(bias_per, 1)) |> 
        mutate(lci = round(lci, 1)) |> 
        mutate(uci = round(uci, 1)) |> 
        mutate(num_within = withins$num_within) |> 
        mutate(per_within = withins$per_within) |> 
        mutate(tot_points = withins$tot_points)

      return(list(rdf = rdf, bias = bias))

    }) -> ans


## Do bias correction and re-plot 
corr_df = map(ans, 1) |> list_rbind()
corr_df = 
  corr_df |> 
  select(-vital) |> 
  rename(vital = vitalc) |> 
  relocate(vital, .before=predicate)

df = 
  df |> 
  left_join(corr_df |> select(disc, assay, vital, biascorr_slope, biascorr_int), by=join_by("disc","assay")) |> 
  rename(vital = vital.x) |> 
  rename(vitalc = vital.y) |> 
  relocate(vitalc, .before=predicate)

  cross2(c("canine","feline"), assay_list) %>%
    map(~{
      print(c(.x[[2]], .x[[1]]))
      make_plots_final(corr_df %>% filter(assay==.x[[2]] & species == .x[[1]]))
      })
  
  pdf_file_list = dir("./results/", pattern = "MethodComparison_Plots_*.*type-3",full.names = T)
  
  analyte = pdf_file_list %>% str_split_i("/",3) %>% str_split_i("_",3) 
  dddf =  as_tibble(data.frame(fL = pdf_file_list, assay = analyte))
          
  
  pdf_combine(dddf$fL, "MC_output_plots_type - With Spectral Correction 30Apr AM - Bias Correction.pdf")
  file.remove(pdf_file_list)

## Now re-caculate all bias metrics by fitting meth comp on corrected data



bias_table = 
  cross2(c("canine","feline"), assay_list) %>%
    map(~{
      print(c(.x[[2]], .x[[1]]))
      ainfo = ainfo |> filter(meas == .x[[2]])  
      pd = compile_tae_limits(ainfo)

      rdf = 
        df |> 
        filter(assay==.x[[2]] & species == .x[[1]]) |> 
        filter(!is.na(vitalc)) |> 
        filter(vitalc > 0) |> 
        filter(predicate > 0) |>
        filter(predicate >= ainfo$amr_low & predicate <= ainfo$amr_high)

      ## get non-regression stats
      
      rdf = 
        rdf |> 
        rowwise() |> 
        mutate(tea_max = predicate*(1 + 0.01*pd$tae$errlims[which.min(abs(pd$tae$xr - predicate))])) |> 
        mutate(tea_min = predicate*(1 - 0.01*pd$tae$errlims[which.min(abs(pd$tae$xr - predicate))])) |> 
        mutate(within_tea = ifelse(is.na(tea_max), T, ifelse(vitalc |> between(tea_min, tea_max), T, F))) |> 
        ungroup()


      withins =
        rdf |> 
        summarise(
          num_within = length(which(within_tea)),
          per_within = round(100*length(which(within_tea))/n(),1),
          tot_points = n()
        )


      if (nrow(rdf)>2) {
        mod = with(rdf |> filter(!outlier), mcreg(predicate, vitalc, error.ratio = 1, method.reg = "WDeming"))
        mdls = ainfo |> filter(meas==.x[[2]]) |> select(contains("ref"), contains("mdl")) |> as.numeric() |> na.omit() |> (\(x) x[x != 0])()
        mdls = c(mdls, prec_levels |> filter(meas == .x[[2]]) |> pull(level) |> as.numeric())
        mdls = mdls[(mdls >= ainfo$amr_low & mdls <= ainfo$amr_high)]

        bias = mod |> calcBias(x.levels = mdls, type="proportional") |> as_tibble() |> magrittr::set_colnames(c("level","bias_per","se","lci","uci"))
        rdf = 
          rdf |> 
          mutate(vitalcc = (vitalc - mod@glob.coef[1])/mod@glob.coef[2]) |> 
          mutate(biascorr_slope = mod@glob.coef[2]) |> 
          mutate(biascorr_int = mod@glob.coef[1])

      } else {
        bias = bind_cols(level = c(NA), bias_per = c(NA), se = c(NA), lci = c(NA), uci = c(NA))
        rdf = 
          rdf |> 
          mutate(vitalcc = (vitalc - 0)/1)
      }
      bias = 
        bias |> 
        mutate(assay = .x[[2]]) |> 
        mutate(species = .x[[1]]) |> 
        mutate(level = round(level, ainfo$rounding)) |> 
        mutate(unit = ainfo$units) |> 
        mutate(bias_per = round(bias_per, 1)) |> 
        mutate(lci = round(lci, 1)) |> 
        mutate(uci = round(uci, 1)) |> 
        mutate(num_within = withins$num_within) |> 
        mutate(per_within = withins$per_within) |> 
        mutate(tot_points = withins$tot_points)

      return(list(rdf = rdf, bias = bias))

    }) -> ans

    corr_bias_df = 
      map(ans, 2) |> 
      list_rbind() |> 
      select(-se) |> 
      relocate(assay, .before=level) |> 
      write_csv("./results/bias_table.csv")
