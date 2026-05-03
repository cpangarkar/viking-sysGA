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
library(broom)
library(ggrepel)
library(ggh4x)


path = "C:/Users/cpang/analysis/viking-sysGA/"
setwd(path)
source("./scripts/config.R")
source("./scripts/library_functions.R")

##Inputs and Config##

assay_list = c("ALT","ALB","ALP","AST","TBIL","BUN","Ca","CREA","GLU","CHOL","TP","TRIG","AMY","CK","LAC","PHOS","Mg")
aldf = as_tibble(data.frame(assay_name = assay_list))
ainfo = read_xlsx("./config/TAE_summary_Mystique.xlsx", sheet = "Sheet1", range = "A1:W24")
prec_levels =
  read_csv("./config/precision_levels.csv")



## Read main method comparison data


fL = "./data/precision-bias-corrected-partial.csv"

df =
  read_csv(fL) |>  
  pivot_longer(cols=-c(disc, sample, timestamp, contains("well")),names_to = "assay", values_to = "value") |> 
  mutate(name=assay) |> 
  mutate(name = assay |> str_split_i("_",2)) |> 
  mutate(name = ifelse(is.na(name), "vital","predicate")) |> 
  rowwise() |> 
  mutate(assay = assay |> str_split_i("_",1) ) |> 
  pivot_wider(id_cols= c(disc, sample, timestamp, assay), names_from=name, values_from=value, values_fn = first) |> 
  filter(assay %in% assay_list) |> 
  filter(sample |> tolower()|> str_detect("control")) |> 
  mutate(date = as.Date(timestamp)) |> 
  mutate(day = timestamp |> day()) |> 
  mutate(instr = disc |> str_split_i("-",1))

df = 
  df |> 
  filter(disc != "E-169") |> 
  filter(disc != "E-170")
# ouliers


.x = filter(df, assay=="CREA" & sample == "Control L")
prec_df = 
  df |> 
    group_by(assay, sample) |> 
    group_modify(~{
      if (.x$assay[1] =="CK") {
        .x = .x |> filter(disc != "E-180")
      }
      mod = anovaVCA(vital ~ instr + day, as.data.frame(.x))
      vc <- mod$aov.tab |> 
            as.data.frame() |> 
            rownames_to_column("component") |> 
            as_tibble() |> 
            select(component, SD, `CV[%]`) |> 
            magrittr::set_colnames(c("component","sd","cv"))

      vc = 
        vc |>
        filter(component != "error") |>       # drop error or keep, your choice
        pivot_wider(
          names_from = component,
          values_from = c(sd, cv),
          names_glue = "{component}_{.value}"
        ) |>
        select(
          total_sd, total_cv,
          instr_sd, instr_cv,
          day_sd, day_cv 
        ) |> 
        mutate(level = mean(.x$vital, na.rm=T)) |> 
        mutate(N = nrow(.x |> filter(!is.na(vital)))) |> 
        relocate(level, .before=total_sd) |> 
        relocate(N, .after=level)
      return(vc)
    },.keep=T)

prec_df |>
  group_by(assay) |> 
  group_modify(~{
    .x = 
      .x |> 
      mutate(sample = factor(sample, levels = c("Control L", "Control M", "Control H"))) |>
      arrange(sample)
  })

prec_df |> write_csv("./results/precision_table.csv")

pdf("./results/precision_charts.pdf", width=8, height=5)

  df |> 
    group_by(assay, sample) |> 
    group_modify(~{
      if (.x$assay[1] =="CK") {
        .x = .x |> filter(disc != "E-180")
      }
      try(varPlot(vital ~ instr + day, as.data.frame(.x), Title = list(paste0(.x$assay[1],"-",.x$sample[1]))))
      return(data.frame())
      
    },.keep=T)

dev.off()



pdf("./results/precision_charts_ggplot.pdf", width=11, height=5)

  plist = 
    df |> 
    group_by(assay, sample) |> 
    group_map(~{
      ylims = c(0.9, 1.1) * range(.x$vital, na.rm=T)
      .m   <- mean(.x$vital, na.rm = TRUE)
      .s   <- sd(.x$vital, na.rm = TRUE)
      .out <- .x |> filter(abs(vital - .m) > 2.5 * .s)
      p = 
        .x |> 
        ggplot(aes(x = 1, y=vital)) + 
        geom_hline(yintercept = .m, col=vpal["vslate_dark"], linewidth=0.5) + 
        geom_hline(
          yintercept = c(
            .m + 1*.s, .m - 1*.s,
            .m + 2*.s, .m - 2*.s
            ), 
          col=vpal["vslate"], linewidth=0.5, linetype="dashed") + 
        geom_text(
            data  = .out,
            aes(label = disc),
            hjust = -0.2,
            size  = 3,
            col   = vpal["vpoppy"]
          ) + 
        geom_jitter(col = vpal["vgreen"], size=1.4, shape = 21, width = 0.02, height = 0, fill = vpal["vleaf"]) + 
        stat_summary(fun = mean, geom = "point", shape=24, size = 1.5, col = vpal["vpoppy"], fill = vpal["vpoppy"]) + 
        scale_y_continuous(limits = ylims) + 
        scale_x_continuous(limits = c(0.5, 1.5)) +
        facet_nested(
          . ~ instr + day,  
          labeller = labeller(
            instr = as_labeller(\(x) paste("Instrument", x)),
            day   = as_labeller(\(x) paste("Day", x))
          )
        ) +
        theme_bw() + 
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) + 
        labs(
        title = paste0(.x$assay[1],"-",.x$sample[1]),
        subtitle = "Individual measurements shown in green, mean of day in poppy triangles. Grey line shows global mean,\n with dashed lines as +/- 2SD. Points outside of +/- 2.5SD are labeled for intvestigation. Day refers to the date in Apr."
        ) + 
        theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.subtitle = element_text(size = 8),
        strip.background = element_rect(fill = "grey80", color = "grey70"),
        panel.background = element_rect(fill = "white"),
        panel.border     = element_rect(color = "grey80")
        )
      return(p)
    },.keep=T)

  plist |> walk(ggplot2:::print.ggplot)

dev.off()
