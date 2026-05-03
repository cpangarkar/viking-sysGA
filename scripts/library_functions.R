get_assay_info <- function(measurand, ainfopath) {
  ainfo = read_xlsx(ainfopath, 1)
  ainfo = filter(ainfo, meas==measurand)
  return(ainfo)
}

compile_tae_limits <- function(ainfo) {
  lims = c(ainfo$plot_lim_min, ainfo$plot_lim_max)
  refr = c(ainfo$ref_low, ainfo$ref_high)
  xr   = seq(0.5*lims[1], 2*lims[2], length.out=300)
  
  ixr = xr |> cut(breaks = c(0.5*lims[1], ainfo$ref_low, ainfo$ref_high, 2*lims[2]), labels=seq(3), include.lowest = T) |> as.numeric()
  errlims = with(ainfo, c(tea_low, tea_ri, tea_high)) |> as.numeric()

  ytae_min = xr*(1 - 0.01*errlims[ixr])
  ytae_max = xr*(1 + 0.01*errlims[ixr])

  atae_min = -xr*0.01*errlims[ixr]
  atae_max =  xr*0.01*errlims[ixr]
  
  tae = data.frame(xr, ytae_min, ytae_max, 
                   atae_min, atae_max, errlims = errlims[ixr])
  return(list(tae=tae, lims=lims, refr=refr))
}

parse_au_predicate_data = function(audf_path) {
  
  au_ts = audf_path %>% str_split_i("/",4) %>% str_split_i("\\.",1) %>% str_split("_")
  au_ts = parse_date_time(paste0(au_ts[[1]][2]," ",au_ts[[1]][3]), "ymd HM")
  au_ts = force_tz(au_ts, "America/New_York")
  
  audf = 
    read_csv(audf_path) %>% 
    select(S.ID, matches("^[0-9]{1,2}\\.[[:upper:]]"))
  
  raudf = 
    audf %>% 
    magrittr::set_names(colnames(audf) %>% str_split_i("\\.",2)) %>% 
    mutate(ID = as.character(ID)) %>% 
    separate(LIH, into=c("Lip","Ict","Hem"), sep = " ") %>% 
    mutate_at(vars(!starts_with("ID")), as.numeric) %>% 
    rename(sampleID = ID) %>% 
    pivot_longer(cols = where(is.numeric)) %>% 
    rename(au_name = name)
  
  au_name_key = read_csv("./config/au_name_key.csv")
  
  raudf = 
    raudf %>% 
    left_join(au_name_key, by="au_name") %>% 
    select(-au_name) %>% 
    filter(assay_name != "exclude") %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    mutate(tstamp = au_ts) %>% 
    mutate(instr = "BC AU480")
  
  
  return(raudf)
}

parse_dxh_predicate_data = function(dxdf_path) {
  
  dx_ts = dxdf_path %>% str_split_i("/",4) %>% str_split_i("\\.",1) %>% str_extract(" 500_2024-*.*m$") %>% str_replace("500_","")
  dx_ts = dx_ts %>% parse_date_time("y-m-d_HM")
  dx_ts = force_tz(dx_ts, "America/New_York")
  
  dxdf = 
    read_csv(dxdf_path) %>% 
    select(Specimen_ID,WBC,RBC,HGB,HCT,MCV,MCH,MCHC,RDW,RDWSD,PLT,LY,MO,NE,EO,BA) %>% 
    rename(sampleID = Specimen_ID) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    pivot_longer(cols = where(is.numeric)) %>% 
    rename(dxh_name = name)
  
  dxh_name_key = read_csv("./config/dxh_name_key.csv")
  
  rdxdf = 
    dxdf %>% 
    left_join(dxh_name_key, by="dxh_name") %>% 
    select(-dxh_name) %>% 
    filter(assay_name != "exclude")%>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    mutate(tstamp = dx_ts) %>% 
    mutate(instr = "BC DxH500")
  
  return(rdxdf)
}

parse_ac_predicate_data = function(acdf_path) {
  
  ac_ts = acdf_path %>% str_split_i("/",4) %>% str_split_i("\\.",1)
  ac_ts = ac_ts %>% parse_date_time("m-d-y-H-M-S")
  ac_ts = force_tz(ac_ts, "America/New_York")
  
  acdf = 
    read_csv(acdf_path) %>% 
    select(`Sample ID`, `Test Name`, Result) %>% 
    magrittr::set_names(c("sampleID","ac_name","value"))
  
  ac_name_key = read_csv("./config/ac_name_key.csv")
  
  racdf = 
    acdf %>% 
    left_join(ac_name_key, by="ac_name") %>% 
    select(-ac_name) %>% 
    filter(assay_name != "exclude")%>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    mutate(tstamp = ac_ts) %>% 
    mutate(instr = "BC Access")
  
  return(racdf)
}

compile_predicate_results = function() {
  
  audf = "./data/predicate/" %>% 
          dir(pattern="^AU)*.*csv$", full.names = T) %>% 
          map_dfr(~{parse_au_predicate_data(.x)})
        
  acdf = "./data/predicate/" %>% 
          dir(pattern="^05-*.*csv$", full.names = T) %>% 
          map_dfr(~{parse_ac_predicate_data(.x)})
        
  dxdf = "./data/predicate/" %>% 
          dir(pattern="^Specimen*.*csv$", full.names = T) %>% 
          map_dfr(~{parse_dxh_predicate_data(.x)})
  
  pred = bind_rows(
                    audf,
                    dxdf,
                    acdf
                  )
  
  pred = 
    pred %>% 
    group_by(assay_name) %>% 
    group_modify(~{
      if (.x$assay_name[1] == "HGB") {
        .x$value = .x$value/10
      }
      if (.x$assay_name[1] == "hsCRP") {
        .x$value = .x$value * 1
      }
      if (.x$assay_name[1] == "NEUT%") {
        .x$value = .x$value * 100
      }
      if (.x$assay_name[1] == "LYMPH%") {
        .x$value = .x$value * 100
      }
      if (.x$assay_name[1] == "MONO%") {
        .x$value = .x$value * 100
      }
      if (.x$assay_name[1] == "EOS%") {
        .x$value = .x$value * 100
      }
      if (.x$assay_name[1] == "HCT") {
        .x$value = .x$value * 100
      }
      if (.x$assay_name[1] == "fT4") {
        .x$value = .x$value * 0.0777
      }
      return(.x %>% select(-assay_name))
    },.keep=T)
      
  return(pred)
}

compile_predicate_tstamps = function() {

 audf = "./data/predicate/" %>%
    dir(pattern="^AU)*.*csv$", full.names = T) %>%
    map_dfr(~{
      au_ts = .x %>% str_split_i("/",4) %>% str_split_i("\\.",1) %>% str_split("_")
      au_ts = parse_date_time(paste0(au_ts[[1]][2]," ",au_ts[[1]][3]), "ymd HM")
      au_ts = force_tz(au_ts, "America/New_York")
      return(data.frame(fname = .x, instr = "BC AU480", tstamp = au_ts))
    })
#
# 
 acdf = "./data/predicate/" %>% 
   dir(pattern="^05-*.*csv$", full.names = T) %>% 
   map_dfr(~{
      ac_ts = .x %>% str_split_i("/",4) %>% str_split_i("\\.",1)
      ac_ts = ac_ts %>% parse_date_time("m-d-y-H-M-S")
      ac_ts = force_tz(ac_ts, "America/New_York")
      return(data.frame(fname = .x, instr = "BC Access", tstamp = ac_ts))
   })
 
  dxdf = "./data/predicate/" %>%
    dir(pattern="^Specimen*.*csv$", full.names = T) %>% 
    map_dfr(~{
      dx_ts = .x %>% str_split_i("/",4) %>% str_split_i("\\.",1) %>% str_extract(" 500_2024-*.*m$") %>% str_replace("500_","")
      dx_ts = dx_ts %>% parse_date_time("y-m-d_HM")
      dx_ts = force_tz(dx_ts, "America/New_York")
      return(data.frame(fname = .x, instr = "BC DxH", tstamp = dx_ts))
    })
  
  return(bind_rows(audf, acdf, dxdf))
}


detect_outliers <- function(res) {

    res = res %>% 
    mutate(dj = (vital-predicate)/predicate) %>% 
    mutate(outlier=F) %>% 
    filter(!is.na(vital)) %>% 
    filter(!is.nan(dj))
  
  if (nrow(res) < 3) {return(res)}
  
  rosTest <- with(res, rosnerTest(dj, k = 2, alpha = 0.01, warn = T))
  nOutliers <- nrow(subset(rosTest$all.stats, Outlier))
  
  if (nOutliers > 0) {
    iOutliers <- subset(rosTest$all.stats, Outlier)$Obs.Num
    res[iOutliers ,"outlier"]=T
  }
  
  return(res)
}


recalibrate = function(assay, corr) {
  if (length(corr)==1) {
    cals =
      caldf %>% 
      filter( assay_name == assay) %>% 
      group_by(assay_name) %>% 
      group_modify(~{
        .x = .x %>% mutate(cal = corr * cal)
        mod = with(.x, lm(calOD ~ cal))
        return(data.frame(t(coef(mod))))
      },.keep=T) %>% 
      magrittr::set_names(c("assay_name","intercept","slope"))
    
    remcdf = 
      mcdf %>% 
      filter(assay_name == assay) %>% 
      mutate(vital = (OD - cals$intercept)/cals$slope) %>% 
      filter(vital > 0) %>% 
      detect_outliers()
  } else {
    
    cals =
      caldf %>% 
      filter( assay_name == assay) %>% 
      group_by(assay_name) %>% 
      group_modify(~{
        mod = with(.x, lm(calOD ~ cal))
        return(data.frame(t(coef(mod))))
      },.keep=T) %>% 
      magrittr::set_names(c("assay_name","intercept","slope"))
    
    remcdf = 
      mcdf %>% 
      filter(assay_name == assay) %>% 
      mutate(vital = (OD - cals$intercept)/cals$slope) %>% 
      mutate(vital = (vital * corr[2]) + corr[1]) %>% 
      filter(vital > 0) %>% 
      detect_outliers()
    if (nrow(remcdf)==0) {stop("all removed")}
  }
  
  pd = compile_tae_limits(ainfo %>% filter(meas == assay))
  lims = pd$lims
  tae  = pd$tae
  mdl = refr = pd$refr
  
  modDem  = with(filter(remcdf, !outlier), mcreg(predicate, vital, method.reg = "WDeming", mref.name="predicate", mtest.name="vital"))
  
  modifier = ifelse(length(corr)==1, ifelse(corr==1, "orig" , "reass"), "r-bcorr")
  
  make_vitalPlot_Hdemo(remcdf, assay, tae, lims, mdl, ainfo %>% filter(meas == assay), modDem, modifier)
  
  return(remcdf)
}


make_vitalPlot_Hdemo <- function(res, measurand, tae, lims, mdl, ainfo, mod, modifier, species="human") {
  
  
  pdf(file=paste0("./results/MethodComparison_Plots_", str_replace(measurand,"\\%",""),"_",species,"_type-", modifier, ".pdf"),  width=7, height=9)
 
  layout(matrix(c(1,2), ncol=1), heights=lcm(c(15.24,7.62)))
  
  if (TRUE) {
    
    par(mar = c(5, 7, 4.9, 5))
    par(family="Vital")
    
    if (ainfo$plot_axis_log) {
      with(res, plot(predicate, vital, col="white", xlim=lims, ylim=lims, log="xy", axes=F, bty="o", pty="s", asp=1,
                     ann=F))
    } else {
      with(res, plot(predicate, vital, col="white", xlim=lims, ylim=lims, axes=F, bty="o", pty="s", asp=1,
                     ann=F))
    }
    box(lwd=0.5, col=vpal["vslate"], bty="l")
    
    if (modifier %in% c(1,3)) {
      with(tae, lines(xr, ytae_max, col=vpal["vsunset"], lwd=1.5))
      with(tae, lines(xr, ytae_min, col=vpal["vsunset"], lwd=1.5))
    }
    
    mtext(side=1, text=paste0("Predicate Method, ",ainfo$units), cex.lab=1, col=vpal["vgreen"], line=1)
    mtext(side=2, text=paste0("Vital Method, ",ainfo$units), cex.lab=1, col=vpal["vgreen"], line=1.2)
    
    xt = axTicks(1)
    xlabels = sapply(xt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=1))))
    if (measurand == "RBC") {
      xlabels = sapply(xt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=2))))
    }
      
    axis(1, at=xt, labels=xlabels, cex.axis=0.7, lwd=0.5, col=vpal["vslate"], lwd.ticks=NA, mgp=c(0,0,0), col.axis=vpal["vgreen"])
    
    yt = xt
    ylabels = sapply(yt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=2))))
    if (measurand == "RBC") {
      ylabels = sapply(yt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=2))))
    }
    axis(2, at=yt, labels=ylabels, cex.axis=0.7, lwd=0.5, col=vpal["vslate"], lwd.ticks=NA, mgp=c(0,0.2,0), col.axis=vpal["vgreen"], las=1)
    
    abline(h=0, col="grey50", lwd=0.5)
    abline(h=yt, col="grey90", lwd=0.5)
    abline(v=xt, col="grey90", lwd=0.5)

    abline(v=c(ainfo$amr_low, ainfo$amr_high), lty=2, col=vpal["vslate_dark"], lwd=1)
    abline(v=mdl, lty=2, col=vpal["vocean"], lwd=1)

    
    #modRes = as.data.frame(calcResponse(mod, x.levels=seq(min(res$predicate), max(res$predicate), length.out=100)))
    xx = seq(min(res$predicate), max(res$predicate), length.out=100)
    modRes = predict(mod, newdata=data.frame(predicate = xx))
    
    if (modifier %in% c(1,2)) {
      lines(xx, modRes, col=vpal["vgreen"], lwd=1.5)
    } else {
      abline(0, 1, col=vpal["vslate"], lwd=0.5)
    }
    
    with(subset(res, uF_QC=="PASS"),
         points(predicate, vital, pch=21, bg=vpal["vleaf"], col=vpal["vgreen"], cex=1))
    with(subset(res, uF_QC=="FAIL"),
         points(predicate, vital, pch=17, col=vpal["vpoppy"], cex=0.9))
    with(subset(res, outlier),
         points(predicate, vital, pch=4, col=vpal["vpoppy"], cex=1.2))
  }
  
  if (TRUE) {
    par(mar = c(5, 4, 0, 1.5))
    par(family="Vital")
    
    if (ainfo$diffplot_type=="P") {
      if (ainfo$plot_axis_log) {
        with(res, plot(predicate, 100*(vital-predicate)/predicate, col="white", xlim=lims, ylim=c(-ainfo$plot_lim_err, ainfo$plot_lim_err), log="x", axes=F, bty="o",
                       ann=F))
      } else {
        with(res, plot(predicate, 100*(vital-predicate)/predicate, col="white", xlim=lims, ylim=c(-ainfo$plot_lim_err, ainfo$plot_lim_err), axes=F, bty="o",
                       ann=F))
      }
      
      
      box(lwd=0.5, col=vpal["vslate"], bty="l")
      
      mtext(side=1, text=paste0("Predicate Method, ",ainfo$units), cex.lab=1, col=vpal["vgreen"], line=1)
      mtext(side=2, text=paste0("Percent Error, (Y-X)/X"), cex.lab=1, col=vpal["vgreen"], line=1.2)
      
      xt = axTicks(1)
      xlabels = sapply(xt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=1))))
      axis(1, at=xt, labels=xlabels, cex.axis=0.7, lwd=0.5, col=vpal["vslate"], lwd.ticks=NA, mgp=c(0,0,0), col.axis=vpal["vgreen"])
      
      yt = 100*seq(-0.01*ainfo$plot_lim_err, 0.01*ainfo$plot_lim_err, length.out=5)
      ylabels = sapply(yt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=1))))
      axis(2, at=yt, labels=ylabels, cex.axis=0.7, lwd=0.5, col=vpal["vslate"], lwd.ticks=NA, mgp=c(0,0.2,0), col.axis=vpal["vgreen"], las=1)
      
      abline(h=0, col="grey50", lwd=0.5)
      abline(h=yt, col="grey90", lwd=0.5)
      abline(v=xt, col="grey90", lwd=0.5)

      abline(v=c(ainfo$amr_low, ainfo$amr_high), lty=2, col=vpal["vslate_dark"], lwd=1)      
      abline(v=mdl, lty=2, col=vpal["vocean"], lwd=1)

      with(tae, lines(xr, 100*atae_max/xr, col=vpal["vsunset"], lwd=1.5))
      with(tae, lines(xr, 100*atae_min/xr, col=vpal["vsunset"], lwd=1.5))
      
      with(subset(res, uF_QC=="PASS"),
           points(predicate, 100*(vital-predicate)/predicate, pch=21, bg=vpal["vleaf"], col=vpal["vgreen"], cex=1))
      with(subset(res, uF_QC=="FAIL"),
           points(predicate, 100*(vital-predicate)/predicate, pch=17, col=vpal["vpoppy"], cex=0.9))
      with(subset(res, outlier),
           points(predicate, 100*(vital-predicate)/predicate, pch=4, col=vpal["vpoppy"], cex=1.2))
    }
    
    if (ainfo$diffplot_type=="A") {
      if (grepl("EOS", measurand)) {ylims = c(-3.5, 3.5)} else {ylims = c(-1.5, 1.5)}
      if (ainfo$plot_axis_log) {
        with(res, plot(predicate, vital-predicate, col="white", xlim=lims, ylim=ylims, log="x", axes=F, bty="o",
                       ann=F))
      } else {
        with(res, plot(predicate, vital-predicate, col="white", xlim=lims, ylim=ylims, axes=F, bty="o",
                       ann=F))
      }
      
      
      box(lwd=0.5, col=vpal["vslate"], bty="l")
      
      mtext(side=1, text=paste0("Predicate Method, ",ainfo$units), cex.lab=1, col=vpal["vgreen"], line=1)
      mtext(side=2, text=paste0("Absolute Error, Y-X"), cex.lab=1, col=vpal["vgreen"], line=1.2)
      
      xt = axTicks(1)
      xlabels = sapply(xt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=1))))
      axis(1, at=xt, labels=xlabels, cex.axis=0.7, lwd=0.5, col=vpal["vslate"], lwd.ticks=NA, mgp=c(0,0,0), col.axis=vpal["vgreen"])
      
      yt = seq(-0.5, 0.5, length.out=5)
      ylabels = sapply(yt, function(x) ifelse(x %% 1 == 0, paste0(floor(x)), paste0(format(x, digits=1))))
      axis(2, at=yt, labels=ylabels, cex.axis=0.7, lwd=0.5, col=vpal["vslate"], lwd.ticks=NA, mgp=c(0,0.2,0), col.axis=vpal["vgreen"], las=1)
      
      abline(h=0, col="grey50", lwd=0.5)
      abline(h=yt, col="grey90", lwd=0.5)
      abline(v=xt, col="grey90", lwd=0.5)
      
      abline(v=c(ainfo$amr_low, ainfo$amr_high), lty=2, col=vpal["vslate_dark"], lwd=1)
      abline(v=mdl, lty=2, col=vpal["vocean"], lwd=1)
      
      with(tae, lines(xr, atae_max, col=vpal["vsunset"], lwd=1.5))
      with(tae, lines(xr, atae_min, col=vpal["vsunset"], lwd=1.5))
      
      with(subset(res, uF_QC=="PASS"),
           points(predicate, vital-predicate, pch=21, bg=vpal["vleaf"], col=vpal["vgreen"], cex=1))
      with(subset(res, uF_QC=="FAIL"),
           points(predicate, vital-predicate, pch=17, col=vpal["vpoppy"], cex=0.9))
      with(subset(res, outlier),
           points(predicate, vital-predicate, pch=4, col=vpal["vgreen"], cex=1.2))
    }
  }
  
  if (TRUE) {
    add_legend <- function(...) {
      opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                  mar=c(0, 0, 0, 0), new=TRUE)
      on.exit(par(opar))
      plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1))
      legend(...)
    }
    
    add_legend2 <- function(...) {
      opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                  mar=c(0, 0, 0, 0), new=TRUE)
      on.exit(par(opar))
      plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1))
      legend(...)
    }
    
    
    add_title <- function(...) {
      opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                  mar=c(0, 0, 0, 0), new=TRUE)
      on.exit(par(opar))
      plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1))
      text(...)
    }
    
    
    
    add_legend(x=0.48,y=1.018, legend=c("valid result", "discordant"),
               pch=c(21, 4),
               col=c(vpal[c("vgreen", "vpoppy")]),
               pt.bg=c(vpal[c("vleaf","vpoppy")]),
               pt.cex=c(0.9,0.9),
               lty = c(NA, NA, NA),
               lwd=c(NA, NA, NA),
               horiz=TRUE, bty='n', cex=0.75, text.col=vpal["vgreen"], xjust=0, x.intersp=0.7, seg.len=c(1,1,1),
               text.width=c(0.116, 0.112))
    # 
    add_legend2(x=0.48, y=0.988, legend=c(ifelse(modifier %in% c(1,2), "WDeming","y = x"),"TEa","MDL","AMR"),
                pch = NA,
                lty = c(1,1,11,11),
                lwd = c(1.25,1.25,1.25,1.25),
                col = c(vpal["vgreen"],vpal["vsunset"],vpal["vocean"], vpal["vslate_dark"]),
                horiz=TRUE, bty='n', cex=0.75, text.col=vpal["vgreen"], seg.len=c(1.5,1.5,1.5,1.5), xjust=0, x.intersp=0.2,
                text.width=c(0.08, 0.08, 0.08, 0.08))
    
    add_title(0.0045, 1.0, "Method Comparison Plots", pos=4, family="VitalBold", cex = 1.2, font=2, col=vpal["vgreen"])
    add_title(0.0045, 0.97, paste0(ainfo$full_name), pos=4, family="VitalBold", cex = 1.2, font=2, col=vpal["vgreen"])
    add_title(0.0045, 0.94, paste0(species), pos=4, family="VitalBold", cex = 1.2, font=2, col=vpal["vgreen"])
    
  }
  
  dev.off()
}


make_composite = function(df) {
  df = 
    df %>% 
    mutate(qc1 = ifelse(uF_QC_1=="PASS", T, F)) %>% 
    mutate(qc2 = ifelse(uF_QC_2=="PASS", T, F)) %>% 
    mutate(qc3 = ifelse(uF_QC_3=="PASS", T, F)) %>% 
    mutate(composite = ifelse(qc1 & is.na(qc2), vital_1,
                            ifelse(is.na(qc1) & qc2, vital_2,
                                ifelse(qc1 & qc2, vital_1, 
                                       ifelse(!qc1 & is.na(qc2), NA, 
                                              ifelse(qc1 & !qc2, vital_1,
                                                     ifelse(!qc1 & qc2, vital_2,
                                                            ifelse(is.na(qc1) & qc2, vital_2,
                                                                   ifelse(is.na(qc1) & !qc2, NA, 
                                                                          ifelse(!qc1 & !qc2, NA, NA)))))))))) %>% 
    mutate(composite = ifelse(is.na(composite), ifelse(!is.na(qc3), ifelse(qc3, vital_3, NA), NA), composite))
    
  
  return(df)
}

make_dtable = function(sample_name, rdf, plot_table = T) {
  df <- 
    rdf %>% 
    select(-OD) %>% 
    filter(sampleID == sample_name) 

  if (df %>% select(run_num) %>% distinct() %>% pull() %>% length() == 1) {
    df = df %>% bind_rows(data.frame(sampleID = sample_name, predicate=NA, vital=NA, assay_name='dummy', run_num=2, uF_QC=NA))
  }
  df = df %>% 
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
      palette = c(vpal["vpoppy"],"gray70"),
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
      container.height = px(600),
      container.padding.y = px(50),
      table.border.top.width = px(1),
      table.border.right.width = px(1),
      table.border.bottom.width = px(1),
      table.border.left.width = px(1)
    )
  
  
  if (plot_table) {dtable %>% gtsave(paste0("./results/table_", sample_name, ".pdf"))}
  
  return(df |> select(assay_name, predicate, composite, units, rounding))
  
}

make_dtable_single_sample = function(sample_name, rdf, assay_list, plot_table = T) {
  tdf <- 
    rdf %>% 
    filter(sampleID == sample_name) %>% 
    filter(assay_name %in% c(assay_list)) %>% 
    left_join(a_order, by="assay_name") %>% 
    arrange(sort_id) %>% 
    select(-sort_id, -value, -timestamp)
  
  if (tdf %>% select(run_num) %>% distinct() %>% pull() %>% length() == 1) {
    if (tdf$run_num[1] ==1) { rn1 = 2; rn2=3}
    if (tdf$run_num[1] ==2) { rn1 = 1; rn2=3}
    if (tdf$run_num[1] ==3) { rn1 = 2; rn2=1}
    
    tdf = tdf %>% bind_rows(tdf[1:2,])
    tdf$run_num[nrow(tdf)-1] = rn1;         tdf$run_num[nrow(tdf)] = rn2
    tdf$vital[nrow(tdf)-1] = NA;            tdf$vital[nrow(tdf)] = NA
    tdf$predicate[nrow(tdf)-1] = NA;        tdf$predicate[nrow(tdf)] = NA
    tdf$assay_name[nrow(tdf)-1] = "dummy";  tdf$assay_name[nrow(tdf)] = "dummy"
    tdf$uF_QC[nrow(tdf)-1] = NA;            tdf$uF_QC[nrow(tdf)] = NA
  }
  
  if (tdf %>% select(run_num) %>% distinct() %>% pull() %>% length() == 2) {
    if (all(unique(tdf$run_num) %in% c(1,2))) { rn1 = 3}
    if (all(unique(tdf$run_num) %in% c(1,3))) { rn1 = 2}
    if (all(unique(tdf$run_num) %in% c(2,3))) { rn1 = 1}
    
    tdf = tdf %>% bind_rows(tdf[1,])
    tdf$run_num[nrow(tdf)] = rn1
    tdf$vital[nrow(tdf)] = NA
    tdf$predicate[nrow(tdf)] = NA
    tdf$assay_name[nrow(tdf)] = "dummy"
    tdf$uF_QC[nrow(tdf)] = NA
  }
  
  
  tdf = tdf %>% 
    select(-run_id, -run_name, -module, -structure_label) %>% 
    mutate(vital = round(vital, 2)) %>% 
    mutate(perErr = round(100*(vital-predicate)/predicate,1)) %>% 
    mutate(absErr = vital-predicate) %>% 
    left_join(select(ainfo, meas, rounding, units, tae_per, tae_abs) , by=join_by(assay_name == meas)) %>% 
#    mutate(within_TEa = ifelse(abs(absErr) <= tae_abs | abs(perErr) <= tae_per, T, F)) %>% 
    select(-perErr, -absErr, -tae_per, -tae_abs) %>% 
    distinct() %>% 
    pivot_wider(names_from = run_num, values_from = c(vital,uF_QC)) %>% 
    make_composite() %>% 
    #mutate(within_TEa_1 = ifelse(is.na(within_TEa_1), T, within_TEa_1)) %>% 
    #mutate(within_TEa_2 = ifelse(is.na(within_TEa_2), T, within_TEa_2)) %>% 
    mutate(predicate = round(predicate, rounding)) %>% 
    mutate(vital_1 = round(vital_1, rounding)) %>% 
    mutate(vital_2 = round(vital_2, rounding)) %>% 
    mutate(vital_3 = round(vital_3, rounding)) %>% 
    mutate(composite = round(composite, rounding)) %>% 
    filter(!is.na(predicate)) %>% 
    rowwise() %>% 
    mutate(
      vital_2 = ifelse(!is.na(uF_QC_1), ifelse(uF_QC_1 == "PASS", NA, vital_2), vital_2)
    ) %>% 
    mutate(
      vital_3 = ifelse(!is.na(uF_QC_2) & !is.na(uF_QC_2), ifelse(uF_QC_2 == "PASS", NA, vital_3), vital_3)
    )
  
  
  dtable = 
    tdf |>
    select(assay_name, predicate, vital_1, uF_QC_1, vital_2, uF_QC_2, vital_3, uF_QC_3, composite, units, rounding) |>
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
    cols_hide(columns = c(rounding)) |> 
    cols_align(
      align = "center"
    )|> 
    tab_header(
      title = paste0("Results for ", sample_name),
      subtitle = paste0("L: archived serum, E: fresh K2-EDTA whole blood | Analyzed: ", tdf$timestamp[1])
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
      vital_3 = "Value",
      uF_QC_1 = "QC flag",
      uF_QC_2 = "QC flag",
      uF_QC_3 = "QC flag",
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
      label = "Vital Run 3",
      columns = c("vital_3","uF_QC_3"),
      id = "num_spanner_22"
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
    )|>
    tab_footnote(
      footnote = "Privileged and Confidential",
      placement = "right"
    )
  
  if (all(is.na(tdf$uF_QC_1))) {
    dtable = 
      dtable|>
      data_color(
        columns = c("uF_QC_2", "uF_QC_1", "uF_QC_3"), 
        target_columns = c("uF_QC_2", "uF_QC_1", "uF_QC_3"),
        method = "factor",
        palette = c("gray70",vpal["vpoppy"]),
        ordered = T,
        apply_to = "text"
      )
  } else {
    dtable = 
      dtable|>
      data_color(
        columns = c("uF_QC_1", "uF_QC_2", "uF_QC_3"), 
        target_columns = c("uF_QC_1", "uF_QC_2", "uF_QC_3"),
        method = "factor",
        palette = c("gray70",vpal["vpoppy"]),
        ordered = T,
        apply_to = "text"
      )
  }
  
  
  
  dtable %>% gtsave(paste0("./results/table_", sample_name, ".pdf"))
  
  return(tdf |> select(assay_name, predicate, composite, units, rounding))
  
}

note_removals = function(test_run_id, orub, fdf) {
  
  ffdf = 
    fdf %>% 
    filter(run_id == test_run_id) 
  
  tvec = ffdf %>% select(matches("rbc")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(13,14,15,16,17,18,19)] = TRUE
    }
  }
  
  tvec = ffdf %>% select(starts_with("hb")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(12)] = TRUE
    }
  }
  
  tvec = ffdf %>% select(starts_with("wbc")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(20:26)] = TRUE
    }
  }
  
  tvec = ffdf %>% select(starts_with("strA")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(30)] = TRUE
    }
  }
  
  tvec = ffdf %>% select(starts_with("strB")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(27)] = TRUE
    }
  }
  
  tvec = ffdf %>% select(starts_with("strC")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(28)] = TRUE
    }
  }
  tvec = ffdf %>% select(starts_with("strD")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(29)] = TRUE
    }
  }
  
  
  # @Luis Barbosa @Amrita On CC, strA = 1:20, strC=1:10, strD = 1:50. Correct?
  tvec = ffdf %>% select(starts_with("str10x")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(4,5,6)] = TRUE
    }
  }
  tvec = ffdf %>% select(starts_with("str20x")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec == "n")) {
      orub$uF_failure[c(1,2,3)] = TRUE
    }
  }
  tvec = ffdf %>% select(starts_with("str50x")) %>% as.character() 
  if (!all(is.na(tvec))) {
    tvec = tvec[!is.na(tvec)]
    if (any(tvec[-5]== "n")) {
      orub$uF_failure[c(7:11)] = TRUE
    }
  }
  
  return(orub)
}

make_plots_final = function(cdf) {
  
  dum = 
    cdf %>% 
    group_map(~{
      assay = .x$assay[1]
      species = .x$species[1]
      pd    = compile_tae_limits(ainfo %>% filter(meas == assay))
      lims  = pd$lims
      tae   = pd$tae
      mdl   = refr = pd$refr
      
      .x = .x %>% 
            filter(!is.na(vital)) %>% 
            filter(!is.na(predicate))
      if (nrow(.x) == 0) {return()}
      if (nrow(.x)<2) {
          mod  = with(filter(.x), lm(seq(10) ~ seq(10)))
      } else {
          mod  = with(filter(.x), lm(vital ~ predicate))
      }
      
      for (modifier in c(3)) {
        print(modifier)
        make_vitalPlot_Hdemo(.x, assay, tae, lims, mdl, ainfo %>% filter(meas == assay), mod, modifier, species)
      }
      return()
    }, .keep=T)
  
  return()
}
  
  
make_summary_table = function(sample_list, assay_list) {
  sdf =
    sample_list %>%
    map_dfr(~{
      print(.x)
      .y = make_dtable_single_sample(.x, rdf, assay_list, F)
      .y = mutate(.y, sampleID = .x)
    })
  
  sdf = 
    sdf %>% 
    filter(assay_name %in% assay_list) %>% 
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
      label = "L1E1",
      columns = matches("L1_E1$"),
      id = "num_spanner_1"
    )|>
    tab_spanner(
      label = "L2E2",
      columns = matches("L2_E2$"),
      id = "num_spanner_2"
    )|>
    tab_spanner(
      label = "L6E4",
      columns = matches("L3_E3$"),
      id = "num_spanner_3"
    )|>
    tab_spanner(
      label = "L6_E4",
      columns = matches("L6_E4$"),
      id = "num_spanner_4"
    )|>
    tab_spanner(
      label = "L7_E5",
      columns = matches("L7_E5$"),
      id = "num_spanner_5"
    )|>
    tab_spanner(
      label = "L9_E6",
      columns = matches("L9_E6$"),
      id = "num_spanner_6"
    )|>
    tab_spanner(
      label = "L10_E7",
      columns = matches("L10_E7$"),
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
    ) |>
    tab_footnote(
      footnote = "Privileged and Confidential",
      placement = "right"
    )
  
  
  
  dtable %>% gtsave(paste0("./results/table_summary.html"))
  
  return(dtable)
}

make_sa_log = function(rdf, pred) {
  vdf = 
    rdf %>% 
    select(run_id, sampleID, timestamp, run_num) %>% 
    separate(sampleID, into=c("sampleLH","sampleED"), sep="_") %>% 
    pivot_longer(matches("sample")) %>% 
    mutate(instr = "VitalOne") %>% 
    mutate(timestamp = as.character(timestamp)) %>% 
    distinct() %>% 
    pivot_wider(id_cols=value, names_from=run_num, values_from=timestamp) %>% 
    magrittr::set_names(c("sampleID","VT1-Run1","VT1-Run2", "VT1-Run3"))
  
  
  
  bdf = 
    pred %>% 
    ungroup() %>% 
    filter(!is.na(sampleID)) %>% 
    select(sampleID, instr, tstamp) %>% 
    mutate(tstamp = as.character(tstamp)) %>% 
    distinct() %>% 
    pivot_wider(id_cols=sampleID, names_from=instr, values_from=tstamp) 
  
  vdf = 
    vdf %>% 
    left_join(bdf, by="sampleID")
  
  
  vtable =
    vdf |>
    gt()|>
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
    cols_align(
      align = "center"
    )|> 
    tab_header(
      title = paste0("Log for time of sample analysis"),
      subtitle = paste0(" ")
    )|>
    sub_values(
      fn = function(x) is.na(x), 
      replacement=""
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
    )|>
    tab_footnote(
      footnote = "Privileged and Confidential",
      placement = "right"
    )
  
  vtable %>% gtsave(paste0("./results/table_sample_analysis.pdf"))
  
}

make_bias_table = function(bias_table) {

vtable =
  bias_table |>
  gt()|>
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
  cols_align(
    align = "center"
  )|> 
  tab_header(
    title = paste0("Summary Table for Method Comp Stats"),
    subtitle = paste0(" ")
  )|>
  sub_values(
    fn = function(x) is.na(x), 
    replacement=""
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
  )|>
  tab_footnote(
    footnote = "Privileged and Confidential",
    placement = "right"
  )
  table_file  = paste0("./results/table_sample_analysis.html")
  vtable %>% gtsave(table_file)

  # html <- as_raw_html(vtable)
  # html <- gsub("overflow-y: scroll;", "", html)
  # html <- gsub("overflow-x: scroll;", "", html)
  # writeLines(html, "table.html")
}