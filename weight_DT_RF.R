# Function to implement weight estimation using disctete time (Random Forest)

DTR_Surv <- function(sel_b, sim_demo=sim_demo, sim_pre_psa=sim_pre_psa, pre_td_orig= pre_td_orig, wtTrunc=F, interval = interval)
{
  ### find Tsalv of regime sel_b for all patinets
  all_ids <- unique(sim_pre_psa$id)
  c_tsalv2 <- NULL # subject i's trt time for regime 1,causal framework 
  for (i in all_ids)
  {
    tmp_sim <- sim_pre_psa[sim_pre_psa$id==i, ] 
    tmp_sim_sorted <- tmp_sim[order(tmp_sim$time),] 
    tmp_time <- tmp_sim_sorted$time
    tmp_psa <- tmp_sim_sorted$psa
    tmp_l <- length(tmp_psa)
    tmp_psa_sl <- tmp_psa[-1] - tmp_psa[-tmp_l] 
    tmp_psa_sl <- c(-1, tmp_psa_sl)
    obsperyear <- 1 
    k <- 50 
    times2 <- 1:50
    pos_sel <- which((tmp_psa > sel_b) & (tmp_psa_sl > 0))[1] #find the position where PSA first go above b
    if (is.na(pos_sel)) {
      c_tsalv2 <- c(c_tsalv2, NA)
    } else {
      if (pos_sel < 2) {
        pos_sel = 2 # do not allow treatment at time 0 
      }
      c_tsalv2 <- c(c_tsalv2, tmp_time[pos_sel]) 
    }
  }
  ### Create time-dependent dataset (compliant to regime sel_b) 
  sim_regb <- sim_demo
  # for those who start treatment after their decision time for regime 2, event=0, survtime=time to decision point for regime 2.
  
  crt20 <- (sim_demo$salvtime > c_tsalv2 & sim_demo$survtime >= c_tsalv2 & (!is.na(c_tsalv2)) & (!is.na(sim_demo$salvtime)))
  
  sim_regb$survtime[crt20] <- c_tsalv2[crt20] 
  sim_regb$event[crt20] <- 0  
  
  sim_regb$censor[crt20] <- 1
  sim_regb$salvtime[crt20] <- NA
  # for those who start treatment before their decision time for regime 2, event=0, survtime=salvtime
  crt21 <- (sim_demo$salvtime < c_tsalv2 & (!is.na(c_tsalv2)) & (!is.na(sim_demo$salvtime)))
  sim_regb$survtime[crt21] = sim_demo$salvtime[crt21]
  sim_regb$event[crt21] <- 0
  sim_regb$censor[crt21] <- 1
  sim_regb$salvtime[crt21] <- NA
  # for those who start treatment at their decision time for regime 2, event=event,survtime=survtime
  crt22 <- ((sim_demo$salvtime == c_tsalv2) & (!is.na(c_tsalv2)) &
              (!is.na(sim_demo$salvtime)))
  sim_regb$survtime[crt22] <- sim_demo$survtime[crt22]
  sim_regb$event[crt22] <- sim_demo$event[crt22]
  sim_regb$censor[crt22] <- 0
  # for those who without decision time for regime 2, if they actually have no treatment,event=event, survtime=survtime.
  crt23 <- (is.na(c_tsalv2) & is.na(sim_demo$salvtime))
  # sim_reg2$survtime[crt23] <- sim_demo$survtime[crt23]
  sim_regb$event[crt23] <- sim_demo$event[crt23]
  sim_regb$censor[crt23] <- 0
  # for those who without decision time for regime 2, if they have a treatment, event=0,survtime=salvtime.
  crt24 <- (is.na(c_tsalv2) & (!is.na(sim_demo$salvtime)))
  sim_regb$survtime[crt24] <-sim_demo$salvtime[crt24]
  sim_regb$event[crt24] <- 0
  sim_regb$censor[crt24] <- 1
  sim_regb$salvtime[crt24] <- NA
  # for those who without treatment, if they die after their decision time for regime 2, event=0, survtime=time to decision point for regime 2.
  crt25 <- (c_tsalv2<sim_demo$survtime & (!is.na(c_tsalv2)) & is.na(sim_demo$salvtime)) 
  sim_regb$survtime[crt25] <- c_tsalv2[crt25]
  sim_regb$event[crt25] <- 0
  sim_regb$censor[crt25] <- 1
  # for those who without treatment, if they die before their decision time for regime 2,event=event(1), survtime=survtime
  crt26 <- (c_tsalv2>=sim_demo$survtime & (!is.na(c_tsalv2)) & is.na(sim_demo$salvtime)) #sim_regb$survtime[crt26] <- sim_demo$survtime[crt26]
  sim_regb$event[crt26] <- 1
  sim_regb$censor[crt26] <- 0
  ### delete the one that reach decision point 2 at the very beginning 
  sim_regb <- sim_regb[sim_regb$survtime>0,]
  ###### create time dependent dataset
  td_sel_reg <- apply(sim_regb,2,rep,times=ceiling(sim_regb$survtime * obsperyear - 1/obsperyear/2))
  td_sel_reg <- as.data.frame(td_sel_reg, stringsAsFactors=F) 
  td_sel_reg$id <- as.numeric(td_sel_reg$id) 
  td_sel_reg$censor <- as.numeric(td_sel_reg$censor) 
  td_sel_reg$survtime <- as.numeric(td_sel_reg$survtime) 
  td_sel_reg$survtime <- round(td_sel_reg$survtime,2) 
  td_sel_reg$salvtime <- as.numeric(td_sel_reg$salvtime) 
  td_sel_reg$salvtime <- round(td_sel_reg$salvtime,2)
  t_end <- unlist(sapply(ceiling(sim_regb$survtime*obsperyear-0.1),function(len) {
    seq(from=(1/obsperyear),by=(1/obsperyear),length.out=len) }))
  t_start <- unlist(sapply(ceiling(sim_regb$survtime*obsperyear-0.1),function(len) 
  {
    seq(from=0,by=(1/obsperyear),length.out=len)
  })) 
  td_sel_reg$t_start <- round(t_start,2)[1:dim(td_sel_reg)[1]];
  td_sel_reg$t_end <- round(t_end,2)[1:dim(td_sel_reg)[1]]
  td_dup <- td_sel_reg[td_sel_reg$censor==1 & td_sel_reg$t_end >= td_sel_reg$survtime, ] 
  td_dup$t_start <- round(td_dup$t_start + (1/obsperyear), 2)
  td_dup$t_end <- round(td_dup$t_end + (1/obsperyear), 2)
  td_sel_reg <- rbind(td_sel_reg, td_dup)
  td_sel_reg <- td_sel_reg[with(td_sel_reg, order(id, t_start)), ]
  td_ind_censor0 <- as.numeric(td_sel_reg$censor==1 & td_sel_reg$t_start>=td_sel_reg$survtime)
  td_ind_event0 <- as.numeric(td_sel_reg$event==1 & td_sel_reg$t_end>=td_sel_reg$survtime)
  td_sel_reg$td_ind_event <- td_ind_event0 
  td_sel_reg$td_ind_censor <- td_ind_censor0
  ### calculate denominator of the weight
  pre_td_orig$intx1 <- pre_td_orig$t_start * pre_td_orig$L2 
  pre_td_orig$intx2 <- pre_td_orig$t_start * pre_td_orig$L1 
  pre_td_orig$intx3 <- pre_td_orig$t_start * pre_td_orig$L2 
  
  ## using random forest with oob
  mod02_den <- tuneRF(x=pre_td_orig[,c("t_start", "L1", "L2", "intx1", "intx2", "intx3",
                                       "t_start", "L1", "L2", "intx1", "intx2", "intx3")],
                      y=(!pre_td_orig$t0), mtryStart=1, ntreeTry=40000/interval, stepFactor=2,
                      doBest=T)
  # print(mod02_den)
  mod02_fit = predict(mod02_den)
  noad <- (is.na(pre_td_orig$t0) | (pre_td_orig$t_start <= pre_td_orig$t0))
  mod02_fit[mod02_fit > 1] <- 1 ## correct for possible overflow by Random Forest due to machine accuracy
  mod02_fit[pre_td_orig$t0 == 1] <- 1 - mod02_fit[pre_td_orig$t0 == 1]
  fit02 <- rep(1, dim(pre_td_orig)[1])
  fit02[noad] <- mod02_fit
  ids2 <- unique(pre_td_orig$id)
  ipcw_den <- NULL 
  for (i in ids2) {
    tmp <- pre_td_orig[pre_td_orig$id==i,] 
    mod02_tmp <- fit02[pre_td_orig$id==i] 
    ipcw_d <- mod02_tmp
    tmp_d <- cumprod(ipcw_d)
    ipcw_den <- c(ipcw_den, tmp_d) }
  pre_td_orig$wt <- 1/ipcw_den # cal denominator first, time varying
  ### calculate numerator of the weight
  pre_td_reg <- td_sel_reg[(is.na(td_sel_reg$td_ind_event) | (td_sel_reg$t_start<=td_sel_reg$td_ind_event)),]
  pre_td_reg$intx1 <- pre_td_reg$t_start * pre_td_reg$L2
  ## using random forest with oob
  
  (!pre_td_reg$td_ind_censor) %>% length
  mod01_num <- tuneRF(x=pre_td_orig[,c("t_start", "L1", "L2", "intx1",
                                       "t_start", "L1", "L2", "intx1")],
                      y=(!pre_td_orig$t0), mtryStart=1, ntreeTry=40000/interval, stepFactor=2,
                      doBest=T)
  # print(mod01_num)
  mod01_fit = predict(mod01_num)
  mod01_fit[mod01_fit > 1] <- 1 ## correct for possible overflow by Random Forest due to machine accuracy
  notrt <- (is.na(td_sel_reg$salvtime) | (td_sel_reg$t_start<=td_sel_reg$salvtime)) 
  fit01 <- rep(1, dim(td_sel_reg)[1])
  fit01[notrt] <- mod01_fit
  noad <- (td_sel_reg$survtime >= td_sel_reg$t_end)
  fit01 <- fit01[noad]
  td_sel_reg <- td_sel_reg[noad,] ## take out the piece after censoring
  reg_ids <- unique(td_sel_reg$id) 
  ipcw_num <- NULL
  wt_n <- rnorm(1000, 1, 0.1)
  for (i in reg_ids) {
    tmp <- td_sel_reg[td_sel_reg$id==i,] 
    mod01_tmp <- fit01[td_sel_reg$id==i] 
    ipcw_n <- mod01_tmp
    tmp_n <- cumprod(ipcw_n)
    ipcw_num <- c(ipcw_num, tmp_n) }
  td_survtime <- numeric(dim(pre_td_orig)[1]) 
  for (j in reg_ids)
  {
    td_survtime[pre_td_orig$id == j] <- sim_regb$survtime[sim_regb$id == j] 
  }
  idx_reg_compliant <- (pre_td_orig$t_start < td_survtime) ## index of time pieces of td_sel_reg within pre_td_orig
  td_sel_reg$wt <- pre_td_orig$wt[idx_reg_compliant] * ipcw_num
  range(td_sel_reg$wt)
  ## weight truncation could set as an option 
  if (wtTrunc) {
    td_sel_reg$wt[td_sel_reg$wt < 0.05] <- 0.05
    td_sel_reg$wt[td_sel_reg$wt > 20] <- 20 }
  ## estimate the regime specific survival function (weighted Nelson-Aalen estimator)
  wtd_haz <- NULL
  nowt_haz <- NULL
  for (i in times2) {
    wtd_n <- sum(td_sel_reg$wt[td_sel_reg$t_end==i])
    wtd_d <- sum(td_sel_reg$wt[td_sel_reg$t_end==i & td_sel_reg$L1==1]) 
    wtd_haz <- c(wtd_haz, (wtd_d/wtd_n))
    nowt_n <- sum(td_sel_reg$t_end==i)
    nowt_d <- sum(td_sel_reg$t_end==i & td_sel_reg$L1==1)
    nowt_haz <- c(nowt_haz, (nowt_d/nowt_n))
  }
  return(wt_n)
}