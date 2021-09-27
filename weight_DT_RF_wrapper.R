# Wrapper function to implement weight estimation using disctete time

weight_JMSM_DT <- function(data, treatment, interval, time_varying_cov){
  data_treatment_cut_interval <- data %>% 
    filter_(paste(treatment,"%in%",1)) %>% 
    arrange(id, time) %>% 
    group_by(id) %>%
    filter(row_number()==1) %>% 
    ungroup %>% 
    mutate(t_start = plyr::round_any(time, interval, floor),
           t_end = plyr::round_any((time+0.01), interval, ceiling)) %>% 
    group_by(id,t_start,  t_end) %>% 
    filter(row_number()==n()) %>% 
    ungroup %>% rename(A = treatment) %>% 
    select(id, t0 = time, Y, A, L1, L2 , t_start, t_end)
  
  data_treatment_cut_interval_demo <- data_treatment_cut_interval %>% 
    group_by(id) %>% 
    filter(row_number()==n()) %>% 
    ungroup %>% 
    select(id,t0, L1, L2,salv = A, salvtime = t0, event = Y, survtime = t_end)
  
  data_treatment_1_dependent <-  data %>% 
    filter_(paste(treatment,"%in%",1)) %>% 
    arrange(id, time) %>% 
    group_by(id) %>%
    filter(row_number()==1) %>% 
    ungroup %>% 
    select(id, time, L = L1)
  
  weight_JMSM_DT_treatment <-
    DTR_Surv(
      sel_b = 1,
      sim_demo = data_treatment_cut_interval_demo,
      sim_pre_psa = data_treatment_1_dependent,
      pre_td_orig = data_treatment_cut_interval,
      wtTrunc = T,
      interval = interval
    )
  return(weight_JMSM_DT_treatment)
}