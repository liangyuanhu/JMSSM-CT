# Wrapper function to implement weight estimation using Cox model

weight_Cox_treatment <- function(data, treatment, time_varying_cov, smooth){
  simulation_data_ragged_treatment_recurrent_format <- reformat_to_recurrent(data = data, treatment)
  
  if (treatment == "A1"){
    weight_treatment_weight_Cox_TV <- weight_Cox_TV(exposure = A1,
                                                    denominator = ~ L1+L2+lag_A1 + lag_A2,
                                                    numerator = ~lag_A1 + lag_A2,
                                                    id = id,
                                                    tstart = start_time ,
                                                    timevar = end_time,
                                                    smooth = smooth,
                                                    data = simulation_data_ragged_treatment_recurrent_format)
  } 
  if (treatment == "A2"){
    weight_treatment_weight_Cox_TV <- weight_Cox_TV(exposure = A2,
                                                    denominator = ~ L1+L2+lag_A1 + lag_A2,
                                                    numerator = ~lag_A1 + lag_A2,
                                                    id = id,
                                                    tstart = start_time ,
                                                    timevar = end_time,
                                                    smooth = smooth,
                                                    data = simulation_data_ragged_treatment_recurrent_format)
  }
  weight_treatment_weight_Cox_TV_result <- simulation_data_ragged_treatment_recurrent_format %>% 
    select(id) %>% 
    mutate(w.numerator = weight_treatment_weight_Cox_TV$w.numerator,
           w.denominator = weight_treatment_weight_Cox_TV$w.denominator) %>% 
    group_by(id) %>% 
    summarise(w.numerator = prod(w.numerator),
              w.denominator = prod(w.denominator)) %>% 
    ungroup %>% 
    mutate(weight = w.numerator/w.denominator) %>% 
    select(id, weight_treatment = weight)
  return(weight_treatment_weight_Cox_TV_result$weight_treatment)
}
