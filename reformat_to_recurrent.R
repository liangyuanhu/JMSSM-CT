# Function to reformat treatment initiation as recurrent events
reformat_to_recurrent <- function(data, treatment){
  
  data_treatment_1 <- data %>% 
    filter_(paste(treatment,"%in%",1)) %>% 
    arrange(id, time) %>% 
    group_by(id) %>%
    filter(row_number()==1) %>% 
    ungroup
  
  data_treatment_recurrent <- data %>% 
    filter_(paste(treatment,"%in%",0)) %>% 
    bind_rows(data_treatment_1) %>% 
    arrange(id, time) %>% 
    group_by(id) %>% 
    mutate(time_lag = lag(time, 1),
           time_lag = ifelse(is.na(time_lag), -1, time_lag),
           A1_lag = lag(A1,1, default = 0),
           A2_lag = lag(A2,1, default = 0)) %>% 
    select(
      id,
      start_time = time_lag,
      end_time = time,
      lag_A1 =  A1_lag,
      lag_A2 = A2_lag,
      everything()
    ) %>% 
    ungroup %>% as.data.frame()
}


