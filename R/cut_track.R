cut_track = function(data,datetime,tmax = NULL,cut.dates = NULL){
    
      mutate(data,t = traipse::track_time(!!sym(datetime))) %>%
      mutate(t = replace_na(t,0)) %>%
      mutate(segment = 1+cumsum(t > tmax*3600))
     
    # date supplied method
     cut.dates = as.POSIXct(cut.dates)
     rowwise(data) %>% mutate(segment = 1+sum(date > cut.dates)) %>% ungroup()
    
  }
  
