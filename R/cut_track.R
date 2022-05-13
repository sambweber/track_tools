cut_track = function(data,datetime,tmax = NULL,cut.dates = NULL){
    
      mutate(data,t = traipse::track_time(!!sym(datetime))) %>%
      mutate(t = replace_na(t,0)) %>%
      mutate(segment = 1+cumsum(t > tmax*3600))
    
  }
  
