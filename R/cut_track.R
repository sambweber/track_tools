cut_track = function(data,datetime,tmax = NULL,cut.dates = NULL){
    
    datetime <- try(as.POSIXct(datetime), silent = TRUE)
    if (is(datetime,"try-error")) stop("Cannot convert 'datetime' to POSIXct")
    
    if(is.null(tmax) & is.null(cut.dates)) stop('Either tmax or cut.dates must be specified')
    if(!is.null(tmax) & !is.null(cut.dates)) stop("'tmax' and 'cut.dates' both provided! Splitting by tmax.")
    
    if(!is.null(tmax)){   
        
        mutate(data,t = traipse::track_time(!!sym(datetime))) %>%
        mutate(t = replace_na(t,0)) %>%
        mutate(segment = 1+cumsum(t > tmax*3600)) %>%
        return()
        
    } else {
       
        cut.dates = try(as.POSIXct(cut.dates),silent=TRUE)
        if(is(cut.dates('try-error')) stop ("Cannot convert 'cut.dates' to POSIXct")
        rowwise(data) %>% mutate(segment = 1+sum(datetime > cut.dates)) %>% 
        ungroup() %>% return()
    
  }
        
}  
