#' ------------------------------------------------------------------------------------------------------------------
#' cut_track
#' ------------------------------------------------------------------------------------------------------------------



cut_track = function(data,datetime,tmax = NULL,cut.dates = NULL,overlap=F){
  
  if (!is(data[[datetime]],"POSIXt")) stop("'datetime' should be of class POSIXt")
  
  if(is.null(tmax) & is.null(cut.dates)) stop('Either tmax or cut.dates must be specified')
  if(!is.null(tmax) & !is.null(cut.dates)) stop("'tmax' and 'cut.dates' both provided! Splitting by tmax.")
  
  if(!is.null(tmax)){   
    
    mutate(data,t = traipse::track_time(!!sym(datetime))) %>%
    mutate(t = replace_na(t,0)) %>%
    mutate(segment = 1+cumsum(t > tmax*3600)) %>%
    dplyr::select(-t) %>% return()
    
  } else {
    
    cut.dates = try(as.POSIXct(cut.dates,tz='GMT'),silent=TRUE)
    if(is(cut.dates,'try-error')) stop ("Cannot convert 'cut.dates' to POSIXct")
    data = rowwise(data) %>% mutate(segment = 1+sum(!!sym(datetime) > cut.dates)) 
    
    if(overlap){
      data = group_by(data,segment) %>% slice(1,n()) %>% 
      mutate(segment = segment+c(-1,1)) %>% 
      subset(!segment %in% c(0,length(cut.dates)+2)) %>%
      bind_rows(.,data) %>% arrange(!!sym(datetime))
    }
      
    return(ungroup(data))
  }
}
  
