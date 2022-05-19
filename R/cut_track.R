# ------------------------------------------------------------------------------------------------------------------
# cut_track
# ------------------------------------------------------------------------------------------------------------------

#' Cuts a movement track into segments based on either a maximum allowable time step (tmax) between successive 
#' fixes or user specified cutting dates.

#' @param datetime A datetime vector of class \code{POSIXt*} containing the timestamps for each location.
#' @param tmax The maximum allowable time step between locations. Tracks will be cut at gaps > tmax.
#' @param cut.dates A vector of user defined datetimes at which tracks should be split. Must either be of class \code{POSIXt*} 
#' or coercible to \code{POSIXt*}
#' @param overlap A logical indicating whether track segments split by cut.dates should be overlapped by one location at the 
#' start and ends. This can be useful to ensure continuous predictions from state space models with no gaps.
#' 
#' @return The original data object with a new 'segment' column added containing the segment ids. 

# -------------------------------------------------------------------------------------------------------


cut_track = function(data,datetime,tmax = NULL,cut.dates = NULL,overlap=FALSE){
  
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
  

# ------------------------------------------------------------------------------------------------------------------
# split_bouts
# ------------------------------------------------------------------------------------------------------------------

#' Seperates a movement track into bouts of consisent behaviour based on a column of state labels. Optionally allows you to remove
#' short bouts of behaviour within larger more consistent blocks using a simple linear interpolation algorithm.

#' @param object A dataframe containing a column of state labels

#' @param state A character vector of length 1 giving the name of the column containing state labels

#' @param dt A character vector of length 1 giving the name of the column containing datetimes (must be in `POSIXt` format). Only 
#' required if `t.min` is specified.

#' @param t.min A numeric of length 1 giving the minimum allowable duration (in hours) for a bout. Bouts that are shorter than 
#' t.min will be assimilated into adjacent bouts by linear interpolation of the state labels implemented using the `zoo::na.fill` function.
#' Note: this won't work well and will give unexpected results if state labels are very changeable over individual timesteps 
#' (which might indicate the model has too many states). It is mainly intended to smooth out 'anomalies' in otherwise consitent blocks
#' of behaviour.
#' 
#' @return The original movement track with a 'bout' column added containing the numerical index of the bout

# -----------------------------------------------------------------------------------------------------------------

split_bouts = function(object, state, dt = NULL, t.min = NULL){
  
  if(!is.null(t.min) & is.null(dt)) stop ('Time column must be specified if t.min is provided')
  if(!state %in% names(object)) stop(cat("Column '",state,"' does not exist"))
  
  object[[state]] <- as.factor(object[[state]])
  labels = levels(object[[state]])
  
  st = sym(state)
  object = mutate(object,bout = 1+cumsum(!!st!=lag(!!st,default = first(!!st)))) 
  
  if(!is.null(t.min)){
    
      object = group_by(object,bout,.add=T) %>%
      mutate(dur = difftime(max(!!sym(dt)),min(!!sym(dt)),units='hours')) %>%
      ungroup(bout) %>% 
      mutate(!!st := ifelse(dur < t.min, NA, as.numeric(!!st))) %>%
      mutate(!!st := labels[round(zoo::na.fill(!!st,'extend'))]) %>%
      dplyr::select(-dur)
      
      object = Recall(object,state=state)

  }

  return(object)
}

