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



# -------------------------------------------------------------------------------------------------------------------------
# split_trips: function for splitting movement track of a central place forager into discrete trips 
# -------------------------------------------------------------------------------------------------------------------------

#' Trips are defined based on movements outside a distance buffer placed around the central place (colony).

#' @params time A POSIXt vector containing times of each track location.
#' @params X,Y A vector of latitudes and longtidues for the track locations
#' @params colony.X,colony.Y A vector of latitudes and longitudes or a single latitude/longitude giving the 
#' location of the colony.
#' @params radius The radius of the colony distance buffer in km
#' @params min.duration The minimum duration in hours that an animal needs to be outside of the colonly buffer before
#' it is counted as a trip

split_trips = function(X, Y, colony.X, colony.Y, time=NULL, radius = 5, colony.radius=NULL, min.duration = 0){
  
    if(is.null(time)) time = Sys.time()
    r = ifelse(is.null(colony.radius), radius, colony.radius)
  
    dist = traipse::track_distance_to(X,Y,colony.X,colony.Y)/1000   
    at.colony = dist < r
    trip = cumsum(at.colony != lag(at.colony,default=TRUE))
    
    tibble(dist,time,trip,at.colony) %>%
    group_by(trip) %>% 
    {if(!is.null(colony.radius)){
        mutate(.,trip = case_when(
             at.colony & row_number() == n()~ trip+1L,
             at.colony & row_number() == 1L ~ trip-1L,
             TRUE ~ trip))
      } else .} %>%
     mutate(duration = suppressWarnings(
                        difftime(max(time[dist>radius]),
                                 min(time[dist>radius]),
                                 unit='hours')))   %>%
     ungroup() %>%
     mutate(returns = !(trip == max(trip,na.rm=T) & !last(dist)<=radius)) %>%
     mutate(trip = ifelse(duration >= min.duration,trip,NA)) %>%
     mutate(returns = ifelse(is.na(trip),NA,returns)) %>%
     mutate(trip = as.numeric(factor(trip))) %>%
     dplyr::select(trip,returns)
  
}

#' @examples This function is designed to be called within a tidy workflow and works with grouping etc. to split trips for multiple animals at the same
#' time. 
#' group_by(data,id) %>% 
#' mutate(split_trips(time = time, X = X, Y = Y, colony.X = 14.212, colony.Y = 9.232, radius = 2, min.duration = 12))

#' @returns A two-column data frame containg the trip ids ('trip') and a logical vector (returns) indicating whether the trip returned to within the colony buffer

