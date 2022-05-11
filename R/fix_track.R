# ----------------------------------------------------------------------------------------------------
# adjust_duplicate_times
# ----------------------------------------------------------------------------------------------------

#' This is just the \code{adjust.duplicateTimes} function from package \code{trip}

#' @param dt A datetime vector of class \code{POSIX*} for correction
#' 
#' @return Returns vector with 1 second added to duplicate timestamps 

# -------------------------------------------------------------------------------------------------------

adjust_duplicate_times <- function (dt,eps=1) {
  dups <- duplicated(dt)
  if (any(dups)) {
    dt[dups] <- dt[dups] + eps
    dt <- Recall(dt)
  }
  return(dt)
}

# ----------------------------------------------------------------------------------------------------
# fix_track
# ----------------------------------------------------------------------------------------------------

#' Fixes common problems in tracking data including fractional seconds, duplicated rows and dulicated timestamps.
#' A tidy/sf implementation of the \code{trip::forceCompliance} function

#' @param x An object of class \code{Spatial} or \code{sf} containing track points for a single individual
#' @param dt The name of the column containing the timestamp (must be of class \code{POSIX*}) 

#' @return Returns the original spatial object with common problems fixed

# -------------------------------------------------------------------------------------------------------

fix_track = function (x, dt = 'datetime', eps = 1){
  
  dt = as.name(dt)
  if(is(x,'Spatial')) x = st_as_sf(x)   
  
  x %>%
    mutate(!!dt := round(!!dt,'secs')) %>%    # deal with fractional seconds from some tags (e.g. FastGPS)
    dplyr::distinct(.keep_all=TRUE) %>%       # remove duplicate rows
    mutate(!!dt := adjust_duplicate_times(!!dt, eps=eps)) %>% # adjust duplicated time stamps
    arrange(!!dt,.by_group=T)              # order by datetime
  
}
  
  
# ----------------------------------------------------------------------------------------------------
# best_location
# ----------------------------------------------------------------------------------------------------

#' Filters tracking dataset to retain only the best location within a burst of locations that are closely spaced in time.
#' Having lots of locations that are very closely spaced in time (< 2 seconds) can lead to problems with model fitting that
#' can't be solved by simply adding a fixed amount to duplicated timestamps. In the event of a tie, the first observation
#' in each burst is retained.

#' @param data
#' @param dt The name of the column containing the timestamp (must be of class \code{POSIX*}) 

#' @return 

# -------------------------------------------------------------------------------------------------------

best_location = function(data,dt,tmin=1,filter_cols){
     
    if(!all(filter_cols %in% names(data))) stop ('named filter_cols do not all appear in data')
    min_na = function(x) {if(all(is.na(x))) is.na(x) else x == min(x,na.rm=T)}
  
    mutate(data,across(where(is.factor),as.numeric)) %>%
    mutate(ok = traipse::track_time(!!as.name(dt)) >= tmin) %>%
    mutate(ok = replace_na(ok,TRUE))   %>%
    mutate(ok = cumsum(ok))  %>%
    group_by(ok) %>% 
    dplyr::filter(across(all_of(filter_cols),min_na)) %>%
    slice(1) %>% ungroup()
  
}
 
# ----------------------------------------------------------------------------------------------------
# spread_times
# ----------------------------------------------------------------------------------------------------

#' Function for adjusting closely spaced times so that they are separated by a user-specified minimum
#' interval. This is an extension of adjust_duplicate_times which ensures that .

#' @param data
#' @param dt The name of the column containing the timestamp (must be of class \code{POSIX*}) 

#' @return 

# -------------------------------------------------------------------------------------------------------  
  
spread_times = function(data,dt,tmin = 1){

data = mutate(data,t = traipse::track_time(!!as.name(dt)))

i = 2

while(i < nrow(data)){
  
  if(data$t[i] >= tmin) {i=i+1;next} else {
    
    w = i:nrow(data)
    x = cumsum(data$t[w])
    s = seq_along(x)*tmin
    data[[dt]][w] <- data[[dt]][w] + pmax((s-x),0)
    i = i + sum(x<s)
    
  }
}

data$t <- NULL  
return(data)  

}
  
  
