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

#' Filters tracking dataset to retain only the best location within a series of locations that are separated by a minimum time interval.
#' This approach can be used to thin a dataset to reduce temporal pseudoreplication (e.g. best hourly location) or as a preparatory step 
#' for State Space Models.
#' Having lots of locations that are very closely spaced in time (< 2 seconds) can lead to problems with SSM model fitting that
#' can't be solved by simply adding a fixed amount to duplicated timestamps. 

#' The function retains the best location within a rolling window of duration tmin, based on quality columns specified in
#' filter_cols. The location with the MAXIMUM value in filter_cols is retained in each hourly window. In the event of a tie, 
#' the first is retained. Multiple filter columns can be specified by wrapping in c() in which case they are applied in order. 
#' The location(s) with the maximum value in the first filter_col is selected first, and subsequent filter_cols are used to break ties.
#' This will work on numeric and factorial filter columns, but factors should be ordered with the best/highest value as the final level. 

#' @param data
#' @param dt The name of the column containing the timestamp (must be of class \code{POSIX*}) 

best_location = function(df,timestamp,tmin='1 hour',filter_cols){

df <- arrange(df,!!sym(timestamp))
selected <- df[0,]

start_time <- min(df[[timestamp]])
end_time <- max(df[[timestamp]])

while (start_time < end_time) {
  
  window_end <- start_time + duration(tmin)
  window_data <- filter(df,!!sym(timestamp) >= start_time, !!sym(timestamp) < window_end)
  
  if (nrow(window_data) > 0) {
  
    best_point <- slice_max(window_data,tibble(!!!syms(filter_cols)),with_ties = FALSE)
    selected <- bind_rows(selected, best_point)
    start_time <- best_point[[timestamp]] + duration(tmin)
    
  } else {
    start_time <- window_end
  }
}

return(selected)

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
  
  
