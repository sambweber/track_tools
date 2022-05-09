# ----------------------------------------------------------------------------------------------------
# adjust_duplicate_times
# ----------------------------------------------------------------------------------------------------

#' This is just the \code{adjust.duplicateTimes} function from package \code{trip}

#' @param dt A datetime vector of class \code{POSIX*} for correction
#' 
#' @return Returns vector with 1 second added to duplicate timestamps 

# -------------------------------------------------------------------------------------------------------

adjust_duplicate_times <- function (dt,eps=10) {
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

fix_track = function (x, dt = 'datetime', drop.duplicate.times = F, eps = 10){
  
  dt = as.name(dt)
  if(is(x,'Spatial')) x = st_as_sf(x)   
  
  x %>%
    mutate(!!dt := round(!!dt,'secs')) %>%    # deal with fractional seconds from some tags (e.g. FastGPS)
    dplyr::distinct(.keep_all=TRUE) %>%       # remove duplicate rows
    {if(drop.duplicate.times) dplyr::filter(.,!duplicated(!!dt)) else .} %>%
    mutate(!!dt := adjust_duplicate_times(!!dt, eps=eps)) %>% # adjust duplicated time stamps
    arrange(!!dt,.by_group=T)              # order by datetime
  
}

