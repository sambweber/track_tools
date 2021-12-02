
# ----------------------------------------------------------------------------------------------------
# Interpolate track
# ----------------------------------------------------------------------------------------------------

#' This function interpolates a track into constant time steps using a range of linear or spline interpolation (smooth)
#' functions.

#' @param spdf a \code{SpatialPoints*} or \code{sf} object containing the point locations along the track
#' 
#' @param datefield a character vector containing the name of the column containing the date/time in POSIXct format
#' 
#' @param id a character vector containing the name of the column containing a unique identifier for each trip or track segment.
#' Interpolation is performed within each id.
#' 
#' @param t the constant time step (in seconds) to interpolate to
#' 
#' @param fun the function used to perform the interpolation. Can be a choice of:
#' 
#' 'linear' (via \code{approx} function),
#' 'pchip' for piecewise cubic hermite interpolation vie the \code{signal::pchip} function
#' 'spline' for cubic spline interpolation vie the \code{stats::spline} function
#' 'smooth.spline' for cubic smoothing spline interpolation vie the \code{stats::smooth.spline} function
#' 'loess' for local polynomial interpolation vie the \code{stats::loess} function
#' 
#' @param data.cols A character vector containing the names of any additional data columns from \code{spdf} object that 
#' should be included in the interpolated output. For factor and character variables, the first value in the field
#' is copied to all interpolated points. For numeric columns, values are interpolated using the function provided in \code{fun}.
#' 
#' @param geom The output geometry for the interpolated track. Can be either 'points' or 'line'.
#' 
#' @param ... Optional additional arguments passed to \code{fun} to control smoothing.
#' 
#' @return Returns a track interpolated to constant timestep \code{t} for each \code{id} in \code{spdf}


# -------------------------------------------------------------------------------------------------------

interpolate_track <- function(spdf,datefield,id,t,fun=c('linear','pchip','spline','smooth.spline','loess'),
                              data.cols=NULL,geom=c('points','line'),...){
  
  require(raster)
  require(signal)
  require(sf)
  require(dplyr)
  
  if(!(is(spdf,'sf') | is(spdf,'SpatialPoints')) stop ("spdf must be an object of class 'sf' or 'SpatialPoints'")
  sf = is(spdf,'sf'); if(sf) spdf = as(spdf,'Spatial')
  fun = match.arg(fun)
  geom = match.arg(geom)
  
  if(!is(spdf[[datefield]],"POSIXt")) stop ("datefield must be an object of class POSIX")
  hasdata<-ifelse(is.null(data.cols),F,T)
  if(hasdata & !is(data.cols,"character")) stop ("data.cols should be a character vector of field names")
  
  spdf[[id]] <- factor(spdf[[id]])
  tracklist<-split(spdf,spdf[[id]])
  singles = sapply(tracklist,nrow)==1
  if(any(singles)) message("Some ids have only a single point associated with them and will be removed")
  tracklist = tracklist[!singles]
  
  output = lapply(tracklist,FUN = function(track) {
    
    x     = coordinates(track)[,1]
    y     = coordinates(track)[,2]
    s     = as.numeric(track[[datefield]])
    tnew  = seq(min(s),max(s),t)
    
    if(fun=='linear'){X = approx(s,x,xout=tnew,method='linear',...)$y;Y = approx(s,y,xout=tnew,method='linear',...)$y}
    if(fun=='pchip'){X = pchip(s,x,tnew);Y = pchip(s,y,tnew)}
    if(fun=='spline'){X = spline(s,x,xout=tnew,...)$y ; Y = spline(s,y,xout=tnew,...)$y} 
    if(fun=='smooth'){X = predict(smooth.spline(s,x,...),tnew)$y; Y = predict(smooth.spline(s,y,...),tnew)$y}
    if(fun=='loess'){X = predict(loess(x~s,...),tnew); Y = predict(loess(y~s,...),tnew)}
    
    result = data.frame("x" = X, "y" = Y)
    result[[id]] = unique(track[[id]]); result[[datefield]] = as.POSIXct(tnew,origin="1970-1-1")
    
    coordinates(result) = c("x","y")
    
    if(hasdata){
      
      for(field in data.cols){
      
      if(!is(spdf[[field]],"numeric")) result[[field]] = track[[field]][1] else {
      if(fun =='pchip'){result[[field]]<-pchip(s,track[[field]],tnew)}
      if(fun =='spline'){result[[field]]<-pchip(s,track[[field]],tnew)}
      if(fun=='linear'){result[[field]] = approx(s,track[[field]],xout=tnew,method='linear',...)$y}
      }
      }
     }  
    
    result
    
  })
  
  output = do.call("rbind",output)  
  if(!is.na(projection(spdf))){ projection(output)=projection(spdf) }
  
  if(geom == 'line'){
  output = st_as_sf(output) %>% group_by(!!as.name(id)) %>% summarise(do_union=F) %>% st_cast('LINESTRING') %>% dplyr::select(-do_union)
  }
  
  if(sf) return(st_as_sf(output)) else return(output)
  
}
