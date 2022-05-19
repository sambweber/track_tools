daynight <- function(object,datetime,twilight = "Civil",moonphase=F){
  
  require(maptools)
  require(lunar)
  require(sf)
  
  if(!datetime %in% names(object)) stop(cat('Column ',datetime,' does not exist'))
  
  spatial = FALSE
  
  if(is(object,'SpatialPoints')) {spdf = st_as_sf(object); spatial = TRUE}
  
  if(is(object,"sf")){
    longlat = st_is_longlat(object)
    if(!longlat){
      crs = st_crs(object)
      object = st_transform(object,4326)
    }
    crds = st_coordinates(object)
  } else {
    stop("Object should be of class 'sf' or 'SpatialPoints'")
  }
    
  dt   = object[[datetime]]
  
  if(!is(dt,"POSIXt")) stop ("datetime must be a vector of class POSIXt")
  
  if(twilight == "Civil"){m = 6} else if (twilight == "Nautical"){m = 12} else stop ("twilight mode should be either 'Civil' or 'Nautical'")
  
  dawn = crepuscule(crds = crds, dt,direction = "dawn",POSIXct.out=T,solarDep = 6)
  rise = sunriset(crds = crds, dt,direction = "sunrise",POSIXct.out=T)
  set  = sunriset(crds = crds, dt, direction = "sunset",POSIXct.out=T)
  dusk = crepuscule(crds = crds, dt,direction = "dusk",POSIXct.out=T,solarDep = 6)
  
  object[["daynight"]] = ifelse(dt >= rise$time & dt < set$time, "Day",ifelse(
                        dt < dawn$time | dt >= dusk$time, "Night","Twilight"))
  
  if(moonphase) object[["moonphase"]] = lunar.phase(dt, name = 8)
  
  if(!longlat) object = st_transform(object,crs)
  
  if(spatial) object = as(object,'Spatial')
  
  return(object)
}
