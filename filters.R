
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#land_filter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Removes or flags points that occur on land

#Takes 2 arguments:

#Points is a SpatialPoints* object
#landmask is a list of one or more SpatialPolygons* objects 

#If remove = F returns a logical vector indicating whether each point occurs on land
#If remove = T returns a SpatialPointsDataFrame minus the filtered points


land_filter<-function(points,landmask,remove=F){
  
require(raster)
require(rgeos)
require(maptools)
  
stopifnot(is(points,"SpatialPoints")|is(points,"SpatialPointsDataFrame"))
stopifnot(is(landmask,"list"))

filterpass <- rep(T,length(points))

for (i in 1:length(landmask)){
  
      l   = landmask[[i]]

      stopifnot(is(l,"SpatialPolygons")|is(l,"SpatialPolygonsDataFrame"))
      
            if(nrow(l)>1){
  
            l = unionSpatialPolygons(l,IDs=rep(1,nrow(l)))
  
            }
      
      
            if(projection(points)!=projection(l)){
        
            l = spTransform(l,CRS(projection(points)))
        
            }

      
      result <- !gIntersects(points,l,byid=T)
      
      filterpass <- ifelse(filterpass!=F,result,filterpass)

}

if(remove==T){output = points[which(filterpass),]} else {output = filterpass}

return(output)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#simple_speed_filter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Simple speed filter which checks if the speed between consecutive track points exceeds max.speed and returns a new column called
# "filterpass" indicating whether each location passed the filter.The process is applied iteratively until the filter passes and then 
# moves on to the next segment i.e. if speed between t and t+1 exceeds max.speed, t+1 is marked as fail, and the speed between t and t+2
# is checked, and so on until a pass is achieved. The filter then restarts from the point of the pass.


simple_speed_filter<-function(points,dt,id,max.speed){
  
  stopifnot(is(points,"SpatialPoints")|is(points,"SpatialPointsDataFrame"))
  
  ll = ifelse(is.projected(points),FALSE,TRUE)
  
  if(is(points[[id]],"factor")){points[[id]] = factor(points[[id]])}
  
  data = split(points,points[[id]])
  
  data1 = lapply(data,FUN = function(x){
  
  x$filterpass = TRUE
  
  t = 1
  
   for (i in 2:nrow(x)){
     
     speed = spDists(x[t,],x[i,],ll)/as.numeric(difftime(x[i,][[dt]],x[t,][[dt]],units="hours"))
     
     if(speed >= max.speed){x$filterpass[i] = FALSE} else {t = i}
     
   }
  
  x
    
    
  })
  
  return(do.call(rbind,data1))
  
}
  


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# speed_filter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Removes locations based on a maximum speed between fixes

# speed_filter uses the speedfilter function from package trip for all tracks where there are enough points. 
# speedfiler is a forwards-backwards speed-averaging algorithm that calculates the velocity required to move between consecutive 
# fixes and removes fixes where velocity exceeds a given threshold. This algorithm requires a minimum number of points to run
# and fails where there are too few. To avoid this happening speed_filter uses a simple, forward only, custom function
# simple_speed_filter to remove unlikely locations for short tracks. 

# A distance threshold can be specified for auto-accepting locations that fall within a given distance of the previous
# and/or subsequent fixes. This allows for high estimated speeds that are due to points closely spaced in time but
# with overlapping errors.

# speed.filter also first removes "Z" classes -  which are the points for which the location process
# failed - before running the filtering algorithms.


speed_filter<-function(points,dt,id,lc,max.speed,min.dist=NULL){
  
  require(trip)
  
  stopifnot(is(points,"SpatialPoints")|is(points,"SpatialPointsDataFrame")|is(points,"sf"))
  
  is.sf = is(points,'sf')
  
  if(is.sf) {points = as(points,'Spatial')}
  
  ll = !is.projected(points)
  
  points = points[!points[[lc]] %in% c("-9","z","Z"),] # Remove z locations before filter
  
  if(is(points[[id]],"factor")){points[[id]] = factor(points[[id]])} #remove unwanted factor levels for split
  
  data = split(points,points[[id]])
  
  output = lapply(data, FUN = function(x){
   
  if(nrow(x)>1) {
    
    if(nrow(x) < 4) { 
      
        x = simple.speed.filter(x,dt=dt,id=id,max.speed=max.speed)
      
      } else {
      
        x$filterpass = trip::speedfilter(trip(x,c(dt,id)),max.speed=max.speed)
      
      }
    
    if(!is.null(min.dist)){
    
    if(!ll) min.dist = min.dist*1000
    accept = spDists(x,longlat=ll,segments=T) < min.dist
    x$filterpass[1] = any(x$filterpass[1],accept[1])
    x$filterpass[nrow(x)] = any(x$filterpass[nrow(x)],accept[length(accept)])
    if(nrow(x)>2){
    accept = accept[-1] | accept[-length(accept)]
    x$filterpass[-c(1,nrow(x))] =  mapply(any,x$filterpass[-c(1,nrow(x))],accept)
    }
    
    }
    
    x
    
  } else {
    
    x$filterpass = TRUE
    x
  }
  
})
      
  output = do.call(rbind,output)
  if(is.sf) output = st_as_sf(output)
  return(output)
      
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remove.early
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Removes points at the start of a track to allow for behavioural equilibration

remove.early = function(points,deploy = NULL, dt,id,adjust.period = 24){
  
  ids = unique(id)
  
  if (is.null(deploy)) {
  deploy = aggregate(points[[dt]] ~ points[[id]], FUN= function(x) {min(x)})
  names(deploy) = c(id,dt)
  }
  
  deploy.times = deploy[match(points[[id]],deploy[[id]]),dt]
  
  t = difftime(points[[dt]],deploy.times,unit="hours")
  
  t == 0 | t > adjust.period
  
 }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dist.filter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# if max.dist is not specified, outliers are classified as fixes that are more than 1.5
# interquartile ranges above the third quartile in interfix distances from all other
# points in the tracking dataset

dist.filter = function(data,id,longlat=TRUE,max.dist=NULL){
  
  if(is.factor(data[[id]])) data[[id]] = factor(data[[id]]) # drop redundant levels
  
  if(longlat == F & !is.null(max.dist)) max.dist = max.dist*1000
  
  do.call('c',lapply(split(data,data[[id]]), function(x) {
    
    d = spDists(x,longlat=longlat)
    if(is.null(max.dist)) {t=d[lower.tri(d,diag=F)]; max.dist = quantile(t,0.75) + IQR(t)*1.5}
    diag(d) = NA
    d = apply(d,1L,min,na.rm=T) < max.dist
    d[1L] = TRUE
    d

  })
  )
}   


# ---------------------------------------------------------------------
# preptrack_sf
# ---------------------------------------------------------------------

# This function is a tidyverse sf implementation of the trip::forceCompliance function.
# It removes duplicate rows, orders points by date-times within ID and adjusts duplicated timestamps by adding 1 second to them

preptrack_sf = function (x, id = 'id', dt = 'datetime', minlocs = 1){

if(is(x,'Spatial')) x = st_as_sf(x)   
  
x %>%
  group_by_(id) %>%
  dplyr::filter(n() >= minlocs) %>%         # remove tracks that are too short
  arrange(!!as.name(dt),.by_group=T) %>%    # order by id and datetime
  dplyr::distinct(.keep_all=TRUE) %>%       # remove duplicate rows
  group_by_(id,dt) %>%
  mutate(replicate = seq(n())-1) %>%        #create 1 second sequence along duplicated timestamps
  ungroup() %>%
  mutate(!!dt := !!as.name(dt) + replicate) %>%  #add sequence to timestamps
  dplyr::select(-replicate)

}



