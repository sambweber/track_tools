require(raster)
require(tidyverse)
library(sf)

# ---------------------------------------------------------------------------------------
# kernelize
# ---------------------------------------------------------------------------------------

# data is a dataframe containing columns labelled 'X' and 'Y' containing the coordinates of the points.
# data may also be an object of class 'sf' in which case coordinates are extracted from the geometry

# id is an optional character vector naming a column of individual identifiers. A seperate kernel density will be
# produced for each individual, but calculated on a common grid. This can assist in later merging of the kernels.

# resolution is the resolution of the grid on which to compute the kernels, given in the units of the coordinate system.
# If coordinates are in degrees, then resolution must be given in degrees.

# h is the bandwidth used for kernelling. Defaults to the reference bandwidth. See ?MASS::kde2d

# crs is an optional coordinate reference system which is applied to the final raster. If data is 'sf' then the crs is taken
# directly from the input.

# fixed.grid is a logical indicating whether all US should be estimated on the same grid. This is useful for merging UDs
# later (i.e. where multiple individuals overlap), but takes much longer when individuals are spread over a large and disparate
# area. It should only be set to TRUE when you may want to calculate population level UDs or overlaps later.

kernelize <- function(data,id,resolution,h=NULL,crs=NULL,extend=1,fixed.grid = FALSE){
    
    if(is(data,'sf')) {
        
        crs = st_crs(data)
        data = dplyr::select(data,-any_of(c('X','Y'))) %>% 
            cbind(st_coordinates(.)) %>%
            st_set_geometry(NULL)
        
    }
    
    if(!all(c('X','Y') %in% names(data))) stop ("data should contain columns 'X' and 'Y' containing coordinates")
    
    # compose grid
    mk.grd = function(data){
          lims = c(range(data$X),range(data$Y))
          r = terra::ext(lims) %>% terra::rast(resolution = resolution) 
          r = terra::extend(r,c(round((nrow(r)*extend-nrow(r))/2), round((ncol(r)*extend-ncol(r))/2)))
          lims = as.numeric(matrix(terra::ext(r)))
          n = dim(r)[2:1]
          list(lims=lims,n=n)
    }

    fit.k = function(d,grd) {
        k = MASS::kde2d(x = d$X, y = d$Y, n=grd$n,lims = grd$lims) %>% raster::raster()
        (k/cellStats(k,'sum'))
    }
    
    if(!missing(id)){
        
        nest(data,crds = -!!id) %>%
          {if(fixed.grid) {
            mutate(.,kernel = purrr::map(crds,~fit.k(.x, grd = mk.grd(data))))
              } else {
            mutate(.,kernel = purrr::map(crds,~fit.k(.x, grd = mk.grd(.x))))
              }
           } %>%
        dplyr::select(-crds)
        
    } else {fit.k(data, grd = mk.grd(data)) }
    
}

# -------------------------------------------------------------------------------------------------------
# volume_contour
# -------------------------------------------------------------------------------------------------------

# input must be a Raster*

# res.out is a numeric indicating how many times to increase the resolution of the input raster by for producing smoother contours and UDs. 
# If res.out = 10 and the input resolution is 100 x 100 m, then the raster will be resampled onto a 10 x 10 m grid for contouring.

# levels are the precentage volume contours to extract

# output is either a classified raster or polygons

# landmask is an optional spatial polygon layer of class 'sf' which is used to clip the resulting volume contours so they don't overlap the coastline. Adding
# landmask can significantly increase computation time if the resolution is high, but it is necessary to get accurate area estimates for the volumne contours.
# The landmask is first used to mask the input raster and remove any pixels on land. If output type is 'polygons' the landmask is also used to clip the resulting
# polygons.

volume_contour <- function(input, res.out = 10, levels = c(95,75,25,50),output = c('raster','polygons'),landmask = NULL){
  
  require(raster); require(rgdal); require(sf); require(rgeos)
  
  if(!is.null(landmask)){
        if(!is(landmask,'sf')) stop('landmask should be of class sf')
        landmask = st_as_sfc(landmask) %>% st_union()  #Make sure it is just a single geom with no data
    }  
    
  output = match.arg(output)
  
  resamp = raster::disaggregate(x=input,fact = res.out,method="bilinear")
  if(!is.null(landmask)) Z = mask(resamp,as(landmask,'Spatial'),inverse=T,updatevalue = 0)
  Z = raster::values(resamp)
  resamp[] = Z/sum(Z,na.rm=T) # normalise so that all values sum to 1
  
  vclevels = sort(levels,decreasing = T) #sort the volume contour levels
  rastvals = sort(raster::values(resamp)) #rank the raster probability values
  breaks   = sapply(1-vclevels/100,FUN=function(x){rastvals[max(which(cumsum(rastvals)<=x))]})
  breaks   = c(breaks,1)
  
  if(output == 'polygons'){
    #create contour lines
    vcs = rasterToContour(resamp,levels = breaks[1:length(levels)],maxpixels = ncell(resamp)) %>%
          st_as_sf()
    vcs$level = vclevels
    
    #convert to polygons 
    vcpolys =  st_cast(vcs,'POLYGON') %>% st_make_valid() 
    if(!is.null(landmask)) vcpolys = st_difference(vcpolys,landmask)
    
    #Calculate areas using equal area projection
    centroid = st_coordinates(st_centroid(st_as_sfc(st_bbox(vcpolys))))
    lamberts = st_crs(paste0("+proj=laea +lat_0=",centroid[2]," +lon_0=",centroid[1]))
    vcpolys$area_km2 = round(st_area(st_transform(vcpolys,lamberts))/1E6,1)
    
    return(vcpolys)
    
  } else {
    
    resamp = reclassify(resamp,matrix(c(min(rastvals),breaks[-length(breaks)],breaks,NA,vclevels),ncol=3),include.lowest=T)
    return(resamp)
    
  }
  
}

# ---------------------------------------------------------------------------------------------------------------------------
# ref_bandwidth
# ---------------------------------------------------------------------------------------------------------------------------

ref_bandwidth = function(data,id){
    
     if(is(data,'sf')) {
        
        crs = st_crs(data)
        data = dplyr::select(-any_of(c('X','Y'))) %>% 
            cbind(st_coordinates(.)) %>%
            st_set_geometry(NULL)
        
    }
    
    if(!all(c('X','Y') %in% names(data))) stop ("data should contain columns 'X' and 'Y' containing coordinates")
    
    bw = function(d) c(MASS::bandwidth.nrd(d$X),MASS::bandwidth.nrd(d$Y))
    
    if(!missing(id)){
        split(data,data[[id]]) %>% purrr::map(bw) %>% bind_rows %>% colMeans
    } else {bw(data)}
        
   }

# ---------------------------------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------------------------------

scaleARS = function(data,scales = 1:24*3600,datetime,id,plot=T){
  
  require(adehabitatLT)
  
  if(!is(data[[datetime]],'POSIXt')) stop(paste(datetime,'should be of class POSIXt'))
  
  if(is(data,'sf')) {
    crds = st_transform(data,3395) %>% st_coordinates
  } else if (is(data,'data.frame') & all(c('lon','lat') %in% names(data))) {
    crds = st_as_sf(data,coords=c('lon','lat'),crs=4326) %>% st_transform(3395) %>% st_coordinates()   
  } else stop("data should be in format 'sf' or a dataframe with columns 'lon' and 'lat'")
    
traj <- as.ltraj(data.frame(crds[,1], crds[,2]), date=data[[datetime]], id=data[[id]], typeII = TRUE)

fp <- fpt(traj, radii = scales, units = "seconds")
var_fp <- varlogfpt(fp,graph=F) #find the variance in log(fpt) for each radius
peaks = scales[apply(var_fp,1,which.max)]
ars = round(mean(peaks)) 

if(plot){
  pl =
  setNames(var_fp,scales) %>% mutate(id = rownames(.)) %>%
  gather(radius,varfpt,-id) %>% 
  mutate(radius = as.numeric(radius)) %>% 
  ggplot(aes(x = radius,y=varfpt,colour=id)) + 
  geom_line()
  print(pl)
}

return(ars)

}
