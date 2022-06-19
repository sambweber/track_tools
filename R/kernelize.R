# ---------------------------------------------------------------------------------------
# kernelize
# ---------------------------------------------------------------------------------------

kernelize = function(data,id,resolution,h=NULL,crs=NULL){

if(is(data,'sf')) {

  crs = st_crs(data)
  data = dplyr::select(-any_of(c('X','Y'))) %>% 
         cbind(st_coordinates(.)) %>%
         st_set_geometry(NULL)
  
}

if(!all(c('X','Y') %in% names(data))) stop ("data should contain columns 'X' and 'Y' containing coordinates")

# compose grid
lims = c(range(data$X),range(data$Y))
r = extent(lims) %>% raster(resolution = resolution) %>% extend(1)
lims = as.numeric(matrix(extent(r)))
n = dim(r)[2:1]

fit.k = function(d) {
   k = MASS::kde2d(x = d$X, y = d$Y, n=n,lims = lims) %>% raster()
   (k/cellStats(k,'sum'))
 }

if(!missing(id)){

    nest(data,crds = -!!id) %>%
    mutate(kernel = map(crds,fit.k)) %>%
    dplyr::select(-crds)
    
} else { fit.k(data) }
 
  
}


# -------------------------------------------------------------------------------------------------------
# volume_contour
# -------------------------------------------------------------------------------------------------------

volume_contour <- function(input, res.out = 10, levels = c(95,75,25,50),output = c('raster','polygons')){
  
  require(raster); require(rgdal); require(sf); require(rgeos)
  
  output = match.arg(output)
  
  resamp = raster::disaggregate(x=input,fact = res.out,method="bilinear")
  Z = raster::values(resamp)
  resamp[] = Z/sum(Z,na.rm=T) # normalise 
  
  vclevels = sort(levels,decreasing = T) #sort the volume contour levels
  rastvals = sort(raster::values(resamp)) #rank the raster probability values
  breaks   = sapply(1-vclevels/100,FUN=function(x){rastvals[max(which(cumsum(rastvals)<=x))]})
  breaks   = c(breaks,1)
  
  if(output == 'polygons'){
    #create contour lines
    vcs = rasterToContour(resamp,levels = breaks[1:length(levels)],maxpixels = ncell(resamp))
    vcs$level = levels
    
    #convert to polygons 
    vcpolys = try(st_as_sf(vcs) %>% st_cast('POLYGON')) 
    if(class(vcpolys) == "try-error") vcpolys = st_as_sf(vcs) %>% st_polygonize() 
    
    vcpolys = st_transform(vcpolys,54034) %>% mutate(area_km2 = round(st_area(.)/1E6,1))
    
    return(vcpolys)
    
  } else {
    
    resamp = reclassify(resamp,matrix(c(min(rastvals),breaks[-length(breaks)],breaks,NA,vclevels),ncol=3),include.lowest=T)
    return(resamp)
    
  }
  
}
