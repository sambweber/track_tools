land_mask = function(land,res){
    
    require('fasterize')
    land <- try(st_as_sf(land))
    if(is(object,'try-error')) stop("'land' should be, or be coercible to, class 'sf'")
 
    st_union(land) %>% st_sf %>%
    fasterize(fasterize::raster(land,res=c(res,res)),background=0)
    
}

# --------------------------------------------------------------------------------------

on_land = function(pts,land,res){
    
    pts <- try(st_as_sf(pts))
    if(is(pts,'try-error')) stop("'pts' should be, or be coercible to, class 'sf'")
    
    land_mask(land,res) %>%
    terra::rast() %>%
    terra::extract(st_coordinates(pts)) %>%
    unlist(use.names=F)
    
}
    
# ---------------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------------
      
track_land_resample = function(pts,land,timestamp,max.iter){

pts <- try(st_as_sf(pts))
land <- try(st_as_sf(land))
if(is(pts,'try-error') | is(land,'try-error')) stop("'pts' and 'land' should be, or be coercible to, class 'sf'")

crs = st_crs(pts)
land = st_union(land)
    
pts %>% cbind(st_coordinates(.)) %>% 
st_set_geometry(NULL) %>%
split(.[[timestamp]]) %>% 
map(st_as_sf,crs=crs,coords=c('X','Y')) %>%
map(land_resample,land) %>%
bind_rows()

}  
    
# ---------------------------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------------------------
    
land_resample = function(pts,land,max.iter=100){
  
  crds = st_coordinates(pts)
  mu = colMeans(crds)
  vcov = cov(crds)
  
  on_land = function(p) {st_intersects(p,land,sparse=F)[,1]}
  on.land = on_land(pts)
  
  iter = 1
  
  while(any(on.land) & iter < max.iter){
    
    new = MASS::mvrnorm(sum(on.land),mu=mu,Sigma=vcov) %>%
          matrix(ncol=2) %>%
          st_multipoint %>% st_sfc %>% st_cast('POINT')
    st_geometry(pts[on.land,]) <- new 
    on.land = on_land(pts)
    iter = iter+1
    
  }
  
  pts[!on.land,]
  
}
    
