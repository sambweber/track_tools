land_mask = function(object,res){
    
    require('fasterize')
    if(is(object,'Spatial')) object = st_as_sf(object)
    st_union(object) %>% st_sf %>%
    fasterize(fasterize::raster(object,res=c(res,res)),background=0)
    
}

# --------------------------------------------------------------------------------------

on_land = function(pts,land,res){
     
    if(!is(pts,'sf') | !is(land,'sf')) stop("'pts' and 'land' should be of class 'sf'")
    land_mask(land,res) %>%
    terra::rast() %>%
    terra::extract(st_coordinates(pts)) %>%
    unlist(use.names=F)
    
}

# ---------------------------------------------------------------------------------------

land_resample = function(pts,land,loc.id,max.iter=1000){
    
    pts = mutate(pts,on.land = on_land(pts,land,res=0.001))
    
    if(any(pts$on.land)){
        
        
    pts %>% cbind(st_coordinates(.)) %>%
    group_by(!!loc.id) %>% 
    summarize(across(c('X','Y'),c(mu=mean,sd=sd)),vcov = list(cov(cbind(X,Y))))
    
    }
