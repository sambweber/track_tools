land_mask = function(object,res){
    
    require('fasterize')
    if(is(object,'Spatial')) object = st_as_sf(object)
    st_union(object) %>% st_sf %>%
    fasterize(fasterize::raster(object,res=c(res,res)))
    
}
