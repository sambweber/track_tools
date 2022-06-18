# ---------------------------------------------------------------------------------------
# kernelize
# ---------------------------------------------------------------------------------------

kernelize = function(data,id,resolution,h=NULL,crs=NULL){

if(is(x,'sf')) {

  crs = st_crs(data)
  data = dplyr::select(-any_of(c('X','Y'))) %>% 
         cbind(st_coordinates(.)) %>%
         st_set_geometry(NULL)
  
}

if(!all(c('X','Y') %in% names(data))) stop ("data should contain columns 'X' and 'Y' containing coordinates)

# compose grid
lims = c(range(data$X),range(data$Y))
r = extent(lims) %>% raster(resolution = resolution) %>% extend(1)
lims = as.numeric(matrix(extent(r)))
n = dim(r)[2:1]

fit.k = function(d) {
   MASS::kde2d(x = d$X, y = d$Y, n=n,lims = lims) %>% 
    raster(crs=crs)
 }

if(!missing(id)){

    nest(data,crds = -!!id) %>%
    mutate(kernel = map(crds,fit.k)) %>%
    dplyr::select(-crds)
    
} else { fit.k(data) }

}

