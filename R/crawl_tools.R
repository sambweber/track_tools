
# ----------------------------------------------------------------------------------------------------
# add_anisotopic_errors
# ----------------------------------------------------------------------------------------------------

#' A convenience wrapper around \code{crawl::argosDiag2Cov} that adds anisotropic errors back to a dataframe
#' of tracking data in a form that \code{crawl::crwMLE} can use.

#' @return The original dataframe with new columns ln.sd.x and ln.sd.y which are log(err/sqrt(2)), where err is
#' Argos errors in metres, along with the covariances between them.

# -------------------------------------------------------------------------------------------------------

add_anisotopic_errors <- function(obj, Major = "Error.Semi.Major.Axis", Minor = "Error.Semi.Minor.Axis",
Orientation = "Error.Ellipse.Orientation"){

cbind(obj, 
      crawl::argosDiag2Cov(
        obj[[Major]], 
        obj[[Minor]], 
        obj[[Orientation]])
      )

}


# --------------------------------------------------------------------------------------------------------------
# crw_effSamp
# --------------------------------------------------------------------------------------------------------------

#' Convenience functions to extract the effective sample size (i.e. approximate number of independent samples)
#' from a posterior simulation object created by \code{crawl::crwSimulator}

#' @return 

# -------------------------------------------------------------------------------------------------------------

crw_effSamp <- function(object){
  if(!is(object,'crwSimulator')) stop ("object must be of class crwSimulator")
  attr(object$thetaSampList[[1]],'effSamp')
}

crw_check_ISW <- function(object,plot=T){
  if(!is(object,'crwSimulator')) stop ("object must be of class crwSimulator")
  
  w = simObj$thetaSampList[[1]]
  d = density(w[,1]*nrow(w),na.rm=T,n=1000)
  dmax = d$x[which.max(d$y)]
      
  if(plot){
    plot(d, xlab='Weight',main='Importance Sampling Weights', sub='(Weights centered on 1 is desirable)')
    polygon(d,col = 'grey80')  
  } else return(dmax)
      
}


# --------------------------------------------------------------------------------------------------------------
# get_sim_track
# --------------------------------------------------------------------------------------------------------------

#' Draws from a posterior simulation object created using \code{crawl::crwSimulator} and formats 
#' simulated tracks into a \code{sf} object.

#' @return 

# -------------------------------------------------------------------------------------------------------------
get_sim_track = function(obj,iter,fullPost=F) {

lapply(1:iter,function(i){
    j = 0
    while(j <= 10){
      sim = try(crwPostIS(obj,fullPost = T))
      if(class(sim)[1] == "try-error") j = j+1 else break
    }
  
    sim$alpha.sim[sim$locType=='p',c(1,3)] %>% as.data.frame() %>%
    mutate(iter = i,datetime = sim$TimeNum[sim$locType=='p']*3600,
           datetime = as.POSIXct(datetime,origin = '1970-01-01',tz='GMT')) 
  }) %>% 
    
  do.call(rbind,.) %>%
  st_as_sf(coords=c('mu.x','mu.y'),crs = attr(obj,'proj4'))
      
 }
