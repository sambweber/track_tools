
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


# ----------------------------------------------------------------------------------------------------
# crw_effSamp
# ----------------------------------------------------------------------------------------------------

#' A convenience function to extract the effective sample size (i.e. approximate number of independent samples)
#' from a posterior simulation object created by \code{crawl::crwSimulator}

#' @return The original dataframe with new columns ln.sd.x and ln.sd.y which are log(err/sqrt(2)), where err is
#' Argos errors in metres, along with the covariances between them.

# -------------------------------------------------------------------------------------------------------

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




