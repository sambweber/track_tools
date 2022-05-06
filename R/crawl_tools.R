
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
#'
#' @param 
#'
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




# --------------------------------------------------------------------------------------------------------------
# as_crwData
# --------------------------------------------------------------------------------------------------------------

#' Constructor function for converting lists of \code{crwFit} and \code{crwPredict} objects for multiple individuals
#' into a \code{momentuHMM::crwData} object for multiple imputation. This is for use with nested tibbles containing
#' list columns with crawl model outputs for multiple individuals and generates the same output as the \code{momentuHMM::crwWrap} function
#'
#' @param ID,crwFit,crwPredict The names of the columns in \code{data} containing the individual track IDs and the lists of crwFit and crwPredict objects
#'
#' @return A formatted crwData object for passing to other \code{momentuHMM} functions

# -------------------------------------------------------------------------------------------------------------

as_crwData = function(data,ID,crwFit,crwPredict) {
 
 if(!all(c(ID,crwFit,crwPredict) %in% names(data))) stop ("Column names provided do not exist in data")
 if(!all(map_lgl(data[[crwPredict]],is,'crwPredict'))) stop ("All crwPredict objects should be of class 'crwPredict'")
 if(!all(map_lgl(data[[crwFit]],is,'crwFit'))) stop ("All crwFit objects should be of class 'crwFit'")
 
 tcol = unique(map_chr(test$crwFit,~.x$Time.name))
 if(length(tcol)>1) stop("Time.name is not the same for all crwFit objects")
 
 data[[crwPredict]] <- map2(data[[ID]],data[[crwPredict]],~mutate(data.frame(.y),ID = as.character(.x)))
     
 crw =
  list(
    crwFits = setNames(data[[crwFit]],data[[ID]]),
    crwPredict = bind_rows(data[[crwPredict]])
  )     
      
 class(crw) <- append('crwData',class(crw))
 attributes(crw$crwPredict)$Time.name <- tcol 
 crw                                  
     
}
      
      
# --------------------------------------------------------------------------------------------------------------
# hmm_inits
# --------------------------------------------------------------------------------------------------------------

#' Uses k-means clustering to generate appropriate starting values for the state-dependent probability distribution parameters of each data stream
#' in hidden Markov models. Currently only works for step length and turning angle, but I plan to extend it to hanlde additional data streams. 
#'
#' @param n.states The number of latent states to initialize
#'
#' @param step,angle The names of columns containing the step lengths and turning angles. By default these are set to the 
#' column names output by \code{momentuHMM::prepData}
#'
#' @return A named list containing the means and standard deviations of step length and turning angle for each state which can be passed straight to
#' the \code{Par0} argument in \code{momentuHMM::fitHMM}
#
#' This function needs enhancing to handle additional data streams (e.g. altitude)

# -------------------------------------------------------------------------------------------------------------

hmm_inits = function(data,nstates,step = 'step',angle = 'angle'){
      
  k         <- na.omit(dplyr::select(data,step,angle)) 
  k1        <- mutate(k,angle=abs(angle))
  clus      <- kmeans(k1,nstates)
  stp_init  <- clus$centers[,1]
  stp_sd    <- aggregate(k[,1]~fitted(clus,method="classes"),FUN='sd')[,2]
  stepPar0  <- c(stp_init,stp_sd)
  anglePar0 <- aggregate(k[,2]~fitted(clus,method="classes"),FUN=CircStats::est.kappa)[,2]
  list(step = stepPar0, angle = anglePar0) 
      
}


