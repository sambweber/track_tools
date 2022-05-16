# ----------------------------------------------------------------------------------------------------
# make_argos_errors
# ----------------------------------------------------------------------------------------------------

#' A function to fill in missing error information in Argos/FastGPS tracking data prior to state space modelling. Missing
#' data is filled in three steps:
#'
#' 1. Location classes with no error information at all are first replaced by named values passed in LC.values.
#  2. Next, observations with no error estimates, but for which isotropic error radii are available for other observations of the same
#  location class are assigned the mean error radius for that location class.
#  3. Finally, observations with no anisotropic errors available have their semi major and semi minor ellipse axes set to the isotropic error
#  radius and the ellipse orientation set to zero (i.e. a circle centred on the point).

#' @param object A dataframe containing Argos locations with errors
#' 
#' @param LC The name of the column containing location classes
#'
#' @param LC.values A named vector giving error radii associated with specific location classes if these are not
#' implicit in the data. This can be useful for passing approximate errors to user entered deployment locations or
#' GPS positions without error information.
#'
#' @param error.radius,minor.axis,major.axis,orientation Character vectors giving the names of the columns containing
#' isotropic error radii, anisotropic error ellipse axes and error ellipse orientation, respectively. By default
#' they are set to the names in Wildlife Computers location summary .csv files.
#
#' @return The original dataframe object with missing location error information filled in.

# -------------------------------------------------------------------------------------------------------

make_argos_errors = function(object,LC = 'Quality',LC.values = NULL ,error.radius = 'Error.Radius',
                             minor.axis = 'Error.Semi.Minor.Axis' ,major.axis = 'Error.Semi.Major.Axis',
                             orientation = 'Error.Ellipse.Orientation'){
  
  if(!is(object,'data.frame')) stop("object must be a dataframe")
  
  object = ungroup(object)  
  
  # fill in missing isotropic errors with named values....
  
 object[[LC]] = factor(object[[LC]])
  
  if(!missing(LC.values)){
    if(is.null(names(LC.values)))  stop("LC.values must be a named vector")
    object[[error.radius]] = ifelse(is.na(object[[error.radius]]),LC.values[as.character(object[[LC]])],object[[error.radius]])
  }
  
  # ... or with means of location classes if values not given
  
  missing = which(is.na(object[[error.radius]]))  
  means = tapply(object[[error.radius]],object[[LC]],mean,na.rm=T)
  object[[error.radius]][missing] = means[object[[LC]][missing]]
  
  #fill in missing anisotropic errors with isotropic error radii
  
  missing = which(is.na(object[[major.axis]]))
  object[[major.axis]][missing] = object[[error.radius]][missing]
  object[[minor.axis]][missing] = object[[error.radius]][missing]
  object[[orientation]][missing] = 0
  
  object
  
}


# ----------------------------------------------------------------------------------------------------
# fastGPS_errors
# ----------------------------------------------------------------------------------------------------

#' Returns a dataframe containing the approximate sqrt(2)*standard deviation error radii of FastlocGPS locations based on
#' number of satellites used. Values digitised from figure 1 of Dujon et al. 2014
#' <https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12286>
#'
#' @return A dataframe of error radii

# -------------------------------------------------------------------------------------------------------

fastGPS_errors <- function(){

c(
  '11' = 16.5,
  '10' = 16.5,
  '9'  = 19.1,
  '8'  = 23.3,
  '7'  = 28.4,
  '6'  = 40.0,
  '5'  = 75.0,
  '4'  = 200.0
) %>% enframe('nsats','err')
  
}


