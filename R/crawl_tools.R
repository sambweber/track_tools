
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
