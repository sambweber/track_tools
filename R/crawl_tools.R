
add_anisotopic_errors <- function(obj, Major = "Error.Semi.Major.Axis", Minor = "Error.Semi.Minor.Axis",
Orientation = "Error.Ellipse.Orientation"){

cbind(obj, 
      crawl::argosDiag2Cov(
        obj[[Major]], 
        obj[[Minor]], 
        obj[[Orientation]])
      )

}
