\name{update}
\docType{methods}
\alias{update}
\alias{update,OutlierDC-method}
\title{ Update a scoring algorithm.}
\description{
This function updates a scoring algorithm using upper and/or lower fences. Using the call stored in the object, the \code{update} function declares outlying observatoins based on the QQ plot. \code{k_s} is used to set the upper cut-off bound.
% and \code{LB} is used to set the lower cut-off bound. 
}
\usage{
    \S4method{update}{OutlierDC}(object, k_s = NA, LB = NA)
}
\arguments{
    \item{object}{
        fitted model object of class \code{\linkS4class{OutlierDC}}.
        }
    \item{k_s}{ 
		cut-off value for the upper fence
	}
    \item{LB}{  
		cut-off value for the lower fence
	}
}
\details{
  This function is a generic function called \code{update} for the S4 class \code{OutlierDC}. Cut-off bounds are added to find outliers on the normal QQ plot. 
}
\seealso{
    \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
}
\keyword{methods}
