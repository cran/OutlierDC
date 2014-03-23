\name{OutlierDC-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierDC-class}
\alias{coef, OutlierDC-class}
\title{ "OutlierDC" class}
\description{
	The \code{OutlierDC} class presents outlier detection algorithms for censored data. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OutlierDC")}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{evaluated function call}
    \item{\code{formula}:}{formula to be used with the type of \code{"Formula"} }
    \item{\code{raw.data}:}{data to be used with the type of \code{"data.frame"} }
    \item{\code{refined.data}:}{the data set after removing outliers}
    \item{\code{outlier.data}:}{the data set containing outliers}
    \item{\code{coefficients}:}{the estimated censored quantile regression coefficient matrix}
    \item{\code{fitted.mat}:}{the censored quantile regression fitted value matrix with the type of \code{"matrix"} }
    \item{\code{score}:}{outlying scores (scoring algorithm) or residuals (residual-based algorithm)}
    \item{\code{cutoff}:}{estimated scale parameter for the residual-based algorithm}
    \item{\code{lower}:}{lower fence vector used for the boxplot and scoring algorithms with the type of \code{"vector"} }
    \item{\code{upper}:}{upper fence vector used for the boxplot and scoring algorithms with the type of \code{"vector"} }
    \item{\code{outliers}:}{logical vector to determine which observations are outliers}
    \item{\code{n.outliers}:}{number of outliers to be used. The object of class \code{"integer"}. }
    \item{\code{method}:}{outlier detection method to be used}
    \item{\code{rq.model}:}{censored quantile regression to be fitted}
    \item{\code{k_r}:}{a value to be used for the tightness of cut-offs in the residual-based algorithm}
    \item{\code{k_b}:}{a value to be used for the tightness of cut-offs in the boxplot algorithm}
    \item{\code{bound}:}{type of fence to be used in the model fittind}
    \item{\code{k_s}:}{a value to be used for the tightness of upper fence cut-offs used for the scoring algorithm with \code{update} function}
    %\item{\code{LB}:}{lower fence used for the scoring algorithm with \code{update} function. The object of class \code{"numeric"} }
    }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "OutlierDC")}: Print the coefficient matrix of censored quantile regression to be used. See \code{\link{coef}}. }
    \item{plot}{\code{signature(x = "OutlierDC", y = "missing")}: See \code{\link{plot}}. }
    \item{show}{\code{signature(object = "OutlierDC")}: See \code{\link{show}}. }
    \item{update}{\code{signature(object = "OutlierDC")}: Update the fitted object to find outliers in scoring algorithm. See \code{\link{update}}. }
	 }
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{update}}
}
\examples{
showClass("OutlierDC")
}
\keyword{class}
%end