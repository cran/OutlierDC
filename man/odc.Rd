\name{odc}
\alias{odc}
\title{Outlier detection using quantile regression for censored data}
\description{
	outlier detection algorithms using quantile regression for censored data 
}
\usage{
     odc(formula, data, 
          method = c("score", "boxplot","residual"), 
          rq.model = c("Wang", "PengHuang", "Portnoy"), 
          k_r = 1.5, k_b =1.5, h = .05)
}

\arguments{
   \item{formula}{a type of \code{Formula} object with a \code{survival} object on the left-hand side of the ~ operator and covariate terms on the right-hand side. The survival object with survival time and its censoring status is constructed by the \code{\link[survival]{Surv}} function in \code{survival} package.}
   \item{data}{a data frame with variables used in the \code{formula}. It needs at least three variables, including survival time, censoring status, and covariates. }
  \item{method}{the outlier detection method to be used. The options \code{"socre", "boxplot"}, and \code{"residual"} conduct the scoring, boxplot, and residual-based algorithm, respectively. The default algorithm is \code{"score"}.
  }
  \item{rq.model}{the type of censored quantile regression to be used for fitting. The options \code{"Wang", "Portnoy"}, and \code{"PengHuang"} conduct Wang and Wang's, Portnoy's, and Peng and Huang's censored quantile regression approaches, respectively. The default is \code{"Wang"}.
  }
  %\item{bound}{the type of bound to be used for fitting. The option \code{"both"} conduct the routine to find both upper and lower outliers, while the options \code{"UB"} and \code{"LB"} work for finding upper outliers and lower outliers, respectively.
  %}
  \item{k_r}{a value to control the tightness of cut-offs for the residual algorithm with a default value of 1.5.
  }
  \item{k_b}{a value to control the tightness of cut-offs for the boxplot algorithm with a default value of 1.5.
  }
   \item{h}{bandwidth for locally weighted censored quantile regression with a default value of 0.05.
   }
}
\details{
	The \code{odc} function conducts three outlier detection algorithms on the basis of censored quantile regression. Three outlier detection algorithms were implemented: residual-based, boxplot, and scoring algorithms. The residual-based algorithm detects outlying observations using constant scale estimates; however, it does not account for the heterogeneity of variability. When the data is extremely heterogeneous, the boxplot algorithm with censored quantile regression is more effective. The residual-based and boxplot algorithms produce cut-offs to determine whether observations are outliers. In contrast, the scoring algorithm provides the outlying magnitude or deviation of each point from the distribution of observations. Outlier detection is achieved by visualising the scores. 
}
\value{
	an object of the S4 class "OutlierDC" with the following slots:\cr
	call: evaluated function call\cr
  formula: formula to be used\cr
  raw.data: data to be used for model fitting\cr
  refined.data: the data set after removing outliers\cr
  refined.data: the data set containing outliers\cr
  coefficients: the estimated censored quantile regression coefficient matrix consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr
  fitted.mat: the censored quantile regression fitted value matrix consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr 
  score: outlying scores (scoring algorithm) or residuals (residual-based algorithm)\cr
  cutoff: estimated scale parameter for the residual-based algorithm\cr
  lower: lower fence vector used for the boxplot and scoring algorithms\cr
  upper: upper fence vector used for the boxplot and scoring algorithms\cr
  outliers: logical vector to determine which observations are outliers\cr
  n.outliers: number of outliers detected\cr
  method: outlier detection method to be used\cr
  rq.model: censored quantile regression to be used\cr
  k_r: a value to be used for the tightness of cut-offs in the residual algorithm \cr
  k_b: a value to be used for the tightness of cut-offs in the boxplot algorithm \cr
  k_s: a value to be used for the tightness of upper fence cut-offs used for the scoring algorithm with the \code{\link{update}} function\cr
  %LB: lower fence used for the scoring algorithm with the \code{\link{update}} function\cr
}
\source{
Eo S-H, Hong S-M Hong, Cho H (2014). Identification of outlying observations with quantile regression for censored data, \emph{Submitted}.

Wang HJ, Wang L (2009) Locally Weighted Censored Quantile Regression. JASA 104:1117--1128. doi: 10.1198/jasa.2009.tm08230
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{update}}
}
\examples{
  \dontrun{
    require(OutlierDC)
    # Toy example 
    data(ebd)
    # The data consists of 402 observations with 6 variables. 
    dim(ebd)
    # To show the first six observations of the dataset,
    head(ebd)
    
    #scoring algorithm
    fit <- odc(Surv(log(time), status) ~ meta, data = ebd)
    fit
    coef(fit)
    plot(fit)

    # Add upper bound for the selection of outleirs
    fit1 <- update(fit, k_s = 4)
    fit1
    plot(fit1)

    # residual-based algorithm
    fit2 <- odc(Surv(log(time), status) ~ meta, data = ebd, method = "residual", k_r = 1.5)
    fit2
    plot(fit2)
    
    # To display all of outlying observations in the fitted object
    fit2@outlier.data
    
    # boxplot algorithm
    fit3 <- odc(Surv(log(time), status) ~ meta, data = ebd, method = "boxplot", k_b = 1.5)
    fit3
    plot(fit3, ylab = "log survival times", xlab = "metastasis lymph nodes")
}
}
\keyword{models}


