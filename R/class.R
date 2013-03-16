setOldClass("Formula")
setOldClass("Surv")
setClass(Class = "OutlierDC", representation(
							call = "language",
							formula = "Formula",
							raw.data = "data.frame",
							refined.data = "data.frame",
							coefficients ="data.frame",
							fitted.mat = "matrix",
							score = "vector",
							cutoff = "vector",
							lower = "vector",
							upper = "vector",
							outliers = "vector",
							n.outliers = "integer",
							method = "character",
							rq.model = "character",
							k = "numeric",
							UB = "numeric",
							LB = "numeric")
)
