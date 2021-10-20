\name{tensorIMMAM-package}
\alias{tensorIMMAM-package}
\alias{tensorIMMAM}
\docType{package}
\title{
  Estimation for integrative multi-view multivariate additive model incorporating tensor decomposition
}
\description{
  For an integrative multi-view multivariate additive model (IMMAM), the B-splines are applied to approximate the component functions. We treat the coefficients as multiple third-order tensors or even a fourth-order tensor in some specific conditions (called structural IMMAM). With the tensor low-rankness, the Tucker decomposition and group sparse penalty (lasso, mcp or scad) reduce the number of parameters. An ADMM-type algorithm based on the coordinate descent strategy is used to estimate the core tensors and factor matrices, and further additive functions.

}
\details{
  This pacakge includes four main functions and four generating functions. \code{immam3} and \code{immam3.dr} yield the estimator of IMMAM. The difference is the former requires the fixed ranks and the latter can search the optimal ranks and regularization parameter simultaneously by BIC or CV method. \code{immam3.sim.fbs} and \code{immam3.sim.fsin} generate data of scenario I and II respectively. Scenario I assumes that the true functions are exactly residing in the space of B-spline basis functions. Scenario II assumes that the true functions are some linear combination of \eqn{sin(2\pi x)} and \eqn{cos(\pi x)}. Similarly, \code{immam4} and \code{immam4.dr} yield the estimator of structural IMMAM. \code{immam4.sim.fbs} and \code{immam4.sim.fsin} are two generating functions of scenario I and II, which have the same assumptions as IMMAM.
}
\author{
Xu Liu

Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
}
\references{
  Estimation for integrative multi-view multivariate additive model incorporating tensor decomposition.
}
\keyword{ Group sparsity; Tensor low-rankness; Tucker decomposition; Multivariate additive regression; Multi-view data. }
