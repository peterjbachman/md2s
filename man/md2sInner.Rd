% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/md2sInner.R
\name{md2sInner}
\alias{md2sInner}
\title{md2sInner}
\usage{
md2sInner(
  X0,
  y0,
  X.c.0 = NULL,
  X.X.0 = NULL,
  X.y.0 = NULL,
  init0 = "svd",
  tol0 = 1e-06
)
}
\arguments{
\item{X0}{double centered matrix of X from md2s passed into md2sinner}

\item{y0}{double centered matrix of X from md2s passed into md2sinner}

\item{X.c.0}{empty matrix of covariates in shared space}

\item{X.X.0}{empty matrix of covariates in X space, alone not shared with that of Y}

\item{X.y.0}{empty matrix of covariates in Y space, alone not shared with that of X}

\item{init0}{option set to singular value decomposition, although pseudoinverse is run if svd is not logical in md2sinner}

\item{tol0}{tolerance param}
}
\value{
a list with the following
\item{z}{See md2s documentation}
\item{z.X}{See md2s documentation}
\item{z.y}{See md2s documentation}
\item{w.Xs}{See md2s documentation}
\item{w.ys}{See md2s documentation}
\item{beta.z}{Coefficient from solution.See md2s documentation}
\item{proportionX}{Proportion of X in y. See md2s documentation}
}
\description{
Determines which method of solution to solve the given system. In md2s, user selects either singular value decomposition or pseudoinverse.
This function computes each method of solving. First, it runs conditionals to check whether pseudoinverse is appropriate. Then, it runs the method for solution,
solves least squares from the solution, and then computes correlation coefficients between normed data sets. It then runs log likelihood for OLS parameters.
In essance, if the system is square, solve the system with singular value decomposition. If the system is not square, compute and solve the system with pseudoinverse.
}
\author{
Peter Bachman \href{mailto:bachman.p@wustl.edu}{bachman.p@wustl.edu}, Patrick Edwards \href{mailto:edwards.p@wustl.edu}{edwards.p@wustl.edu}, and Zion Little \href{mailto:l.zion@wustl.edu}{l.zion@wustl.edu}
}
