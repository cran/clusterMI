% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chooseB.intern.R
\name{chooseB.intern}
\alias{chooseB.intern}
\title{Tune the number of iterations for variable selection using varselbest}
\usage{
chooseB.intern(
  res.varselbest,
  linewidth = 1,
  linetype = "dotdash",
  xlab = "B",
  ylab = "Proportion",
  graph = TRUE,
  title = NULL
)
}
\arguments{
\item{res.varselbest}{an output from the varselbest function}

\item{linewidth}{a numerical value setting the widths of lines}

\item{linetype}{what type of plot should be drawn}

\item{xlab}{a title for the x axis}

\item{ylab}{a title for the y axis}

\item{graph}{a boolean. If FALSE, no graphics are ploted. Default value is TRUE}

\item{title}{the main title}
}
\description{
Tune the number of iterations for variable selection using varselbest
}
\references{
Bar-Hen, A. and Audigier, V. An ensemble learning method for variable selection: application to high dimensional data and missing values. ArXiv e-prints <arXiv:1808.06952>
}
\keyword{internal}
