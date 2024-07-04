\name{Rcpp_modelobject-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_modelobject-class}

\title{Class \code{"Rcpp_modelobject"}}
\description{
A class required when using the \code{imputedata} function with argument \code{method = "JM-DP"}
}
\section{Extends}{
Class "C++Object", directly.

All reference classes extend and inherit methods from "envRefClass".

}
\author{
Kim, Hang Joon. University of Cincinnati.
}
\note{
From the R package \code{DPImputeCont} \url{https://github.com/hang-j-kim/DPImputeCont}.
}

\examples{
showClass("Rcpp_modelobject")
}
\keyword{classes}
\keyword{internal}
\section{Fields}{
  \describe{
    \item{\code{.Mu}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.Z_vec}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.a_Phi}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.a_alpha}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.alpha}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.b_Phi}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.b_alpha}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.cube_UT_cholSigma}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.h_0}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.logpi}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.missing_flag_mat}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.vec_Phi}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.where_we_are}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{Y_mat}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{max_Y_obs}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{min_Y_obs}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{msg_level}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{vec_HyperParameters}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{finalize()}:}{ ~~ }
    \item{\code{Run(...)}:}{ ~~ }
    \item{\code{Iterate()}:}{ ~~ }
  }
}
