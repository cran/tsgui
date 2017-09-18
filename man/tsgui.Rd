\name{tsgui}
\alias{tsgui}
\title{
Gui for times series
}
\description{
  This GUI allows real-time simulations of ARMA and GARCH processes.
  Further models can be added by the user.
}
\usage{
tsgui(user = NULL, wait = 1000, included.models =c("ARMA", "GARCH"))	  
}
\arguments{
  \item{user}{a further class of model (or a list of such models),
    additional to ARCH(3,3) and GARCH(3,3)
  }
  \item{wait}{non-negative integer. Users will usually not change this value.
    if the value is greater than 0, the control over the command line of R
    is returned to the user only after the tcltk window has been left.
    \command{tsgui} checks every \code{wait} /  1000 seconds whether the
    tcltk window has been left.
  }
  \item{included.models}{\code{NULL}
    or string vector giving the internally defined
    models that shall be available. 
  }
}

\value{
  A list with the latest chosen paremeters is returned.
}
%\references{}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{\command{\link{ARMA.model}}, \command{\link{GARCH.model}}}

\examples{  
  if (interactive()) tsgui()
}

\keyword{ts}
\keyword{models}
