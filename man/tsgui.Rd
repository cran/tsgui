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

\details{
  Sliders and Buttons of the GUI are:
  \describe{
    \item{\sQuote{New innovations}}{The random seed is fixed when model
      parameters are changed. \sQuote{New Innovations} draws a new random
      seed}
    \item{\sQuote{Return}}{leaves the program}
    \item{\sQuote{MODELS}}{choose one of the given models to be displayed}
    \item{\sQuote{time}}{length of the time series}
    \item{\sQuote{burn.in}}{the length of the burn in period}
    \item{\sQuote{repetitions}}{number of simulations (with different seeds) shown
      in the graphics}
    \item{\sQuote{INNOVATIONS}}{the distribution families of the independent
      errors. If more than one family is chosen the graphs are put on
      top of each other.       
    }
    \item{\sQuote{green button cpy}}{
      allows to freeze the current choice of the parameter values.
      If it is pressed a new column appears with the frozen values. The
      effect becomes apparent when now the sliders are moved.
      The new column has a \code{d} for deleting that column and a
      \code{1} to switch this column with the current columns of the
      sliders. 
    }
    \item{\sQuote{MODEL PARAMETERS}}{
      The ARMA model is given by
      \deqn{Z_t - \phi_1 Z_{t-1}- \phi_2 Z_{t-2}- \phi_3 Z_{t-3}
	= \varepsilon_t + \theta_1 \varepsilon_{t-1}
	+ \theta2 \varepsilon_{t-2}
	+ \theta3 \varepsilon_{t-3}}
      \cr
      The GARCH model is given by \eqn{Z_t = \sigma_t \varepsilon_t} and
      \deqn{
	\sigma_t^2 - \beta_1 \sigma^2_{t-1} - \beta_2 \sigma^2_{t-2}
	- \beta_3 \sigma^2_{t-3} =
	\alpha_0 + \alpha_1 Z_{t-1}^2 + \alpha_2 Z_{t-2}^2 
      }
    }
    \item{\sQuote{INNOVATION PARAMETERS}}{the parameters corresponding to the
      distributions of the INNOVATIONS}
    \item{\sQuote{STARTING VALUES}}{starting values for the burn in period}
  }

  The argument \code{user} allows to define
  ARCH and GARCH models with different numbers of parameters,
  see the examples in \command{\link{ARMA.model}}
  and \command{\link{GARCH.model}}.

  The argument \code{user} also allows to define one's own model.
  The example of such a definition is returned by \command{\link{ARMA.model}}.
  While  \code{ARMA.model()} returns an overview and
  \code{str(ARMA.model(), give.attr=TRUE)} gives some more details,
  the code \code{print(ARMA.model(), details=TRUE)} gives the full
  definition of a time series:
  \describe{
    \item{\code{time}, \code{burn.in}, \code{repetitions}}{the default values of
      the  respective sliders, see above}
    \item{\code{*.min}, \code{*.max}}{The lower and upper bound of the
  slider values.} 
    \item{\code{phi}}{named vector; the names will be the titles of the
      sliders; phi contains all model parameters}
    \item{\code{distr}}{note that all the three distribution functions have access
      to the values of all three INNOVATION PARAMETERS in the gui, since
      \code{param} is a vector of three elements}
    \item{\code{distr.param}}{named vector of parameters. It could be also a
      shorter or longer vector than the standard vector of length 3}
    \item{\code{distr.show}}{at least one of the components must be \code{TRUE},
      indicating the default distribution}
    \item{\code{titles}}{a list of 2 functions
      \describe{
	\item{\code{titles[[1]]}}{
	  returns the title of the model in dependence of the current
	  values of the model parameters
	}
	\item{\code{titles[[2]]}}{
	  returns a string of additional information on the model
	}
      }
    }
    \item{\code{update.function}}{
      for given \code{innovations}, i.e. realisations of the errors
      \eqn{\varepsilon_i} including burn in,
      the parameter values \eqn{\phi} of the model and
      for given starting values \code{start.ts} the function calculates
      the time series, including burn in. Note that \code{start.ts}
      is a matrix of \code{length{param}} rows and \eqn{n r}
      columns, where \eqn{n} is the number of chosen
      distribution families and  \eqn{r} is the number of repetitions.
      Similar, \code{innovations} is also 
    }
    \item{\code{starting.ts}}{the starting values of the process in the burn-in
      period of the simulation
    }
    \item{\code{starting.matrix.ts}}{
      This function turns the starting values \code{starting.ts}
      into a matrix where is essentially the
      starting values \code{starting.ts} are written in each column.
      The number of columns are the number of chosen distribution family
      time the number of repetitions. Note that in the matrix the
      repetitions are grouped to gether, i.e. if \eqn{r} is the number
      of repetitions the first \eqn{r} rows refer to the \eqn{r}
      repetitions of the first chosen distribution.

      The function is introduced to allow also for transformations of
      \code{starting.ts} for different distributions.
    }
    \item{\code{starting.innovations}}{the starting values of the innovations
      in the burn-in 
      period of the simulation
    }
    \item{\code{starting.matrix.innovations}}{
     function analogue to \code{starting.matrix.ts}
    }
      
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
