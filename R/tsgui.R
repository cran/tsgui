
print.tsgui <- function(x,...) {
  n <- names(x) 
  idx <- sapply(strsplit(n, ".mi"), length) > 1 |
         sapply(strsplit(n, ".ma"), length) > 1 |
         sapply(strsplit(n, ".intege"), length) > 1 
  name <- attr(x, "model")
  if (!is.null(name)) cat("'", name, "' is a ")
  class(x) <- NULL
  utils::str(x[!idx], ..., give.attr=FALSE) ##
}

str.tsgui <- function(object, ..., give.attr=FALSE) {
  name <- attr(object, "model")
  if (!is.null(name)) cat("'", name, "' is a ")
  class(object) <- NULL
  utils::str(object,..., give.attr=give.attr) ##
}


ARMA.model <- function(p.max=3, q.max=3) {
  nphi <- p.max
  ntheta <- q.max
  stopifnot(nphi >= 0, ntheta >= 0, nphi + ntheta > 0)
  arma.phi <- c(if (nphi > 0) c(0.5, rep(0, nphi - 1)),
		if (ntheta > 0) c(0.5, rep(0, ntheta -1))
		)
  names(arma.phi) <- c(if (nphi > 0) paste("phi", 1:nphi),
		       if (ntheta > 0) paste("theta", 1:ntheta))
  starting.ts <- rep(0, nphi)
  if (nphi > 0) names(starting.ts) <- paste("process value at", (1-nphi):0)
  starting.innovations <- rep(0, ntheta)
  if (ntheta > 0)
    names(starting.innovations) <- paste("innovation at", (1-ntheta):0)
  
  ARMA <-
  list(
    time = 2000,
    burn.in=20L,
    repetitions = 1L,

    phi = arma.phi,
    phi.min = c(rep(-2, nphi),
		rep(-2, ntheta)),
    phi.max = c(rep(2, nphi),
		rep(2, ntheta)),

    distr = c(norm=function(p, param) qnorm(p, 0, param[1]),
	      exponential=function(p, param) qexp(p, rate=1/param[2])-param[2],
	      cauchy=function(p, param) qcauchy(p, scale=param[3])
	      ),
    distr.param = c("normal sigma" = 1,
                    "exponential inverse rate"=0.5,
                    "cauchy scale"=0.01
                    ),
    distr.param.min = c(0.01, 0.2, 0.001),
    distr.param.max = c(10, 10, .1),
    distr.show = c("normal"=TRUE, "exponential"=FALSE, "cauchy"=FALSE),
    
    titles = list(function(phi) {
      theta <- phi[-1:-nphi]
      i <- nphi
      while (i>0 && phi[i] == 0) i <- i - 1;
      j <- ntheta
      while (j>0 && theta[j] == 0) j <- j - 1;
      if (j == 0) {
	if (i == 0) "white noise"
	else paste("AR(", i, ")", sep="")
      } else {
	if (i==0) paste("MA(", j, ")", sep="")
	else paste("ARMA(", i, ", ", j, ")", sep="")
      }
    },
    function(phi) {
      i <- nphi
      while (i>0 && phi[i] == 0) i <- i - 1;
      if (i == 0) "stationary process"
      else {
	r <- polyroot(c(1, -phi[1:i]))
	r.re <- round(Re(r) * 1000) / 1000
	r.im <- round(Im(r) * 1000) / 1000
	r.txt <- format(r.re + r.im * 1i, dig=2)
	if (any(r.im==0)) r.txt[r.im==0] <- format(r.re[r.im==0], dig=2)
	outside <- all(Mod(r)>1)
	paste(if (!outside) "not all ",
	      "roots of phi (", paste(r.txt, collapse=";"),
	      ") outside unit circle, i.e. ", if (!outside) "non-",
	      "stationary", sep="")
      }
    }),

    update.function = function(innovations, phi, start.ts) {
      theta <- phi[-1:-nphi]
      phi <- phi[1:nphi]
      n.innovations <- nrow(innovations) - length(theta)
      X <- rbind(start.ts,
		 matrix(NA, ncol=ncol(start.ts), nrow=n.innovations))
      for (i in 1:n.innovations) {
	segX <- i + nphi
	segI <- i + ntheta
	u <- innovations[segI, ]
	if (nphi > 0) u <- u + (phi %*% X[segX + (-1 : - nphi), ])
	if (ntheta > 0) u <- u+(theta %*% innovations[segI + (-1 : - ntheta), ])
	X[segX, ] <- u
      }
      X
    },

    starting.ts = starting.ts,
    starting.ts.min = rep(-10, nphi),
    starting.ts.max = rep(10, nphi),
    starting.matrix.ts = function(repetitions, n.distr, param) {
      matrix(param, nrow=length(param), ncol=n.distr * repetitions)
    },
    starting.innovations = starting.innovations,
    starting.innovations.min = rep(-10, ntheta),
    starting.innovations.max = rep(10, ntheta),
    starting.matrix.innovations = function(repetitions, n.distr, param) {
      matrix(param, nrow=length(param), ncol=n.distr * repetitions)
    }
    )
  attr(ARMA, "multivariate") <- 1
  attr(ARMA, "model") <- "ARMA"
  class(ARMA) <- "tsgui"
  ARMA
}


GARCH.model <- function(a.max=3, b.max=3) {
  na <- a.max
  nb <- b.max
  stopifnot(na > 1, nb >= 0)
  garch.a <- c(c(0.005, 0.1, rep(0, na - 2)),
	       if (nb > 0) c(0.5, rep(0, nb -1)))

  names(garch.a) <- c(paste("a", 0:(na-1)),  ## x
		      if (nb > 0) paste("b", 1:nb)) ## sigma^2
  starting.ts <- rep(0, nb)
  if (nb > 0) names(starting.ts) <- paste("process value at", (1-nb):0)
  starting.innovations <- rep(0, na)
  names(starting.innovations) <- paste("innovation at", (1-na):0)
  
  GARCH <-
  list(
    time = 2000,
    burn.in=20L,
    repetitions = 1L,

    phi = garch.a,
    phi.min = c(rep(0, na),
		rep(0, nb)),
    phi.max = c(0.01, rep(10, na-1) / na,
		if (nb > 0) c(3, rep(1.5, nb - 1)) / b.max),
#
    ## button whether they should be used
    distr = c(norm=function(p, param) qnorm(p, 0, param[1]),
	      exponential=function(p, param) qexp(p, rate=1/param[2])-param[2],
	      cauchy=function(p, param) qcauchy(p, scale=param[3])
	      ),
    distr.param = c("normal sigma" = 1,
                    "exponential inverse rate"=1,
                    "cauchy scale"=0.01
                    ),
    distr.param.min = c(0.01, 0.2, 0.001),
    distr.param.max = c(10, 10, .1),
    distr.show = c("normal"=TRUE, "exponential"=FALSE, "cauchy"=FALSE),
    
    titles = list(function(a) {
      b <- a[-1:-na]
      i <- na
      while (i>0 && a[i] == 0) i <- i - 1;
      j <- nb
      while (j>0 && b[j] == 0) j <- j - 1;
      if (j == 0) {
	if (i <= 1) c("white noise", "constant")
	else c(paste("ARCH(", i-1, ")", sep=""), "sigma")
      } else {
	if (i==0) c("independent variables", paste("deterministic", sep=""))
	else c(paste("GARCH(", i-1, ", ", j,")", sep=""), "sigma")
      }
    }),

    update.function = function(innovations, a, start.ts) {
      b <- a[-1:-na]
      a0 <- a[1]
      a <- a[2:na]
      n.innovations <- nrow(innovations) - na + 1
      X <- innovations
      S2 <- rbind(start.ts,
		  matrix(NA, ncol=ncol(start.ts), nrow=n.innovations))
      for (i in 1:n.innovations) {
	segS <- i + nb
	segX <- i + na - 1
	u <- a0 + (a %*% X[segX + 1 + (-2 : -na), ]^2)
##	if (i == 1) Print(u, a, a0, nb)
	if (nb > 0) u <- u + (b %*% S2[segS + (-1 : -nb), ])
	S2[segS, ] <- u	
	X[segX, ] <- sqrt(S2[segS, ]) * X[segX, ] 
      }
      list(X=X, sigma=sqrt(S2))
    },

    starting.ts = starting.ts,
    starting.ts.min = rep(0, nb),
    starting.ts.max = rep(10, nb),
    starting.matrix.ts = function(repetitions, n.distr, param) {
      matrix(param, nrow=length(param), ncol=n.distr * repetitions)
    },
    starting.innovations = starting.innovations,
    starting.innovations.min = rep(-10, na),
    starting.innovations.max = rep(10, na),
    starting.matrix.innovations = function(repetitions, n.distr, param) {
      matrix(param[-1], nrow=length(param)-1, ncol=n.distr * repetitions)
    }
    )
  attr(GARCH, "multivariate") <- 2
  attr(GARCH, "model") <- "GARCH"
  class(GARCH) <- "tsgui"
  GARCH
}


tsgui <- function(user = NULL,
		  wait = 1000,		    
		  included.models =c("ARMA", "GARCH")
		  ) {
  if (length(included.models) > 0) {
    included.models <- match.arg(included.models, several.ok = TRUE)
    all.included.models <- eval(as.list(args(tsgui))$included.models)
  }
  
#### see https://www.tcl.tk/man/tcl/TkCmd/label.htm
  
  if (!interactive()) {
    #warning("'tsgui' can be used only in an interactive mode")
    #return(NULL)
  }
  wait <- as.integer(wait)
  Env <- if (wait >= 0) environment() else .GlobalEnv
  if (exists(".tsgui.exit", .GlobalEnv))
    rm(".tsgui.exit", envir=.GlobalEnv)

  order.show <- list("time",
		     "burn.in",
		     "repetitions", ## how often is simulation study repeated?
		     ##     Then it gives the time instance at which it is shown
		     "GLOBAL:INNOVATIONS" = c(1),
		     "distr.show",
		     "MODEL PARAMETERS" = c(1, 1, 1),
		     "phi",
		     "INNOVATION PARAMETERS" = c(2, 2),
		     "distr.param",
		     "STARTING VALUES" = c(3, 3),
		     "starting.ts",
		     "starting.innovations"
		     )
  
  constraints <- list(
      time.min = 2L, ## fraction of totaltime
      time.max = 10000L,
      burn.in.min = 0L, ## fraction of totaltime
      burn.in.max = 50L,
      repetitions.min = 1L,
      repetitions.max = 5L,
      distr.show.min = FALSE,
      distr.show.max = TRUE,
      phi.min = -1,
      phi.max = 1,
      distr.param.min = 0,
      distr.param.max = 10,
      starting.ts.min = -10,
      starting.ts.max = 10,
      starting.innovations.min = -10,
      starting.innovations.max = 10
    )

  if (length(user) > 0) {   
    if (!is.list(user[[1]])) user <- list(user = user)
    modelclasses <- user
  } else modelclasses <-list()

  for (i in included.models)
    modelclasses[[i]] <- switch(pmatch(i, all.included.models),
			   "1" = ARMA.model(),
			   "2" = GARCH.model()
			   )

  tsgui.intern(currentClass = 1, 
	      order.show = order.show,
	      constraints = constraints,
              parent.ev = Env,
	      modelclasses = modelclasses)

  if (wait >= 0) {
    while (!exists(".tsgui.exit", envir=Env))
      RandomFieldsUtils::sleep.micro(wait)
    res <- get(".tsgui.exit", envir=Env)
    rm(".tsgui.exit", envir=Env)
    if (is.null(res)) return(res)
    for (i in 1:length(res)) {
      if (is.function(res[[i]])) environment(res[[i]]) <- .GlobalEnv
    }
    class(res) <- "tsgui"
    return(res)
  } else invisible(NULL)
}

tsgui.intern <- function(currentClass,
			order.show,
			constraints,
			parent.ev=NULL,		
			colour=c("black",  "darkred", "darkblue",
				      "orange", "darkgreen"),
#			subcol=c("grey",  "red", "lightblue",
#				     "yellow", "green", "orange"),
			fgcol=c("white", "white", "white", "white", "white"),
			modelclasses = NULL
			) {
 # set.mixcol <- c(colour[1], subcol[-1])
  max.sets <- length(colour)
  numberSteps <- 256L ## even 
  image.rowspan <- 15L
  image.colspan <- 5L
  fst.col <- 1L
  half.col <- fst.col + 0L
  snd.col <- half.col + 2L
  trd.col <- snd.col + c(3L, 5L)
  fst.row <- 1L
  snd.row <- 17L
  trd.row <- 32L ## position of sliders of global variables
  col.sl <- fst.col + 2 * image.colspan + 1L ## position of slides for model
  row.last <- image.rowspan + 4L
  row.sl <- 1L # giving the position of the first label
  width.entry <- 4L
  width.slider <- 15L ## 18L
  length.slider <- 170L
  length.slider.main <- 130L
  plothscale <- 0.8    # Horizontal scaling
  plotvscale <- 0.8    # Vertical scaling
  importance <- rep(2, 10)
  fg <- c("black", "gray25", "gray40", "gray70")
  col.alert <- "darkred"

  ENVIR <- environment()

  L <- (wait.simulation <- innovations <- Len <- isinteger <- minall 
	<- maxall <- islogical <- imgLU <- sumLen <- strictpos
	<- currentNrVariab <- tt <- sets <- buttonNewSimu <- classTitle
	<- col_bg <- copyFctn <- buttonReturn
	<- NULL)
  
  tkDestroy <- tcltk::tkdestroy
  tkValue <- tcltk::tclvalue
  "tkValue<-" <- do.call("::", list("tcltk", "tclvalue<-"))
  tkLabel <- tcltk::tklabel
  tkEntry <- tcltk::tkentry
  tkScale <- tcltk::tkscale
  tkBind <- tcltk::tkbind
  tkGridConf <- tcltk::tkgrid.configure
  tkVar <- tcltk::tclVar
  tkGrid <- tcltk::tkgrid
  tkPlot <- tkrplot::tkrplot
  tkCheckbutton <- tcltk::tkcheckbutton
  tkCheckbutton <- tcltk::tkcheckbutton
  tkConfigure <- tcltk::tkconfigure
  Tcl <- tcltk::tcl
  tkRreplot <- tkrplot::tkrreplot
  tkButton <- tcltk::tkbutton
  tkDelete <- tcltk::tkgrid.forget

  non.show.loc <- which(names(order.show) != "")
  titles.val <- order.show[non.show.loc]
  idx <- sapply(titles.val, is.logical)
  do.not.show <- names(titles.val[idx])
  titles.loc <- non.show.loc[!idx]

  titles.val <- order.show[titles.loc]
  rawtitles <- strsplit(names(titles.val), "GLOBAL:")
  titles <- sapply(rawtitles, function(x) if (length(x) == 1) x else x[2])
  globaltitles <- sum(rawtitles != titles)
  globalParams <- titles.loc[1 + globaltitles] - globaltitles - 1

  titles.guiloc <- c(titles.loc, 9999)
  for (i in 1:length(titles.loc))
    titles.guiloc[i] <- titles.guiloc[i] - sum(non.show.loc < titles.loc[i])
  
  Names <- unlist(order.show[-non.show.loc])
  min.names <- paste(Names, "min", sep=".")
  max.names <- paste(Names, "max", sep=".")
  integer.names <- paste(Names, "integer", sep=".")
  all.names <- c(Names, do.not.show)

  Laengen <- sapply(modelclasses, function(x) sapply(x[all.names], length))
  row.names(Laengen) <- all.names
  orig.constraints <- rep(list(constraints), length(modelclasses))
  
  multivariate <- sapply(modelclasses, function(x) {
    v <- attr(x, "multivariate")
    if (length(v) == 0) 1 else v
  })
				       
  storedvariables <- c("L", "Len", "minall", "maxall", "isinteger",
		       "islogical", "strictpos")
  
  setClassLists <- function(s) {
    assign("currentClass", s, envir = ENVIR)
    assign("currentNrVariab", multivariate[s], envir=ENVIR)

    assign("Len", apply(Laengen, 1, max), envir=ENVIR)
    ##
    assign("Len", Laengen[, s], envir=ENVIR)
    assign("sumLen", sum(Len), envir=ENVIR)

    constraints <- orig.constraints[[s]]
    for (i in 1:length(Len)) {
      m <- max(Len[i], length(constraints[[min.names[i]]]))
      constraints[[min.names[i]]] <- rep(constraints[[min.names[i]]],
					 length.out=m)
      constraints[[max.names[i]]] <- rep(constraints[[max.names[i]]],
					 length.out=m)
    }
    assign("constraints", constraints, envir=ENVIR)

    Lorig <- modelclasses[[s]]
    for (n in all.names)
      Lorig[[n]] <- c(Lorig[[n]], rep(NaN, Len[n] - length(Lorig[[n]])))
    modelclasses[[s]] <- Lorig
    assign("modelclass", envir=ENVIR, modelclasses)

    L <- Lorig[all.names]
    Len <- sapply(L, length)
    names(Len) <- all.names

    L[do.not.show] <- Lorig[do.not.show]
    
    is_fctn <- sapply(Lorig, function(x) is.function(x) ||
                      (is.list(x) && is.function(x[[1]])))
    
    fctnNames <- names(Lorig)[is_fctn]
    L[fctnNames] <- Lorig[fctnNames]

    maxall <- constraints[max.names]
    tmp <- Lorig[max.names]
    tmp <- tmp[sapply(tmp, length) > 0]
    for (n in names(tmp)) { ## wegen strictpos
      x <- c(tmp[[n]], rep(if (is.logical(tmp[[n]])) FALSE else -1,
			   length(maxall[[n]]) - length(tmp[[n]])))      
      maxall[[n]] <- x
    }
    ## maxall[names(tmp)] <- tmp

    val <- titles.val
    tmp <- Lorig[titles]
    tmp <- tmp[sapply(tmp, length) > 0]
    val[names(tmp)] <- tmp
    assign("titles.val", val, envir=ENVIR)
  
    minall <- constraints[min.names]
    islogical <- lapply(minall, function(x) rep(is.logical(x), length(x)))
    for (i in 1:length(Names))
      if ( #Len[i] == 0 ||
	  length(islogical[[i]]) < Len[i] || any(is.na(islogical[[i]]))) {
	stop("'", Names[i], "' has an incomplete definition.\nSee 'order.show', 'constraints' and the model definition itself.\n")
      }
    isinteger <- lapply(minall, function(x) rep(is.integer(x), length(x)))
    tmp <- Lorig[min.names]
    tmp <- tmp[sapply(tmp, length) > 0]
    for (n in names(tmp)) { ## wegen strictpos
      x <- c(tmp[[n]], rep(if (is.logical(tmp[[n]])) FALSE else -1,
			   length(minall[[n]]) - length(tmp[[n]])))
      minall[[n]] <- x
    }
 
    names(isinteger) <- integer.names
    for (cond in list(constraints, Lorig)) {
      tmp <- cond[integer.names]
      names(tmp) <- integer.names
      tmp.idx <- sapply(tmp, length) > 0
      isinteger[tmp.idx] <- tmp[tmp.idx]
    }

    strictpos <- lapply(minall, function(m) m > 0)
    names(strictpos) <- Names

    for (v in storedvariables) assign(v, get(v), envir=ENVIR)
   
  }


  ERROR <- function(txt) {
    cat(txt, "\n")
    OnReturn()
  }

  Cat <- function(...) cat(...)
  Cat <- function(...) {}
 
  isScalar <- function(i) Len[i] == 1
  noeffect <- "next slider has no effect"
  labName <- function(i, j) {
    n0 <- names(L[[i]])
    if (!is.finite(L[[Names[i]]][j])) return(noeffect)
    if (length(n0) > 0) n0[j] else 
    paste(Names[i], if (!isScalar(i)) paste("[", j, "]", sep=""), sep="")
  }
  labObj <- function(i, j) paste(Names[i], "_", j, sep="")
  basename <- function(name, j, set=1)
     paste(name, set, if (!missing(j) && length(j)>0) j else 1, sep="")

  BaseName <- function(i, j, set=1) basename(Names[i], j, set)
  slValue <- function(i, j, set=1)
    paste("sl_", BaseName(i, j, set), sep="")
  entryValue <- function(i, j, set=1) 
    paste("entry_", sep="",
	  if (is.numeric(i)) BaseName(i, j, set) else basename(i, j, set))
  slWidget <- function(i, j, set=1)
    paste("wiS_", BaseName(i, j, set), sep="")
  entryWidget <- function(i, j, set=1)
    paste("wiE_", BaseName(i, j, set), sep="")
  buttonWidget <- function(i, j, set=1)
    paste("wiB_", BaseName(i, j, set), sep="")
  Value <- function(name, j, set=1) { ## sehr langsam !!
    prefix <- "entry_"
    if (missing(j)) {
      len <- Len[name]
      if (len == 0) return(logical(0))
      ans <- numeric(len)
      variab <- paste(prefix, basename(name, 1:len, set), sep="")
      for (j in 1:len) {
	ans[j] <- tkValue(get(variab[j], envir=ENVIR))
      }
    } else {
      variab <- paste(prefix, basename(name, j, set), sep="")
      ans <- tkValue(get(variab, envir=ENVIR))
    }
    as.numeric(ans)
  }
  simuValue <- function(set) paste("simu_", set, sep="")
  deleteFctn <- function(set) paste("delete_", set, sep="")
  firstFctn <- function(set) paste("first_", set, sep="")
  classButton <- function(s) paste("buttonClass", s, sep="")

  entryAssign <- function(i, j, set=1, value) { ## loest EntryChanges aus
    entry <- get(entryValue(i, j, set=set), envir=ENVIR)
    if (isinteger[[i]][j]) value <- as.integer(round(value))
    tkValue(entry) <- value
    assign(entryValue(i, j, set=set), entry, envir=ENVIR)
  }

  slAssign <- function(i, j, value) { ## loest SliderChanges aus
    sl <- get(slValue(i, j), envir=ENVIR)
    tkValue(sl) <- value
    assign(slValue(i, j), sl, envir=ENVIR)
  }
 
  
  GetL <- function(set=1) {
    L <- modelclasses[[currentClass]]
    for (i in 1:length(all.names)) {      
      if (Len[i] > 0) 
	for (j in 1:Len[i]) {
	  if (is.finite(L[[ all.names[i] ]][j])) {
	    L[[all.names[i]]][j] <- 
	    as.numeric(tkValue(get(entryValue(all.names[i],j, set),
				   envir=ENVIR)))
	  }
	}
    }
    L[names(minall)] <- minall
    L[names(maxall)] <- maxall

    modelclasses[[currentClass]] <- L
    assign("modelclasses", modelclasses, envir=ENVIR)
    orig.constraints[[currentClass]] <- constraints
    assign("orig.constraints", orig.constraints, envir=ENVIR)
    
    L[names(isinteger)] <- isinteger
    L[names(islogical)] <- islogical
    attr(L, "model") <- names(modelclasses)[[currentClass]]
    return(L)
  }

  first <- function(set, changeEntry=TRUE) { #, do.simu=TRUE, do.plot=TRUE) {
   ## assign("wait.simulation", c(first=sumLen - 1), envir=ENVIR)
    simu <- get(simuValue(set))
    for (v in storedvariables) assign(v, simu[[v]], envir=ENVIR)

    for (i in 1:length(Names)) {
      if (Len[i] > 0)
	for (j in 1:Len[i]) {
	  a <-  as.numeric(tkValue(get(entryValue(i, j), envir=ENVIR)))
	  b <-  as.numeric(tkValue(get(entryValue(i, j, set), envir=ENVIR)))
	  entryAssign(i, j, set=set, value=a)
	  entryAssign(i, j, value=b)
	  if (changeEntry) EntryChangesShort(i, j, value=b)
	}
    }
    assign("wait.simulation", c(first=0), envir=ENVIR)
    entryAssign(i, j, value=b)

    a <- get(simuValue(1), envir=ENVIR)
    for (v in storedvariables) a[[v]] <- get(v, envir=ENVIR)
    b <- get(simuValue(set), envir=ENVIR)
    assign(simuValue(set), a, envir=ENVIR)
    assign(simuValue(1), b, envir=ENVIR)
    for (v in storedvariables) assign(v, b[[v]], envir=ENVIR)
    setClass()
  }
  
  setClass <- function() {
    starting.values(sets=1)    
    start_simu(update="time")
    assign("wait.simulation", c(position=sumLen -1), envir=ENVIR)
    position()
    tkRreplot(imgLU)     
  }

  
  copy <- function(...) {
    sets <- get("sets", envir=ENVIR)
    set <- which(!sets)
    if (length(set) == 0) return()
    set <- set[1]
    sets[set] <- TRUE
    assign("sets", sets, envir=ENVIR)

    for (i in 1:length(Names)) {
      if (Len[i] > 0)
	for (j in 1:Len[i])
	  entryAssign(i, j, set,
		      as.numeric(tkValue(get(entryValue(i, j), envir=ENVIR))))
    }
    simu <- get(simuValue(1), envir=ENVIR)
    for (v in storedvariables) simu[[v]] <- get(v, envir=ENVIR)

    assign(simuValue(set), simu, envir=ENVIR)

    assign("wait.simulation", c(copy=0), envir=ENVIR)
    position()
    tkRreplot(imgLU)
  }

  
  delete <- function(set) {
    stopifnot(set > 1)
    sets <- get("sets", envir=ENVIR)
    sets[set] <- FALSE
    assign("sets", sets, envir=ENVIR)
    tkDelete(get(deleteFctn(set), envir=ENVIR))
    tkDelete(get(firstFctn(set), envir=ENVIR))

    for (i in 3:length(Names)) {
      if (Len[i]>0)
	for (j in 1:Len[i]) {
	  widget <- (if (islogical[[i]][j]) buttonWidget(i, j, set)
		     else entryWidget(i, j, set))
	  tkDelete(get(widget, envir=ENVIR))
	}
    }
    assign("wait.simulation", c(delete=0), envir=ENVIR)
    position()
    tkRreplot(imgLU)
  }
 
   
  EntryChangesShort <- function(i, j, set=1, factor=3, value) {
    if (islogical[[i]][j]) return;
    is.int <- isinteger[[i]][j]
    name <- Names[i]
    if (is.int) value <- as.integer(base::round(value))
    s <- strictpos[[i]][j]

    if (missing(value)) stop("missing value in EntryChangesShort")
    
    if (set == 1) {
      if (value == 0) {
	from <- minall[[i]][j]
	to <- maxall[[i]][j]
      } else if (!is.logical(minall[[i]][j])) {	
	from <- if (value < 0) value * factor else value / factor
	to <- if (value > 0) value * factor else value / factor
	minall[[i]][j] <<- from
 	maxall[[i]][j] <<- to
      }

      if (s) {
	value <- sqrt(value)
	from <- sqrt(from)
	to <- sqrt(to)
      }
     					      
      tkConfigure(get(slWidget(i, j, set=set), envir=ENVIR),
      		  to=to, from=from,
		  resolution=if (is.int) -1 else (to-from)/numberSteps)
      slAssign(i, j, value=value) ## zwingend nach tkconfigure
   #   v <- as.numeric(tkValue(get(slValue(i,j,set),envir=ENVIR)))      
  #    if (s) v <- v^2
  #    if (isinteger[[i]][j]) v <- as.integer(base::round(v))
   #   entryAssign(i, j, set, value=v)
    }
  }

  EntryChanges <- function(i, j, set=1, factor=3) {
    value <- as.numeric(tkValue(get(entryValue(i, j, set=set), envir=ENVIR)))
    EntryChangesShort(i, j, set=set, factor=factor, value=value)
 
    if (set != 1) {
      for (i in 1:2) first(set, changeEntry=FALSE)
    }

    return(NULL)
  }

 
  SliderChanges <- function(i, j, set=1) {

    Evalue <- value  <-
       as.numeric(tkValue(get(slValue(i,j,set=set), envir=ENVIR)))
    
    name <- Names[i]
    if (strictpos[[i]][j]) Evalue <- Evalue^2
    if (isinteger[[i]][j]) Evalue <- as.integer(base::round(Evalue))
    entryAssign(i, j, set=set, value=Evalue)

    assign("wait.simulation", get("wait.simulation", envir=ENVIR) - 1 , envir=ENVIR)
    if (get("wait.simulation", envir=ENVIR) < 0) {
      start_simu(update=name)
      tkRreplot(imgLU)
    }
    return(Evalue)
  }

  
  get.RS <- function() .Random.seed[length(.Random.seed)]
  
  start_simu <- function(update, all=FALSE) {
    Cat("start_simu ... ", if (missing(update)) "" else update, "...")
   
    repetitions <- Value("repetitions")
    time <- Value("time")
    burn.in <- Value("burn.in")
    totaltime <- time + burn.in
#    if (!exists("innovations", envir=ENVIR))
#      assign("innovations", matrix(NA, nrow=totaltime, ncol=repetitions))

    if (!exists("all.rs", envir=ENVIR)) {
      runif(1)
      all <- TRUE
    }
    
    if (ts.unif <- burn.in.unif <- do.distr <- all) {
      assign("all.rs",
             if (exists("all.rs", envir=ENVIR)) get("all.rs", envir=ENVIR)+100
	     else
	       base::round(runif(1,1,100000)), envir=ENVIR)
    } else {
      if (missing(update)) ERROR("update missing")
      switch(update,
	     "time" = ts.unif <- do.distr <- TRUE,
	     "burn.in" = burn.in.unif <- do.distr <- TRUE,
	     "repetitions" = ts.unif <- burn.in.unif <- do.distr <- TRUE,
	     "phi" = NULL,
	     "distr.param" = do.distr <- TRUE,
	     "distr.show" = do.distr <- TRUE,
	     "starting.ts" = NULL,
	     "starting.innovations" = do.distr <- TRUE,
	     {
	       ERROR("unknown parameter")
	     }
	     )
    }

    show <- Value("distr.show")
    if (ts.unif) {
      if (all) assign("ts.seed", get.RS(), envir=ENVIR)
      ts.unif <- matrix(nrow=time, ncol = repetitions)
      for (i in 1:repetitions) {
	set.seed(get("ts.seed", envir=ENVIR) + (i-1) * 100)
	ts.unif[, i] <- runif(time)
      }
      assign("ts.unif", ts.unif, envir=ENVIR)
    }

    if (burn.in.unif) {
      if (all) assign("burn.in.seed", get.RS(), envir=ENVIR)
      burn.in.unif <- matrix(nrow=burn.in, ncol = repetitions)
      for (i in 1:repetitions) {
	set.seed(get("burn.in.seed", envir=ENVIR) + (i-1) * 100)
	burn.in.unif[, i] <- rev(runif(burn.in))
      }
      assign("burn.in.unif", burn.in.unif, envir=ENVIR)
    }

    if (do.distr) {
      innovations <- matrix(NA, nrow=totaltime, ncol=repetitions * length(show))
      param <- Value("distr.param")
      unif <- rbind(get("burn.in.unif", envir=ENVIR),
		    get("ts.unif", envir=ENVIR))
      for (d in 1:length(L$distr))
	if (show[d]) {
	  idx <-(d-1) * repetitions + (1:repetitions)
#	  Print(idx, unif, param,  L$distr[[d]](unif, param))
	  innovations[, idx] <- L$distr[[d]](unif, param)
	}
      start.inno <- L$starting.matrix.innovations(repetitions, length(show),
						  Value("starting.innovations"))
      assign("innovations", rbind(start.inno, innovations), envir=ENVIR)
    }
    Cat("ende\n")
  }

  
  plotvoid <- function()
    plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")


  plotSimulation <- function() {
    ##    print("plotSimulation"); print(get("wait.simulation", envir=ENVIR))
    
     assign("wait.simulation",
	   get("wait.simulation", envir=ENVIR) - 1 , envir=ENVIR)
    if (get("wait.simulation", envir=ENVIR) >= 0) return()

    show <- which(Value("distr.show") != 0)
    ls <- length(show)
    if (!exists("all.rs", envir=ENVIR) || ls == 0) {
      plotvoid()
      return()
    }

    Cat("simu ... ")
    start.ts <-Value("starting.ts")
    start.inno <-Value("starting.innovations")
    phi <- Value("phi")
    repet <- Value("repetitions")
    col.idx <- (show - 1) * repet + (1:repet)
    inno <- innovations[, col.idx, drop=FALSE]
    burn.in <- Value("burn.in")
    ntime <- Value("time")
    burn.in <- Value("burn.in")
    row.idx.inno <- nrow(inno) + 1 - (ntime:1)
    start.ts.matrix <- L$starting.matrix.ts(repet, ls, start.ts)
    row.idx.ts <- function(M) nrow(M) + 1 - (ntime:1)
    stopifnot(all(row.idx.inno > 0))
    
    if (is(L$update.function, "RegisteredNativeSymbol")) {
      stop("not programmed yet")
    } else {
      X <- L$update.function(inno, phi, start.ts.matrix)
    }
    
    assign(simuValue(1), envir=ENVIR,
	   list(X=X, dim.i=dim(inno), dim.s=dim(start.ts) ))
    
    titles <- list()
    n.titles <- length(L$titles)
    for (i in 1:n.titles) titles[[i]] <- L$titles[[i]](phi)

    n.sets <- sum(sets)
    par(mfcol=c(currentNrVariab + 1, 1), cex=0.8, mar=c(3.2, 4.2, 0.2, 0.2))
   
    s.sets <- rev(which(sets))

    lty <- rep(show, each=repet)
    pch <- 1:repet
    if (simple.repet <- length(s.sets)==1 &&
	xor(repet > 1, ls > 1))
      rbcol <- c(colour, rainbow(repet * show))[1:(repet * ls)]
    type <- (if (sum(c(length(s.sets), repet, ls) != 1) <= 1) "l"
	     else "b")
    matplot(1:ntime, inno[row.idx.inno, ],
	    cex.lab=1.4, cex.axis=1.4,  xlab="", #"time",
	    ylab="innovations",
	    type=type,
	    col=if (simple.repet) rbcol else colour[1],
	    pch = pch, lty=lty)
    
    max.ylim <- 1e100
    for (m in 1:currentNrVariab) {
      M <- if (is.list(X)) X[[m]] else X	
      row.idx <- row.idx.ts(M)     
      y <-try(range(M[row.idx, ], na.rm=TRUE), silent=TRUE)
      if (is(y, "try-error") || any(is.na(y)) ||
	  y[2] > max.ylim || y[1] < -max.ylim) {
	plot(Inf, Inf, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
	text(0.5,0.5, label="no finite value with these innovations and these model parameters")
	next
      }
      
      plot(1:ntime, rep(NA, ntime), ylim=y,
	   cex.lab=1.4, cex.axis=1.4, xlab="", # "time",
	   ylab= if (n.titles > 0) titles[[1]][m] else ""
	   )
      for (set in s.sets) {
	simu <- get(simuValue(set), envir=ENVIR)
	M <- if (is.list(simu$X)) simu$X[[m]] else simu$X
	row.idx <- row.idx.ts(M)
	stopifnot(all(row.idx > 0))
	try(silent= TRUE,
	    matlines(1:ntime,
		     M[row.idx, ],
		     type=type,
		     col=if (simple.repet) rbcol else colour[set],
		     pch = pch, lty=lty))
    }
      if (n.titles > 1)
	for (i in 2:n.titles)
	  text(0, y[2] - diff(y) * (i-2), cex = 1.2,
	       adj=c(0,1), titles[[i]][m], col="darkred")
    }
   
    Cat("ende\n")
    
  }

  OnReturn <- function(...) {
    #RFoptions(LIST=get("RFopt.old", envir=ENVIR))
    assign(".tsgui.exit", GetL(), envir=parent.ev)
    tkDestroy(tt)    
  }

  
  position_j <- function(i, row, col) {
    if (Len[i] == 0) return(row)
    for (j in 1:Len[i]) {
      tkGridConf(get(labObj(i, j), envir=ENVIR), pady=1, 
		 column=col, row=row, sticky="w")
      if (i==globalParams + 1 && j==1) {
	tkGridConf(copyFctn, column=col + 0, row=row, sticky="e", pady=1)
      }
      if (islogical[[i]][j]) {
	tkGridConf(get(buttonWidget(i,j), envir=ENVIR), column=col + 0,
		   pady=1, row=row, sticky="e")
	if (sum(sets) > 1) {
	  for (set in which(sets)[-1]) {
	    tkGridConf(get(buttonWidget(i,j,set), envir=ENVIR),
		       pady=1, column=col + 2 * set + 1, row=row)
	  }
	  }	  
      } else {
	row <- row+1
	tkGridConf(get(slWidget(i,j), envir=ENVIR), column=col, row=row,
		   pady=1, sticky="w")
	tkGridConf(get(entryWidget(i,j), envir=ENVIR), column=col + 0,
		   pady=1, row=row, sticky="e")
	if (sum(sets) > 1) {
	  for (set in which(sets)[-1]) {
	    tkGridConf(get(entryWidget(i,j,set), envir=ENVIR),
		       pady=1, column=col + 2 * set + 1, row=row)
	  }
	}
      }
      row <- row+1
    }
    return(row)
  }
  
  position <- function(...) {
   # tkGridConf(emptyColumn, column=col.sl - 1, row=1, sticky="w")  
    ##--- Parameterwaehler ----------------------------------------------
    title.count <- 1

    row <- trd.row
    for (i in 1:globalParams) {      
      if (titles.guiloc[title.count] == i) {
	row <- trd.row
	tkGridConf(get(titles[title.count], envir=ENVIR),
		   column=trd.col[title.count + 1], row=row)
	row <- row + 1
	title.count <- title.count + 1
       }
      row <- position_j(i, row, trd.col[title.count])
    }
   
    row <- -1
    orig.title.count <- title.count
    for (i in (globalParams + 1):length(Names)) {
      if (titles.guiloc[title.count] == i) {
	row <- row + 2
	if (sum(sets) > 1 && title.count == orig.title.count) {
	  for (set in which(sets)[-1]) {
	    tkGridConf(get(deleteFctn(set), envir=ENVIR),
		       column=col.sl + 2 * set,
		       row=row, sticky="w")
	    tkGridConf(get(firstFctn(set), envir=ENVIR),
		       column=col.sl + 2 * set + 1,
		       row=row, sticky="w")
	  }
	} else {
	  tkGridConf(get(titles[title.count], envir=ENVIR),
		     column=col.sl, row=row)
	}
	row <- row + 1
	title.count <- title.count + 1
       }

      row <- position_j(i, row, col.sl)      
    }  
 
    ##--- Buttons - new simulation (new seed), return ---------------
    tkGridConf(buttonNewSimu,column=half.col, row=row<-trd.row+1, sticky="e")
    tkGridConf(buttonReturn, column=half.col, row=row <- row + 2, sticky="e")

    ##--- PLOT  ---------------------------------------------------------
    tkGridConf(imgLU, rowspan= 2 * image.rowspan, columnspan=2 * image.colspan,
	       column=fst.col, row=fst.row,sticky="w")

    if (length(modelclasses) > 1) {
      row <- trd.row
      tkGridConf(classTitle, column=snd.col, row=row, sticky="e")
      Namen <- names(modelclasses)      
      for (s in 1:length(modelclasses)) {
	tkGridConf(get(classButton(s), envir=ENVIR),
		   column=snd.col, row=row <- row + 1,
		   sticky="e")
      }
    }  


  } ## end fct position

    
  classButtons <- function() {
    Namen <- names(modelclasses)
    for (s in 1:length(modelclasses)) {
#      text <- paste("function(...) { delete.titles(); setClassLists(", s,
#		    ");classButtons();setClass()}")
      text <- paste("function(...) { ",
		    "GetL();",
		    "setClassLists(", s, "); setClass(); restart()}")
      Delete(X <- classButton(s))
      assign(X, envir=ENVIR,
	     tkButton(tt, text=Namen[s], command=eval(parse(text=text)),
		      fg=if (s == currentClass) col.alert else "black"
		      ))      
    }
  }

  Delete <- function(X)
    if (exists(X, envir=ENVIR, inherits=FALSE))
      do.call("rm", list(X, envir=ENVIR))
  
  SetSlider <- function(i, j, set=1, value, length.slider) {
    from <- minall[[i]][j]
    to <- maxall[[i]][j]
    s <- strictpos[[i]][j]    

    if (s) {
      value <- sqrt(value)
      from <- sqrt(from)
      to <- sqrt(to)
    }
    sltext <- paste("function(...) SliderChanges(i=", i, ", j=", j,
		    ", set=", set, ")")
    Delete(X <- slWidget(i, j, set))
    Delete(Y <- slValue(i, j, set))
    assign(Y, tkVar(value), envir=ENVIR)
    assign(X,
	   tkScale(tt, command = eval(parse(text=sltext)),
		   from= from, 
		   to = to,
		   showvalue=FALSE,
		   variable=get(Y),
		   ## neg value needed to get precise bounds
#		   bd = 1,
		   resolution=(if (isinteger[[i]][j]) -1
			       else (to - from) / numberSteps), 
		   orient="horizontal",
		   length=length.slider,
		   width=width.slider),
	   envir=ENVIR)
  }

  SetEntry <- function(i, j, set=1, value) {	  
    Delete(X <- entryValue(i, j, set))
    if (isinteger[[i]][j]) value <- as.integer(round(value))
    assign(X, tkVar(value), envir=ENVIR)
    text <- paste("function(...) EntryChanges(i=", i, ", j=", j,
		  ", set=", set, ")")
    Delete(X <- entryWidget(i, j, set))
    assign(X, tkEntry(tt, width=width.entry, borderwidth=1,
		      selectborderwidth=1,
		      textvariable=get(entryValue(i,j,set), envir=ENVIR)),
	   envir=ENVIR)
    tkBind(get(entryWidget(i, j, set), envir=ENVIR), "<Return>",
	   eval(parse(text=text)))
  }

  SetButton <- function(i, j, set=1, value) {
    Delete(X <- entryValue(i, j, set))
    assign(X, tkVar(value), envir=ENVIR)
    Delete(X <- buttonWidget(i, j, set))
    assign(X, envir=ENVIR,
	   tkCheckbutton(tt, onvalue = TRUE, offvalue=FALSE,
			 variable=get(entryValue(i,j,set), envir=ENVIR),
			 command=function(...) {
			   start_simu(all=TRUE)
			   tkRreplot(imgLU)}))
  }

  delete.titles <- function() {
    for (i in 1:length(Names)) {
      if (Len[i] > 0)
	for (j in 1:Len[i]) {
	  Delete(n <- labObj(i, j))
	  assign(n, tkLabel(tt, text = labName(i, j),fg=col_bg
			    ##			  ,bd=1
			    ),
		 envir=ENVIR)
	}
    }
    assign("wait.simulation", c(delete=0), envir=ENVIR)
    position()
  }


  starting.values <- function(sets=1:max.sets) {
    sub.n <- 1
    title.count <- 1
    for (i in 1:length(Names)) {
      ##      cat("starting.values ", i, "\n")
      if (titles.guiloc[title.count] == i) {
	importance <- titles.val[[title.count]]
	sub.n <- 1
	title.count <- title.count + 1
      }
      if (Len[i] > 0) {
	val <- L[[i]]
	val[!is.finite(val)] <- 1 ## dummy
	
	for (j in 1:Len[i]) {
	  for (set in sets) {
	    if (islogical[[i]][j]) {
	      SetButton(i, j, set, value=val[j])
	    } else {
	      SetSlider(i, j, set, value=val[j],
			if (title.count<2) length.slider.main
			else length.slider)
	      v <- as.numeric(tkValue(get(slValue(i,j,set),envir=ENVIR)))
	      if (strictpos[[i]][j]) v <- v^2
	      SetEntry(i, j, set, value=v)
	    }
	  }
	  Delete(n <- labObj(i, j))
	  ln <- labName(i, j)
	  assign(n, tkLabel(tt, text = ln,
			    fg=if (ln==noeffect) fg[length(fg)]
			    else fg[importance[sub.n]]
					#			  ,bd=1
			    ),
		 envir=ENVIR)
	}
      }
      sub.n <- sub.n + 1
    }

    if (length(do.not.show) > 0) {
      for (i in 1:length(do.not.show)) {
	if (Len[i] > 0) {
	  val <- L[[do.not.show[i]]]
	  for (j in 1:Len[i]) {
	    for (set in sets) {
	      Delete(X <- entryValue(do.not.show[i], j, set))
	      assign(X, tkVar(val[j]), envir=ENVIR)
	    }
	  }
	}
      }
    }
  }


  main <- function() {
    assign("sets", c(TRUE, rep(FALSE, max.sets - 1)), envir=ENVIR)
    setClassLists(currentClass)
    tt <- tcltk::tktoplevel()
    tcltk::tktitle(tt) <- "Time Series Gui"
    tcltk::tkwm.protocol(tt, "WM_DELETE_WINDOW", OnReturn)
    ##  tkGrid(tkLabel(tt, text="", width=1, height=0), column=0, row=0)
    tkGrid(tkLabel(tt, text="", width=1), column=fst.col+image.colspan,row=1)
    tkGrid(tkLabel(tt, text="BB", width=1), column=snd.col+image.colspan,row=1)
    tkGrid(tkLabel(tt, text="", width=1), column=col.sl+image.colspan, row=1)
    ##
    assign("col_bg", tkValue(tcltk::tkcget(tt, "-bg")), envir=ENVIR)
    ##  col_fg <- tkValue(tcltk::tkcget(tt, "-fg"))
    assign("tt", tt, envir = ENVIR)
    
    ## checkbuttion setzt die variable und fuehrt dann noch zusaetzlich
    ## command aus.
    
    ##--- Buttons ----------------------------
    assign("buttonNewSimu", envir=ENVIR,
	   tkButton(tt,text="New innovations",
		    command=function(...) {
		      start_simu(all=TRUE)
		      tkRreplot(imgLU)}))
    assign("buttonReturn", envir=ENVIR,
	   tkButton(tt, text="      Return       ",
		    command=OnReturn, fg = col.alert))
    assign("copyFctn", envir=ENVIR,
	   tkButton(tt, text="cpy", command=copy, fg="green"))
    assign("classTitle", envir=ENVIR, tkLabel(tt, text="MODELS", fg=col.alert))
    assign("emptyColumn", envir=ENVIR, tkLabel(tt, text="  "))
    for (i in 1:length(Names))
      assign(titles[i], tkLabel(tt, text=titles[i], fg=col.alert, bd=1),
	   envir=ENVIR)
    
    assign(simuValue(1), NULL, envir=ENVIR)
    for (set in 2:max.sets) {
      assign(simuValue(set), NULL, envir=ENVIR)
      text <- paste("function(...) delete(", set, ")")
      assign(deleteFctn(set),
	     tkButton(tt, text="d", command=eval(parse(text=text)),
		      bg=colour[set], fg=fgcol[set]), envir=ENVIR)
      text <-  paste("function(...) first(", set, ")")
      assign(firstFctn(set),
	     tkButton(tt, text="1", command=eval(parse(text=text)),
		      bg=colour[set], fg=fgcol[set]), envir=ENVIR)
    }
    
    assign("wait.simulation", c(start = 1 + sumLen - sum(unlist(islogical))),
				envir=ENVIR)
    classButtons()
    assign("imgLU", envir = ENVIR,
	   tkPlot(tt, 
		  fun = plotSimulation,
		  hscale=2 * plothscale, vscale=2 * plotvscale))
    starting.values()
    start_simu(all=TRUE)
        
    position()
    tkRreplot(imgLU)
  }
  
  restart <- function() {
    tkDestroy(tt)
    main()
  }
  
  main()
}
