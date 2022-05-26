#' @title Lee-Carter Estimation
#' @description Estimate Lee-Carter parameters (Lee and Carter 1992).
#' @details 
#' The function estimates parameters of \eqn{\log(m_x(t)) = a_x + b_x k(t) + \epsilon_x(t)} (Lee and Carter 1992). 
#' The argument \code{ax.index} determines which time periods to use to 
#' estimate the \eqn{a_x} parameter, while \code{ax.smooth} controls if 
#' the resulting \eqn{a_x} should be smoothened over ages (see Sevcikova et al. 2016 for details).
#' 
#' @param mx A matrix of age-specific mortality rates where rows correspond to age groups
#'     and columns correspond to time periods. Rownames define the starting ages of the age groups.
#' @param ax.index A vector of column indices of \code{mx} to be used to estimate the \eqn{a_x} parameter.
#'     By default all time periods are used.
#' @param ax.smooth Logical allowing to smooth the \eqn{a_x} over ages.
#' @param ax.smooth.df Degree of freedom for smoothing if \code{ax.smooth} is \code{TRUE}. 
#' Default is half the length of \eqn{a_x}.
#' @param bx.postprocess Logical determining if numerical anomalies in \eqn{b_x} should be dealt with.
#' @param nx Size of age groups. By default ages are determined by rownames of \code{mx}. This argument is only used if 
#'     \code{mx} has no rownames. If \code{nx} is 5, the age groups are interpreted as 0, 1, 5, 10, \dots. For \code{nx} equals 1, 
#'     the age groups are interpreted as 0, 1, 2, 3, \dots.
#' @return List with elements \code{ax}, \code{bx} and \code{kt} corresponding to the estimated parameters.
#' @export
#' @seealso \code{\link{mortcast}}, \code{\link{lileecarter.estimate}}
#' @references 
#' Lee, R. D. and Carter, L. (1992). Modeling and forecasting the time series of 
#' US mortality. Journal of the American Statistical Association, 87, 659-671.
#' 
#' Sevcikova H., Li N., Kantorova V., Gerland P., Raftery A.E. (2016). 
#' Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. 
#' In: Schoen R. (eds) Dynamic Demographic Analysis. The Springer Series on Demographic Methods
#' and Population Analysis, vol 39. Springer, Cham
#' 
#' @examples
#' data(mxM, package = "wpp2017")
#' mx <- subset(mxM, name == "Netherlands")[,4:16]
#' rownames(mx) <- c(0,1, seq(5, 100, by=5))
#' lc.ax.avg <- leecarter.estimate(mx)
#' lc.ax.last <- leecarter.estimate(mx, ax.index=ncol(mx))
#' plot(lc.ax.avg$ax, type="l")
#' lines(lc.ax.last$ax, col="blue")
#' 
leecarter.estimate <- function(mx, ax.index = NULL, ax.smooth = FALSE, 
                               ax.smooth.df = NULL, bx.postprocess = TRUE, nx = 5) {
    if(length(dim(mx)) < 2 || ncol(mx) < 2) stop("mx must be a matrix (age x time) of at least two columns.")
    lmx <- log(mx)
    if(length(dim(lmx))==0)
        lmx <- as.matrix(lmx)
    if(is.null(rownames(lmx))) # default ages
        rownames(lmx) <- if(nx > 1) c(0,1, seq(nx, by=nx, length=nrow(lmx)-2)) else 0:(nrow(lmx)-1)
    nest <- ncol(lmx)
    if(is.null(ax.index)) ax.index <- 1:nest
    ax <- apply(lmx[,ax.index, drop=FALSE], 1, sum, na.rm=TRUE) / length(ax.index)
    if(ax.smooth) {
        if (is.null(ax.smooth.df)) ax.smooth.df <- ceiling(length(ax)/2)
        ax.sm <- smooth.spline(ax, df = ax.smooth.df)$y
        ax[-1] <- ax.sm[-1] # keep value of the first age group
    }
    kt <- apply(lmx, 2, sum) - sum(ax)
    bx <- bx.estimate(lmx, ax, kt, postprocess = bx.postprocess)
    return(list(ax = ax, bx = bx, kt = kt, nx = nx))
}

bx.estimate <- function(lmx, ax, kt, postprocess=TRUE) {
    x2 <- sum(kt*kt)
    x1 <- rep(NA, nrow(lmx))
    for (i in 1:nrow(lmx)) 
        x1[i] <- sum((lmx[i,]-ax[i])*kt)
    bx <- x1/x2
    if(postprocess)
        bx <- .finish.bx(bx)
    names(bx) <- rownames(lmx)
    return(bx)
}

.finish.bx <- function(bx) {
    negbx <- which(bx <= 0)
    lnegbx <- length(negbx)
    if(lnegbx > 0 && negbx[1] == 1) {
        bx[1] <- 0
        negbx <- if(lnegbx > 1) negbx[2:lnegbx] else c()
        lnegbx <- length(negbx)
    }
    while(lnegbx > 0) {
        bx[negbx] <- 0.5 * bx[negbx-1]
        negbx <- which(bx < 0)
        lnegbx <- length(negbx)
    }
    lbx <- length(bx)
    for (i in 1:(lbx-1)) { 
        if (bx[lbx+1 - i] == 0) bx[lbx+1 - i] <- bx[lbx - i]
    }
    bx <- bx/sum(bx) # must sum to 1
    return(bx)
}

#' @title Coherent Lee-Carter Estimation
#' @description Estimate coherent Lee-Carter parameters (Li and Lee 2005).
#' @details The coherent Lee-Carter parameters for male and female mortality rates 
#'     share the same \eqn{b_x} which is the average of the age-specific 
#'     \eqn{b_x} parameters. 
#'     
#'     The function in addition computes the ultimate \eqn{b^u_x} as defined in 
#'     Li et al. (2013) based on the coherent \eqn{b_x}.
#' @param mxM A matrix of male age-specific mortality rates where rows correspond to age groups
#'     and columns correspond to time periods. For 5-year age groups, the first and second rows corresponds to 
#'     0-1 and 1-5 age groups, respectively. Rownames define the starting ages of the respective groups.
#' @param mxF A matrix of female mortality rates of the same shape as \code{mxM}.
#' @param nx Size of age groups. Should be either 5 or 1.
#' @param ... Additional arguments passed to \code{\link{leecarter.estimate}}.
#' @return List containing elements \code{bx} (coherent \eqn{b_x} parameter), 
#'     \code{ultimate.bx} (ultimate \eqn{b^u_x} parameter), \code{ages} (age groups), \code{nx} (age group interval), and
#'   lists \code{female} and \code{male}, each with the Lee-Carter parameters.
#' @export
#' 
#' @references 
#' Li, N. and Lee, R. D. (2005). Coherent mortality forecasts for a group of populations: 
#' An extension of the Lee-Carter method. Demography, 42, 575-594.
#' 
#' Li, N., Lee, R. D. and Gerland, P. (2013). Extending the Lee-Carter method to model the rotation 
#' of age patterns of mortality decline for long-term projections. Demography, 50, 2037-2051.
#' 
#' @examples
#' data(mxM, mxF, package = "wpp2017")
#' country <- "Germany"
#' mxm <- subset(mxM, name == country)[,4:16]
#' mxf <- subset(mxF, name == country)[,4:16]
#' rownames(mxm) <- rownames(mxf) <- c(0,1, seq(5, 100, by=5))
#' lc <- lileecarter.estimate(mxm, mxf)
#' plot(lc$bx, type="l")
#' lines(lc$ultimate.bx, lty=2)
#' 
lileecarter.estimate <- function(mxM, mxF, nx = 5, ...) {
    lc.male <- leecarter.estimate(mxM, nx = nx, ...)
    lc.female <- leecarter.estimate(mxF, nx = nx, ...)
    if(length(lc.female$ax) != length(lc.male$ax) || length(lc.female$kt) != length(lc.male$kt))
        stop("Mismatch in dimensions of male and female mortality. Check the mxM and mxF arguments.")
    bx <- (lc.male$bx + lc.female$bx)/2
    ages <- as.integer(names(lc.female$bx))
    # correct nx if wrong input
    nx <- ages[4] - ages[3]
    return(list(bx = bx, ultimate.bx = ultimate.bx(bx), ages = ages, nx = nx, 
                female=list(ax=lc.female$ax, bx=bx, kt=lc.female$kt, sex.bx = lc.female$bx),
                male=list(ax=lc.male$ax, bx=bx, kt=lc.male$kt, sex.bx = lc.male$bx)
            ))
}


#' @title Rotated Lee-Carter
#' @description Rotate the Lee-Carter parameter \eqn{b_x} over time to reach an ultimate \eqn{b^u_x}, 
#'     as described in Li et al. (2013).
#' 
#' @param bx A vector of the Lee-Carter \eqn{b_x} parameter
#'    (from e.g. \code{\link{lileecarter.estimate}} or \code{\link{leecarter.estimate}}).
#' @param ultimate.bx A vector of the ultimate \eqn{b^u_x} parameter as defined in Li, Lee, Gerland (2013)
#'    (obtained using \code{\link{lileecarter.estimate}} or \code{\link{ultimate.bx}}).    
#' @param e0 A time series of life expectancies.
#' @param e0l Level of life expectancy at which the rotation starts.
#' @param e0u Level of life expectancy at which the rotation finishes.
#' @param p Exponent of the smooth function.
#' @return Function \code{rotate.leecarter} returns a matrix of rotated \eqn{B_x(t)} where rows correspond to age groups and columns 
#'    correspond to time periods (given by the vector \code{e0}). 
#' @export
#' 
#' @references
#' Li, N., Lee, R. D. and Gerland, P. (2013). Extending the Lee-Carter method to model the rotation 
#' of age patterns of mortality decline for long-term projections. Demography, 50, 2037-2051.
#'
#' @examples
#' data(mxF, mxM, e0Fproj, e0Mproj, package = "wpp2017")
#' country <- "Japan"
#' mxm <- subset(mxM, name == country)[,4:16]
#' mxf <- subset(mxF, name == country)[,4:16]
#' e0f <- as.numeric(subset(e0Fproj, name == country)[-(1:2)])
#' e0m <- as.numeric(subset(e0Mproj, name == country)[-(1:2)])
#' rownames(mxm) <- rownames(mxf) <- c(0,1, seq(5, 100, by=5))
#' lc <- lileecarter.estimate(mxm, mxf)
#' rotlc <- rotate.leecarter(lc$bx, lc$ultimate.bx, (e0f + e0m)/2)
#' plot(lc$bx, type="l")
#' lines(lc$ultimate.bx, col="red")
#' for(i in 1:ncol(rotlc)) lines(rotlc[,i], col="grey")
#'   
rotate.leecarter <- function(bx, ultimate.bx, e0, e0l = 80, e0u = 102, p = 0.5) {
    npred <- length(e0)
    wt <- (e0 - e0l)/(e0u-e0l)
    wst <- (0.5*(1+(sin(pi/2*(2*wt-1)))))^p
    Bxt <- matrix(NA, nrow=length(bx), ncol=npred)
    for(t in 1:npred) {
        Bxt[,t] <- switch(cut(e0[t], c(0, e0l, e0u, 9999), labels=FALSE, right=FALSE),
                          bx, 
                          (1-wst[t])*bx + wst[t]*ultimate.bx,
                          ultimate.bx)
    }
    dimnames(Bxt) <- list(names(bx), names(e0))
    return(Bxt)
}

#' @rdname rotate.leecarter
#' @return Function \code{ultimate.bx} returns a vector of the ultimate \eqn{b^u_x}.
#' @export
ultimate.bx <- function(bx) {
    # ultimate bx (Li, Lee, Gerland 2013)
    bux <- bx
    lbx <- length(bx)
    idx15 <- which(names(bx) == 15)
    idx60 <- which(names(bx) == 60)
    bux[1:idx60] <- mean(bux[idx15:idx60])
    bux[(idx60+1):lbx] <- bux[(idx60+1):lbx] * (bux[idx60]/bux[idx60+1]) # adjust so that b(65-69)=b(60-64)
    return(bux/sum(bux)) # must sum to 1
}

.kranges <- function(Bxt, axm, axf) {
    lmin <- -12
    lmax <- 0
    machine.max <- log(.Machine$double.xmax)
    machine.min <- log(.Machine$double.xmin)
    npred <- ncol(Bxt)
    kranges <- list(male=list(), female=list())
    ax <- list(male=axm, female=axf)

    bx.lim <- rbind(apply(Bxt, 2, function(x) min(x[x>0])), 
                    apply(Bxt, 2, function(x) max(x[x>0])))
    for(sex in c("male", "female")) {
        this.lmin <- floor(min(lmin, min(ax[[sex]])-0.1))
        kranges[[sex]]$kl <- pmin((lmax - apply(ax[[sex]], 2, max))/bx.lim[2,], machine.max)
        kranges[[sex]]$ku <- pmax((this.lmin - apply(ax[[sex]], 2, min))/bx.lim[1,], machine.min)
    }
    return(kranges)
}

#' @title Coherent Rotated Lee-Carter Prediction
#' @description Predict age-specific mortality rates using the coherent rotated Lee-Carter method.
#' @details This function implements Steps 6-9 of Algorithm 2 in Sevcikova et al. (2016). 
#'     It uses the abridged or unabridged life table function to find the level of mortality that coresponds to the given 
#'     life expectancy. Thus, it can be used for both, mortality for 5- or 1-year age groups. 
#' @param e0m A time series of future male life expectancy.
#' @param e0f A time series of future female life expectancy.
#' @param lc.pars A list of coherent Lee-Carter parameters with elements \code{bx}, \code{ultimate.bx},
#'     \code{ages}, \code{nx}, 
#'     \code{female} and \code{male} as returned by \code{\link{lileecarter.estimate}}. 
#'     The \code{female} and \code{male} objects are again lists that should contain a vector
#'     \code{ax} and optionally a matrix \code{axt} if the \eqn{a_x} parameter 
#'     needs to be defined as time dependent. In such a case, rows are age groups and columns are 
#'     time periods corresponding to the length of the \code{e0f} and \code{e0m} vectors.
#' @param rotate If \code{TRUE} the rotation method of \eqn{b_x} is used as described in Li et al. (2013).
#' @param keep.lt Logical. If \code{TRUE} additional life table columns are kept in the 
#'     resulting object.
#' @param constrain.all.ages By default the method constrains the male mortality to be above female 
#'     mortality for old ages if the male life expectancy is below the female life expectancy. Setting 
#'     this argument to \code{TRUE} causes this constraint to be applied to all ages.
#' @param \dots Additional life table arguments.
#' @return List with elements \code{female} and \code{male}, each of which contains a matrix \code{mx}
#'     with the predicted mortality rates. If \code{keep.lt} is \code{TRUE}, it also 
#'     contains matrices \code{sr} (survival rates), and life table quantities \code{Lx} and \code{lx}.
#' @export
#' 
#' @seealso \code{\link{rotate.leecarter}}, \code{\link{leecarter.estimate}}, \code{\link{lileecarter.estimate}},
#'     \code{\link{mortcast.blend}}
#' @references
#' Li, N., Lee, R. D. and Gerland, P. (2013). Extending the Lee-Carter method to model the rotation 
#' of age patterns of mortality decline for long-term projections. Demography, 50, 2037-2051.
#' 
#' Sevcikova H., Li N., Kantorova V., Gerland P., Raftery A.E. (2016). 
#' Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. 
#' In: Schoen R. (eds) Dynamic Demographic Analysis. The Springer Series on Demographic Methods
#' and Population Analysis, vol 39. Springer, Cham
#' 
#' @examples
#' # estimate parameters from historical mortality data (5-year age groups)
#' data(mxM, mxF, e0Fproj, e0Mproj, package = "wpp2017")
#' country <- "Brazil"
#' mxm <- subset(mxM, name == country)[,4:16]
#' mxf <- subset(mxF, name == country)[,4:16]
#' rownames(mxm) <- rownames(mxf) <- c(0,1, seq(5, 100, by=5))
#' lc <- lileecarter.estimate(mxm, mxf)
#' 
#' # project into future for given levels of life expectancy
#' e0f <- as.numeric(subset(e0Fproj, name == country)[-(1:2)])
#' e0m <- as.numeric(subset(e0Mproj, name == country)[-(1:2)])
#' pred <- mortcast(e0m, e0f, lc)
#' 
#' # plot first projection in black and the remaining ones in grey 
#' plot(lc$ages, pred$female$mx[,1], type="b", log="y", ylim=range(pred$female$mx),
#'     ylab="female mx", xlab="Age", main=paste(country, "(5-year age groups)"), cex=0.5)
#' for(i in 2:ncol(pred$female$mx)) lines(lc$ages, pred$female$mx[,i], col="grey")
#'
#' # similarly for 1-year age groups
#' # derive toy 1-year mx using model life tables at given level of e0
#' mxm1y <- mlt(seq(65, 71, length = 4), sex = "male", nx = 1)
#' mxf1y <- mlt(seq(73, 78, length = 4), sex = "female", nx = 1)
#' 
#' # estimate parameters
#' lc1y <- lileecarter.estimate(mxm1y, mxf1y, nx = 1)
#'  
#' # project into the future 
#' pred1y <- mortcast(e0m, e0f, lc1y)
#' 
#' # plot first projection in black and the remaining ones in grey 
#' plot(lc1y$ages, pred1y$female$mx[,1], type="b", log="y", ylim=range(pred1y$female$mx),
#'     ylab="female mx", xlab="Age", main="1-year age groups", cex=0.5)
#' for(i in 2:ncol(pred1y$female$mx)) lines(lc1y$ages, pred1y$female$mx[,i], col="grey")
#' 

mortcast <- function (e0m, e0f, lc.pars, rotate = TRUE, keep.lt = FALSE, 
                      constrain.all.ages = FALSE, ...) {
    get.a0rule <- function(a0rule = c("ak", "cd"), ...)
        list(ak = 1, cd = 2)[[match.arg(a0rule)]]
    # if e0 is a data.frame, convert to vector (it would not drop dimension without as.matrix)
    if(length(dim(e0m)) > 0) e0m <- drop(as.matrix(e0m)) 
    if(length(dim(e0f)) > 0) e0f <- drop(as.matrix(e0f)) 
    
    # prepare for computation
    e0  <- list(female=e0f, male=e0m)
    npred <- length(e0f)
    nage <- length(lc.pars$ages)
    if(lc.pars$nx > 1) {
        resnage <-  nage-1 # number of age groups of the resulting matrices 
        age.groups <- lc.pars$ages[-2] # group 0-1 collapsed into 0-5
    } else {
        resnage <-  nage # all ages
        age.groups <- lc.pars$ages
    }
    a0cat <- get.a0rule(...)
    zeromatsr <- matrix(0, nrow=resnage, ncol=npred)
    zeromatmx <- matrix(0, nrow=nage, ncol=npred)
    # in an abridged case, lx, Lx and sr will be returned for 5-year intervals
    ressex <- list(mx=zeromatmx, lx=zeromatsr, sr=zeromatsr, Lx=zeromatsr)
    result <- list(female = ressex, male = ressex)
    # rotate bx if needed
    if(rotate)
        Bxt <- rotate.leecarter(lc.pars$bx, lc.pars$ultimate.bx, (e0f + e0m)/2)
    else
        Bxt <- matrix(lc.pars$bx, nrow=nage, ncol=npred)
    
    # allow for time-dependent ax
    for(sex in c("female", "male")) 
        if(is.null(lc.pars[[sex]]$axt)) 
            lc.pars[[sex]]$axt <- matrix(lc.pars[[sex]]$ax, nrow=nrow(Bxt), ncol=npred)
    
    # compute ranges for k(t)
    kranges <- .kranges(Bxt, lc.pars$male$axt, lc.pars$female$axt)

    #Get the projected kt from e0, and make projection of Mx
    for (sex in c("female", "male")) { # iterate over female and male (order matters because of the constrain)
        LCres <- .C("LC", as.integer(npred), as.integer(c(female=2, male=1)[sex]), as.integer(nage), as.integer(lc.pars$nx),
                  as.numeric(lc.pars[[sex]]$axt), as.numeric(Bxt), as.numeric(e0[[sex]]), 
                  Kl=as.numeric(kranges[[sex]]$kl), Ku=as.numeric(kranges[[sex]]$ku), 
                  # 1 for constraining old ages only; 2 for constraining all ages
                  constrain=as.integer((sex == "male") * ((sex == "male") + (constrain.all.ages == TRUE))), 
                  FMx=as.numeric(result$female$mx), FEop=as.numeric(e0$female), a0rule = as.integer(a0cat),
                  LLm = as.numeric(result[[sex]]$Lx), Sr=as.numeric(result[[sex]]$sr), 
                  lx=as.numeric(result[[sex]]$lx), Mx=as.numeric(result[[sex]]$mx))
        result[[sex]]$mx <- matrix(LCres$Mx, nrow=nage, 
                                   dimnames=list(lc.pars$ages, names(e0m)))
        result[[sex]]$nx <- lc.pars$nx
        if(keep.lt) {
            result[[sex]]$sr <- matrix(LCres$Sr, nrow=resnage,
                                       dimnames=list(age.groups, names(e0m)))
            result[[sex]]$Lx <- matrix(LCres$LLm, nrow=resnage,
                                       dimnames=list(age.groups, names(e0m)))
            result[[sex]]$lx <- matrix(LCres$lx, nrow=resnage, 
                                       dimnames=list(age.groups, names(e0m)))
        } else {
            result[[sex]]$sr <- NULL
            result[[sex]]$Lx <- NULL
            result[[sex]]$lx <- NULL
        }
    }
    return(result)    
}
