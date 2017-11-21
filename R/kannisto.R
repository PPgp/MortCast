#' @title Kannisto Method
#' @description Extrapolate given mortality rates using the original Kannisto method.
#' @details 
#'     The function first estimates the original Kannisto parameters 
#'     by passing mortality rates for age groups \code{est.ages} into 
#'     the \code{\link{kannisto.estimate}} function.
#'     The estimated parameters are then passed to the projection function
#'     \code{\link{kannisto.predict}} to extrapolate into ages \code{proj.ages}.
#'     Lastly, the input mortality object is extended by results for the extrapolated ages. 
#'     If \code{proj.ages} contains age groups that are included in \code{mx},
#'     values for those age groups are overwritten.
#' @param mx A vector or matrix of mortality rates. If it is a matrix,
#'    rows correspond to age groups with rownames identifying the corresponding age as integers.
#'    By default five-years age groups are assigned to rows if rownames are not given.
#' @param est.ages A vector of integers identifying the ages to be used 
#'    for estimation. It should be a subset of rownames of \code{mx}.
#' @param proj.ages A vector of integers identifying the age groups for which mortality rates 
#'    are to be projected.
#' @return A vector or matrix containing the input mortality object \code{mx}
#'    extended by the extrapolated age groups.
#' @export
#' @references
#' Thatcher, A. R., Kannisto, V. and Vaupel, J. W. (1998). The Force of Mortality at Ages 80 to 120, 
#' volume 5 of Odense Monographs on Population Aging Series. Odense, Denmark: Odense University Press.
#' 
#' @seealso \code{\link{kannisto.estimate}}, \code{\link{kannisto.predict}}, \code{\link{cokannisto}}
#' @examples 
#' data(mxM, package = "wpp2017")
#' mx <- subset(mxM, name == "Burkina Faso")[,-(1:3)]
#' rownames(mx) <- c(0,1, seq(5, 100, by=5))
#' mxnew <- kannisto(mx)
#' ages <- as.integer(rownames(mxnew))
#' plot(ages, mxnew[,"2095-2100"], type="l", log="y", 
#'     xlab="age", ylab="mx", col="red")
#' lines(ages, mxnew[,"2010-2015"])
#' 
kannisto <- function(mx, est.ages = seq(80, 95, by=5), 
                     proj.ages = seq(100, 130, by=5)) {
    ages <- if(length(dim(mx)) == 0) names(mx) else rownames(mx) 
    Mxe <- as.matrix(mx)
    if(is.null(ages)) # default ages
        ages <- as.character(c(0,1, seq(5, by=5, length=nrow(Mxe)-2)))
    ages.num <- as.integer(ages)
    rownames(Mxe) <- ages
    est.ages.char <- as.character(est.ages)
    if(any(!est.ages.char %in% rownames(Mxe)))
        stop("est.ages are not included in mx.")
    est.data <- Mxe[est.ages.char,]
    kann.pars <- apply(est.data, 2, kannisto.estimate, ages = est.ages)
    res <- lapply(kann.pars, kannisto.predict, ages=proj.ages)
    mres <- sapply(res, cbind)
    proj.ages.char <- as.character(proj.ages)
    notincl <- which(!proj.ages.char %in% rownames(Mxe))
    resMx <- rbind(Mxe, matrix(NA, nrow=length(notincl), ncol=ncol(Mxe)))
    rownames(resMx)[(nrow(Mxe)+1):nrow(resMx)] <- proj.ages[notincl]
    resMx[proj.ages.char,] <- mres
    if(length(dim(mx)) == 0) { # convert to vector
        resMxv <- resMx[,1]
        names(resMxv) <- rownames(resMx)
        resMx <- resMxv
    }
    return(resMx)
}

#' @title Coherent Kannisto Method
#' @description Extrapolate given mortality rates into higher ages 
#'     using the Coherent Kannisto method as described in Sevcikova et al. (2016).
#' @details 
#'     The function first estimates the coherent Kannisto parameters 
#'     by passing mortality rates for age groups \code{est.ages} into 
#'     the \code{\link{cokannisto.estimate}} function.
#'     The estimated parameters are then passed to the projection function
#'     \code{\link{kannisto.predict}} to extrapolate into ages \code{proj.ages}.
#'     Lastly, the input mortality objects are extended by results for the extrapolated ages. 
#'     If \code{proj.ages} contains age groups that are included in \code{mxM} and \code{mxF},
#'     values for those age groups are overwritten.
#' @param mxM A vector or matrix of male mortality rates. If it is a matrix,
#'    rows correspond to age groups with rownames identifying the corresponding age as integers.
#'    By default five-years age groups are assigned to rows if rownames are not given.
#' @param mxF A vector or matrix of female mortality rates. Its length or dimension 
#'    should be the same \code{mxM}.
#' @param est.ages A vector of integers identifying the ages to be used 
#'    for estimation. It should be a subset of rownames of \code{mxM}.
#' @param proj.ages A vector of integers identifying the age groups for which mortality rates 
#'    are to be projected.
#' @return A list of two vectors or matrices (for male and female) containing the input motality 
#'    objects extended by the extrapolated age groups.
#' @export
#' 
#' @seealso \code{\link{cokannisto.estimate}}, \code{\link{kannisto.predict}}
#' @references
#' Sevcikova H., Li N., Kantorova V., Gerland P., Raftery A.E. (2016). 
#' Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. 
#' In: Schoen R. (eds) Dynamic Demographic Analysis. The Springer Series on Demographic Methods
#' and Population Analysis, vol 39. Springer, Cham
#' 
#' @examples 
#' data(mxM, mxF, package = "wpp2017")
#' country <- "South Africa"
#' mxm <- subset(mxM, name == country)[,-(1:3)]
#' mxf <- subset(mxF, name == country)[,-(1:3)]
#' rownames(mxm) <- rownames(mxf) <- c(0,1, seq(5, 100, by=5))
#' mxnew <- cokannisto(mxm, mxf)
#' ages <- as.integer(rownames(mxnew$male))
#' plot(ages, mxnew$male[,"2095-2100"], type="l", log="y", 
#'     xlab="age", ylab="mx", col="blue", main=country)
#' lines(ages, mxnew$female[,"2095-2100"], col="red")
#' lines(ages, mxnew$male[,"2010-2015"], lty=2, col="blue")
#' lines(ages, mxnew$female[,"2010-2015"], lty=2, col="red")
#' legend("bottomright", legend=c("male 2010-2015", "female 2010-2015",
#'     "male 2095-2100", "female 2095-2100"), bty="n",
#'     col=rep(c("blue", "red"),2), lty=c(2,2,1,1))
#' 
cokannisto <- function(mxM, mxF, 
                     est.ages = seq(80, 95, by=5), 
                     proj.ages = seq(100, 130, by=5)) {
    ages <- if(length(dim(mxM)) == 0) names(mxM) else rownames(mxM) 
    MxeM <- as.matrix(mxM)
    MxeF <- as.matrix(mxF)
    if(is.null(ages)) # default ages
        ages <- as.character(c(0,1, seq(5, by=5, length=nrow(MxeM)-2)))
    ages.num <- as.integer(ages)
    rownames(MxeM) <- rownames(MxeF) <- ages
    est.ages.char <- as.character(est.ages)
    if(any(!est.ages.char %in% rownames(MxeM)))
        stop("est.ages are not included in mxM.")
    est.data <- rbind(MxeM[est.ages.char,], MxeF[est.ages.char,])
    nest <- length(est.ages)
    kann.pars <- apply(est.data, 2, 
                       function(x) cokannisto.estimate(x[1:nest], x[(nest+1):length(x)], 
                                                              ages = est.ages))
    male.pars <- lapply(kann.pars, function(x) x[["male"]])
    resM <- lapply(male.pars, kannisto.predict, ages=proj.ages)
    matresM <- sapply(resM, cbind)
    female.pars <- lapply(kann.pars, function(x) x[["female"]])
    resF <- lapply(female.pars, kannisto.predict, ages=proj.ages)
    matresF <- sapply(resF, cbind)
    proj.ages.char <- as.character(proj.ages)
    notincl <- which(!proj.ages.char %in% rownames(MxeM))
    resMxM <- rbind(MxeM, matrix(NA, nrow=length(notincl), ncol=ncol(MxeM)))
    resMxF <- rbind(MxeF, matrix(NA, nrow=length(notincl), ncol=ncol(MxeF)))
    rownames(resMxM)[(nrow(MxeM)+1):nrow(resMxM)] <- rownames(resMxF)[(nrow(MxeF)+1):nrow(resMxF)] <- proj.ages[notincl]
    resMxM[proj.ages.char,] <- matresM
    resMxF[proj.ages.char,] <- matresF
    if(length(dim(mxM)) == 0) { # convert to vector
        resMxMv <- resMxM[,1]
        resMxFv <- resMxF[,1]
        names(resMxMv) <- names(resMxFv) <- rownames(resMxM)
        resMxM <- resMxMv
        resMxF <- resMxFv
    }
    return(list(male=resMxM, female=resMxF))
}

#' @title Kannisto Estimation
#' @description Estimate the Kannisto parameters (Thatcher et al. 1998).
#' @details Given the Kannisto equation \eqn{mlogit(m_x) = \log(c) + dx}{{mlogit(mx) = log(c) + dx}},
#'     the function estimates the \eqn{c} and \eqn{d} parameters using 
#'     values of \code{ages} as the covariate \eqn{x}.
#' @param mx A vector of mortality rates.
#' @param ages A vector of ages corresponding to \code{mx}.
#' @return List with Kannisto coefficients \eqn{c} and \eqn{d}.
#' @export
#' @references
#' Thatcher, A. R., Kannisto, V. and Vaupel, J. W. (1998). The Force of Mortality at Ages 80 to 120, 
#' volume 5 of Odense Monographs on Population Aging Series. Odense, Denmark: Odense University Press.
#' 
#' @seealso \code{\link{kannisto.predict}}, \code{\link{kannisto}}, \code{\link{cokannisto.estimate}}
#' 
#' @examples 
#' data(mxM, package = "wpp2017")
#' mx <- subset(mxM, name == "Canada")[,"2010-2015"]
#' kannisto.estimate(mx[18:21], ages = 18:21)
#' 
kannisto.estimate <- function(mx, ages){
    y <- log(mx) - log(1-mx)
    x <- ages
    coefs <- coefficients(lm(y ~ x))
    return(list(c=exp(coefs[[1]]), d=coefs[['x']]))
}

#' @title Coherent Kannisto Estimation
#' @description Estimate the coherent Kannisto parameters as described in Sevcikova et al. (2016).
#' @details Given the Kannisto equation \eqn{mlogit(m_x) = \log(c) + dx}{{mlogit(mx) = log(c) + dx}},
#'     the Coherent Kannisto method estimates the \eqn{d} parameter jointly for male and female 
#'     data, in order to prevent mortality crossovers in higher ages. 
#' @param mxM A vector of male mortality rates.
#' @param mxF A vector of female mortality rates.
#' @param ages A vector of ages corresponding to \code{mxM} and \code{mxF}.
#' @return List of two lists, each with coherent Kannisto coefficients \eqn{c} and \eqn{d}, 
#'   one for male and one for female. The \eqn{d} values are the same in both lists.
#' @export
#' 
#' @seealso \code{\link{cokannisto}}, \code{\link{kannisto.predict}}, \code{\link{kannisto}}
#' 
#' @references
#' Sevcikova H., Li N., Kantorova V., Gerland P., Raftery A.E. (2016). 
#' Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. 
#' In: Schoen R. (eds) Dynamic Demographic Analysis. The Springer Series on Demographic Methods
#' and Population Analysis, vol 39. Springer, Cham
#' 
#' @examples 
#' data(mxM, mxF, package = "wpp2017")
#' country <- "Brazil"
#' mxm <- subset(mxM, name == country)[,"2010-2015"]
#' mxf <- subset(mxF, name == country)[,"2010-2015"]
#' cokannisto.estimate(mxm[18:21], mxf[18:21], ages = 18:21)
#' 
cokannisto.estimate <- function(mxM, mxF, ages){
    mx <- c(mxM, mxF)
    y <- log(mx) - log(1-mx)
    x <- c(ages, ages)
    g <- c(rep(1, length(mxM)), rep(0, length(mxF)))
    coefs <- coefficients(lm(y ~ g + x))
    return(list(female=list(c=exp(coefs[[1]]), d=coefs[['x']]),
                male=list(c=exp(coefs[[1]] + coefs[['g']]),
                            d=coefs[['x']]))
                )
}

#' @title Kannisto Prediction 
#' @description Given estimated Kannisto parameters (coherent or original), 
#'     it predicts mortality rates for given ages.
#' @details Given parameters \eqn{c} and \eqn{d} in \code{pars}, 
#'     the function uses the Kannisto equation \eqn{mlogit(m_x) = \log(c) + dx}{{mlogit(mx) = log(c) + dx}},
#'     to predict mortality rates for age groups \eqn{x} given by \code{ages}.
#' @param pars A list with Kanisto coefficients \eqn{c} and \eqn{d} 
#'     (e.g. result of \code{\link{kannisto.estimate}} or \code{\link{cokannisto.estimate}}).
#' @param ages A vector of ages to make prediction for.
#' @return Vector of predicted mortality rates.
#' @export
#' @seealso \code{\link{cokannisto}}, \code{\link{kannisto.estimate}}, \code{\link{cokannisto.estimate}}
#' @references
#' Thatcher, A. R., Kannisto, V. and Vaupel, J. W. (1998). The Force of Mortality at Ages 80 to 120, 
#' volume 5 of Odense Monographs on Population Aging Series. Odense, Denmark: Odense University Press.
#'
#' @examples 
#' data(mxM, mxF, package = "wpp2017")
#' mxm <- subset(mxM, name == "Germany")[,"2010-2015"]
#' ages <- c(0, 1, seq(5, 130, by=5))
#' 
#' # using original Kannisto parameters
#' pars <- kannisto.estimate(mxm[18:21], ages = ages[18:21])
#' mxm.pred <- kannisto.predict(pars, ages = ages[22:28])
#' plot(ages, c(mxm[1:21], mxm.pred), type="l", log="y", 
#'     xlab="age", ylab="mx")
#'     
#' # Coherent Kannisto
#' mxf <- subset(mxF, name == "Germany")[,"2010-2015"]
#' copars <- cokannisto.estimate(
#'    mxm[18:21], mxf[18:21], ages = ages[18:21])
#' cmxm.pred <- kannisto.predict(copars[["male"]], ages = ages[22:28])
#' cmxf.pred <- kannisto.predict(copars[["female"]], ages = ages[22:28])
#' plot(ages, c(mxm[1:21], cmxm.pred), type="l", log="y", 
#'     xlab="age", ylab="mx", col="blue")
#' lines(ages, c(mxf[1:21], cmxf.pred), col="red")
#' 
kannisto.predict <- function(pars, ages){
    numer <- pars$c * exp(pars$d * ages)
    return(numer/(1 + numer))
}

