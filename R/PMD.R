#' @title Pattern of Mortality Decline Prediction
#' @description Predict age-specific mortality rates using the Pattern of mortality decline (PMD) method.
#' @details This function implements the PMD method published in Andreev et al. (2013). 
#'     The method assumes that the future decline in age-specific mortality will follow a certain pattern 
#'     with the increase in life expectancy at birth. 
#' @param mx0 A vector with starting age-specific mortality rates.
#' @param e0 A time series of target life expectancy.
#' @param sex Either "male" or "female".
#' @param interp.rho Logical controlling if the \eqn{\rho} coefficients should be interpolated 
#'     (\code{TRUE}) or binned (\code{FALSE}).
#' @param kranges A vector of size two, giving the min and max of the \eqn{k} parameter which is 
#'     estimated to match the the target \code{e0}. 
#' @param keep.lt Logical. If \code{TRUE} additional life table columns are kept in the 
#'     resulting object.
#' @param keep.rho Logical. If \code{TRUE} the \eqn{\rho} coefficients are included in the resulting object.
#' @return List with elements a matrix \code{mx}
#'     with the predicted mortality rates. If \code{keep.lt} is \code{TRUE}, it also 
#'     contains matrices \code{sr} (survival rates), and life table quantities \code{Lx} and \code{lx}.
#'     If \code{keep.rho} is \code{TRUE}, it contains a matrix \code{rho} where columns correpond 
#'     to the values in the \code{e0} vector and rows correspond to age groups.
#' @export
#' 
#' @seealso \code{\link{mortcast}}
#' @references
#' Andreev, K. Gu, D., Gerland, P. (2013). Age Patterns of Mortality Improvement by Level of Life Expectancy at Birth with Applications to Mortality Projections. Paper presented at the Annual Meeting
#' of the Population Association of America, New Orleans, LA. \url{http://paa2013.princeton.edu/papers/132554}.
#' 
#' 
#' @examples
#' data(mxF, e0Fproj, package = "wpp2017")
#' country <- "Hungary"
#' # get initial mortality for the current year
#' mxf <- subset(mxF, name == country)[,"2010-2015"]
#' names(mxf) <- c(0,1, seq(5, 100, by=5))
#' # get target e0
#' e0f <- as.numeric(subset(e0Fproj, name == country)[-(1:2)])
#' # project into future
#' pred <- pmd(mxf, e0f, sex = "female")
#' # plot first projection in black and the remaining ones in grey 
#' plot(pred$mx[,1], type="l", log="y", ylim=range(pred$mx),
#'     ylab="female mx", xlab="Age", main=country)
#' for(i in 2:ncol(pred$mx)) lines(pred$mx[,i], col="grey")
#' 

pmd <- function(mx0, e0, sex = c("male", "female"), interp.rho = FALSE,
                kranges = c(0.01, 25), keep.lt = FALSE, keep.rho = FALSE) {
    sex <- match.arg(sex)
    if(length(dim(e0)) > 0) e0 <- drop(as.matrix(e0)) # if it's a data.frame, it would not drop dimension without as.matrix
    if(length(dim(mx0)) > 0) mx0 <- drop(as.matrix(mx0)) 
    
    npred <- length(e0)
    nage <- length(mx0)
    # initialize results
    zeromatsr <- matrix(0, nrow=nage-1, ncol=npred)
    zeromatmx <- matrix(0, nrow=nage, ncol=npred)
    res <- list(mx=zeromatmx, lx=zeromatmx, sr=zeromatsr, Lx=zeromatsr)
    env <- new.env()
    data("rhoPMD", envir = env)
    rho <- if(sex == "male") env$RhoMales else env$RhoFemales
    rhocols <- colnames(rho)
    rho.mids <- as.numeric(rhocols) # mid points
    brks <- c(0, rho.mids + 2.5, 200)
    # find rho for all e0
    this.rho <- matrix(0, nrow = nage, ncol = npred)
    for(time in 1:npred) {
        irho <- min(findInterval(e0[time], brks, left.open = FALSE), length(rhocols))
        rho.level <- rhocols[irho]
        this.rho[,time] <- rho[,rho.level]
        if(interp.rho) { # interpolate coefficients
            if(e0[time] <= rho.mids[irho]) {
                irho1 <- irho-1
                irho2 <- irho
            } else {
                irho1 <- irho
                irho2 <- irho+1
            }
            if(irho1 > 0 && irho2 <= length(rho.mids)) {
                e0grid <- seq(rho.mids[irho1], rho.mids[irho2], length = 50)
                iorde0 <- which.min(abs(e0[time]-e0grid))
                this.rho[,time] <- apply(rho[,rhocols[c(irho1, irho2)]], 1, 
                    function(x) return(seq(x[1],x[2], length = 50)[iorde0]))
            }
        }
    }
    PMDres <- .C("PMD", as.integer(npred), as.integer(c(female=2, male=1)[sex]), as.integer(nage),
                as.numeric(mx0), as.numeric(this.rho), as.numeric(e0), 
                Kl=as.numeric(kranges[1]), Ku=as.numeric(kranges[2]), 
                LLm = as.numeric(res$Lx), Sr=as.numeric(res$sr), 
                lx=as.numeric(res$lx), Mx=as.numeric(res$mx))
    ages <- names(mx0)
    if(is.null(ages)) ages <- c(0, 1, seq(5, length = nage - 2, by = 5))
    res$mx <- matrix(PMDres$Mx, nrow=nage, dimnames=list(ages, names(e0)))
    if(keep.lt) {
        res$sr <- matrix(PMDres$Sr, nrow=nage-1, dimnames=list(ages[-2], names(e0)))
        res$Lx <- matrix(PMDres$LLm, nrow=nage-1, dimnames=list(ages[-2], names(e0)))
        res$lx <- matrix(PMDres$lx, nrow=nage, dimnames=list(ages, names(e0)))
    } else {
        res$sr <- NULL
        res$Lx <- NULL
        res$lx <- NULL
    }
    if(keep.rho)
        res$rho <- this.rho
    return(res)
}