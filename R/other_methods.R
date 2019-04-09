#' @title Pattern of Mortality Decline Prediction
#' @description Predict age-specific mortality rates using the Pattern of mortality decline (PMD) method (Andreev et al. 2013).
#' @details These functions implements the PMD method introduced in Andreev et al. (2013). 
#'     It assumes that the future decline in age-specific mortality will follow a certain pattern 
#'     with the increase in life expectancy at birth (e0): 
#'     \deqn{\log mx(t) = \log mx(t-1) - k(t) \rho_x(t)}
#'     
#'     Here, \eqn{\rho_x(t)} is the age-specific pattern of mortality decline between \eqn{t-1}
#'     and \eqn{t}. Such patterns for each sex and various levels of e0 
#'     are stored in the dataset \code{\link{PMDrho}}. The \code{pmd} function can be instructed 
#'     to interpolate between neighboring levels of e0 by setting the argument \code{interp.rho} 
#'     to \code{TRUE}. The \eqn{k} parameter is estimated to match the e0 level using the bisection 
#'     method.
#'     
#'     Function \code{pmd} evaluates the method for a single sex, while  \code{copmd} does it
#'     coherently for both sexes. In the latter case, the same \eqn{\rho_x} 
#'     (namely the average over sex-specific \eqn{\rho_x}) is used 
#'     for both, male and female.
#' @param e0 A vector of target life expectancy, one element for each predicted time point. 
#' @param mx0 A vector with starting age-specific mortality rates.
#' @param sex Either "male" or "female".
#' @param interp.rho Logical controlling if the \eqn{\rho} coefficients should be interpolated 
#'     (\code{TRUE}) or if the raw (binned) version should be used (\code{FALSE}), as stored in 
#'     the dataset \code{\link{PMDrho}}.
#' @param kranges A vector of size two, giving the min and max of the \eqn{k} parameter which is 
#'     estimated to match the target \code{e0} using the bisection method.
#' @param keep.lt Logical. If \code{TRUE} additional life table columns are kept in the 
#'     resulting object.
#' @param keep.rho Logical. If \code{TRUE} the \eqn{\rho} coefficients are included in the resulting object.
#' @return Function \code{pmd} returns a list with the following elements: a matrix \code{mx}
#'     with the predicted mortality rates. If \code{keep.lt} is \code{TRUE}, it also 
#'     contains matrices \code{sr} (survival rates), and life table quantities \code{Lx} and \code{lx}.
#'     If \code{keep.rho} is \code{TRUE}, it contains a matrix \code{rho} where columns correpond 
#'     to the values in the \code{e0} vector and rows correspond to age groups.
#' @export
#' 
#' @seealso \code{\link{mortcast}}, \code{\link{mortcast.blend}}, \code{\link{PMDrho}}
#' @references
#' Andreev, K. Gu, D., Gerland, P. (2013). Age Patterns of Mortality Improvement by Level of Life Expectancy at Birth with Applications to Mortality Projections. Paper presented at the Annual Meeting
#' of the Population Association of America, New Orleans, LA. \url{http://paa2013.princeton.edu/papers/132554}.
#' 
#' Gu, D., Pelletier, F. and Sawyer, C. (2017). Projecting Age-sex-specific Mortality: A Comparison of the Modified Lee-Carter and Pattern of Mortality Decline Methods, UN Population Division, 
#' Technical Paper No. 6. New York: United Nations. \url{https://population.un.org/wpp/Publications/Files/WPP2017_TechnicalPaperNo6.pdf}
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
#' pred <- pmd(e0f, mxf, sex = "female")
#' # plot first projection in black and the remaining ones in grey 
#' plot(pred$mx[,1], type="l", log="y", ylim=range(pred$mx),
#'     ylab="female mx", xlab="Age", main=country)
#' for(i in 2:ncol(pred$mx)) lines(pred$mx[,i], col="grey")
#'
#' @rdname pmdgroup

pmd <- function(e0, mx0, sex = c("male", "female"), interp.rho = FALSE,
                kranges = c(0.01, 25), keep.lt = FALSE, keep.rho = FALSE) {
    sex <- match.arg(sex)
    if(length(dim(e0)) > 0) e0 <- drop(as.matrix(e0)) # if it's a data.frame, it would not drop dimension without as.matrix
    if(length(dim(mx0)) > 0) mx0 <- drop(as.matrix(mx0)) 
    
    npred <- length(e0)
    nage <- length(mx0)
    # initialize results
    env <- new.env()
    data("PMDrho", envir = env)
    rho <- .find.pmd.rho(if(sex == "male") env$RhoMales else env$RhoFemales, 
                         e0, nage, npred, interp.rho = interp.rho)
    mx0l <- list(mx0)
    e0l <- list(e0)
    names(mx0l) <- names(e0l) <- sex
    res <- .do.copmd(e0l, mx0l, rho = rho, npred = npred, kranges = kranges,
                     keep.lt = keep.lt)
    if(keep.rho)
        res[[sex]]$rho <- rho
    return(res[[sex]])
}

.find.pmd.rho <- function(rho, e0, nage, npred, interp.rho = FALSE) {
    rhocols <- colnames(rho)
    rho.mids <- as.numeric(rhocols) # mid points
    brks <- c(0, rho.mids + 2.5, 200)
    # find rho for all e0
    this.rho <- matrix(0, nrow = nage, ncol = npred)
    for(time in 1:npred) {
        irho <- min(findInterval(e0[time], brks, left.open = FALSE), length(rhocols))
        rho.level <- rhocols[irho]
        this.rho[,time] <- unlist(rho[,rho.level])
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
    return(this.rho)
}

.do.copmd <- function(e0l, mx0l, rho, npred, kranges = c(0.01, 25), keep.lt = FALSE, 
                      sexratio.adjust = FALSE, adjust.with.mxf = FALSE) {
    # e0l and mx0l should be named lists of e0 and mx0 arrays with names being male and/or female. 
    # PMD is performed on all elements of the list using the same rho
    sexes <- c("female", "male")
    if(! ((all(names(mx0l) %in% sexes)) | all(names(e0l) %in% sexes)))
        stop("Names of the e0l and mx0l lists must be 'male' and/or 'female'.")
    nage <- length(mx0l[[1]])
    default.ages <- c(0, 1, seq(5, length = nage - 2, by = 5))
    # initialize results
    zeromatsr <- matrix(0, nrow=nage-1, ncol=npred)
    zeromatmx <- matrix(0, nrow=nage, ncol=npred)
    ressex <- list(mx=zeromatmx, lx=zeromatmx, sr=zeromatsr, Lx=zeromatsr)
    result <- list(female = ressex, male = ressex)
    constraint <- -1
    nconstr <- 0
    # iterate over sexes - rho stays the same
    for(sex in sexes) { # important that female is processed first because of a possible sex constraint
        if(!sex %in% names(mx0l)) next
        PMDres <- .C("PMD", as.integer(npred), as.integer(c(female=2, male=1)[sex]), 
                     as.integer(nage),
                 as.numeric(mx0l[[sex]]), as.numeric(rho), as.numeric(e0l[[sex]]), 
                 Kl=as.numeric(kranges[1]), Ku=as.numeric(kranges[2]), 
                 Constr = constraint, Nconstr = as.integer(nconstr),
                 LLm = as.numeric(result[[sex]]$Lx), Sr=as.numeric(result[[sex]]$sr), 
                 lx=as.numeric(result[[sex]]$lx), Mx=as.numeric(result[[sex]]$mx))
        ages <- names(mx0l[[sex]])
        if(is.null(ages)) ages <- default.ages
        result[[sex]]$mx <- matrix(PMDres$Mx, nrow=nage, 
                                   dimnames=list(ages, names(e0l[[sex]])))
        if(keep.lt) {
            result[[sex]]$sr <- matrix(PMDres$Sr, nrow=nage-1, dimnames=list(ages[-2], names(e0l[[sex]])))
            result[[sex]]$Lx <- matrix(PMDres$LLm, nrow=nage-1, dimnames=list(ages[-2], names(e0l[[sex]])))
            result[[sex]]$lx <- matrix(PMDres$lx, nrow=nage, dimnames=list(ages, names(e0l[[sex]])))
        } else {
            result[[sex]]$sr <- NULL
            result[[sex]]$Lx <- NULL
            result[[sex]]$lx <- NULL
        }
        if(sex == "female" && sexratio.adjust && "male" %in% names(mx0l)) { # both sexes must be present if applying constraint
            # compute minimum male mx
            if(adjust.with.mxf) { # using female mx
                minmx <- result$female$mx
            } else { # using Danan's regression
                env <- new.env()
                data("PMDadjcoef", envir = env)
                minmx <- matrix(-1, nrow = nrow(env$PMDadjcoef), ncol = npred)
                for(iage in 1:nrow(minmx)) {
                    coef <- env$PMDadjcoef[iage, ]
                    minmx[iage,] <-  10^(coef[,"intercept"] + coef[,"lmxf"]*log10(result[[sex]]$mx[iage,]) + 
                                         coef[,"e0f"]*e0l$female + coef[,"e0f2"]*e0l$female^2 + coef[,"gap"]*(e0l$female - e0l$male))
                }
            }
            minmx[,e0l$male > e0l$female] <- -1 # apply only if e0F >= e0M
            constraint <- as.numeric(minmx)
            nconstr <- nrow(minmx)
        }
    }
    return(result)
}

#' @export
#' @rdname pmdgroup
#' @param e0m A time series of target male life expectancy.
#' @param e0f A time series of target female life expectancy.
#' @param mxm0 A vector with starting age-specific male mortality rates.
#' @param mxf0 A vector with starting age-specific female mortality rates.
#' @param \dots Additional arguments passed to the underlying function. For \code{copmd}, in addition to
#'      \code{kranges} and \code{keep.lt}, it can be \code{sexratio.adjust} which is 
#'      a logical controlling if a sex-ratio adjustment should be applied to prevent crossovers 
#'      between male and female mx. In such a case it uses coefficients from the \code{\link{PMDadjcoef}} dataset. 
#'      However, if the argument \code{adjust.with.mxf} is set to \code{TRUE} (in addition to \code{sexratio.adjust}),
#'      the adjustment is done using the 
#'      female mortality rates as the lower constraint for male mortality rates.
#' @return Function \code{copmd} returns a list with one element for each sex 
#'     (\code{male} and \code{female}) where each of them is a list as described above.
#'     In addition if \code{keep.rho} is \code{TRUE}, element \code{rho.sex} 
#'     gives the sex-dependent (i.e. not averaged) \eqn{\rho_x} coefficient.
#' 
copmd <- function(e0m, e0f, mxm0, mxf0, interp.rho = FALSE, keep.rho = FALSE, ...) {
    e0  <- list(female=e0f, male=e0m)
    mx0  <- list(female=mxf0, male=mxm0)
    for(sex in names(e0)) {
        # convert to vectors if needed
        if(length(dim(e0[[sex]])) > 0) 
            e0[[sex]] <- drop(as.matrix(e0[[sex]])) # if it's a data.frame, it would not drop dimension without as.matrix
        if(length(dim(mx0[[sex]])) > 0) 
            mx0[[sex]] <- drop(as.matrix(mx0[[sex]])) 
    }
    npred <- length(e0$male)
    if(length(e0$female) != npred)
        stop("Mismatch in length of the e0 vectors.")
    nage <- length(mx0$male)
    if(length(mx0$female) != nage)
        stop("Mismatch in length of the mx0 vectors.")
    # derive rho as an average over male and female
    env <- new.env()
    data("PMDrho", envir = env)
    if(nage != nrow(env$RhoMales)) {
        warning("Mismatch in length of mx0 and the coefficient dataset. mx truncated to ",
                nrow(env$RhoMales), " age categories.")
        nage <- nrow(env$RhoMales)
        for(sex in names(e0)) mx0[[sex]] <- mx0[[sex]][1:nage]
    }
    rho.male <- .find.pmd.rho(env$RhoMales, e0$male, nage, npred, interp.rho = interp.rho)
    rho.female <- .find.pmd.rho(env$RhoFemales, e0$female, nage, npred, interp.rho = interp.rho)
    rho <- (rho.male + rho.female)/2
    res <- .do.copmd(e0, mx0, rho = rho, npred = npred, ...)

    if(keep.rho) {
        res$male$rho.sex <- rho.male
        res$female$rho.sex <- rho.female
        res$male$rho <- res$female$rho <- rho
    }
    return(res)
}

#' @title Model Life Tables Mortality Patterns
#' @description Predict age-specific mortality rates using Coale-Demeny and UN model life tables.
#' @details Given a level of life expectancy (e0), sex and a type of model life table, the function 
#'     extracts the corresponding mortality pattern from \code{\link{MLTlookup}}, 
#'     while interpolating between neighboring e0 groups.
#'     Function \code{mlt} is for one sex, while \code{mltj} can be used for both sexes.
#' @param e0 A time series of target life expectancy.
#' @param sex Either "male" or "female".
#' @param type Type of the model life table. Available options are \dQuote{CD_East}, \dQuote{CD_North}, 
#' \dQuote{CD_South}, \dQuote{CD_West}, \dQuote{UN_Chilean}, \dQuote{UN_Far_Eastern}, 
#' \dQuote{UN_General}, \dQuote{UN_Latin_American}, \dQuote{UN_South_Asian}. 
#' @return A matrix with the predicted mortality rates. Columns correpond 
#'     to the values in the \code{e0} vector and rows correspond to age groups.
#' @export
#' 
#' @seealso \code{\link{mortcast}}, \code{\link{mortcast.blend}}, \code{\link{pmd}}, \code{\link{MLTlookup}}
#' 
#' @references 
#' \url{https://population.un.org/wpp/Download/Other/MLT}
#' 
#' Coale, A., P. Demeny, and B. Vaughn. 1983. Regional model life tables and stable 
#' populations. 2nd ed. New York: Academic Press.
#' 
#' @examples
#' data(e0Fproj, package = "wpp2017")
#' country <- "Uganda"
#' # get target e0
#' e0f <- as.numeric(subset(e0Fproj, name == country)[-(1:2)])
#' # project into future using life table Cole-Demeny North
#' mx <- mlt(e0f, sex = "female", type = "CD_North")
#' # plot first projection in black and the remaining ones in grey 
#' plot(mx[,1], type="l", log="y", ylim=range(mx),
#'     ylab="female mx", xlab="Age", main=country)
#' for(i in 2:ncol(mx)) lines(mx[,i], col="grey")
#' 
#' @rdname mltgroup

mlt <- function(e0, sex = c("male", "female"), type = "CD_West") {
    sex <- match.arg(sex)
    sexcode <- c(female=2, male=1)[sex]
    if(length(dim(e0)) > 0) e0 <- drop(as.matrix(e0)) # if it's a data.frame, it would not drop dimension without as.matrix
    env <- new.env()
    data("MLTlookup", envir = env)
    conv.type <- gsub(" |_", "", type)
    conv.mlttypes <- gsub(" |_", "", env$MLTlookup$type)
    if(! conv.type %in% conv.mlttypes) {
        stop("Wrong MLT type. Available types:\n", paste(unique(env$MLTlookup$type), collapse = ", "))
    }
    mlt <- env$MLTlookup[conv.mlttypes == conv.type & 
                             env$MLTlookup$sex == sexcode, c("age", "e0", "mx")]

    mltw <- reshape(mlt, direction = "wide", timevar = "e0", idvar = "age", sep = "_")
    mlt.mat <- mltw[, -1]
    colnames(mlt.mat) <- sapply(strsplit(colnames(mlt.mat), "_"), function(x) x[2])
    rownames(mlt.mat) <- mltw[, 1]
    npred <- length(e0)
    nage <- nrow(mlt.mat)

    # initialize results
    zeromatsr <- matrix(0, nrow=nage-1, ncol=npred)
    zeromatmx <- matrix(0, nrow=nage, ncol=npred)
    res <- list(mx=zeromatmx, lx=zeromatmx, sr=zeromatsr, Lx=zeromatsr)

    mltcols <- colnames(mlt.mat)
    mlt.mids <- as.numeric(mltcols) # mid points
    brks <- c(0, mlt.mids, 200)
    this.mx <- matrix(0, nrow = nage, ncol = npred, 
                      dimnames = list(rownames(mlt.mat), names(e0)))
    for(time in 1:npred) {
        ie01 <- findInterval(e0[time], brks, left.open = FALSE) - 1
        ie02 <- ie01 + 1
        if(ie01 < 1) this.mx[, time] <- mlt.mat[, mltcols[1]] # first category
        else {
            if(ie02 > length(mltcols)) 
                this.mx[, time] <- mlt.mat[, mltcols[length(mltcols)]] # last category
            else {
                e0grid <- seq(mlt.mids[ie01], mlt.mids[ie02], length = 50)
                iorde0 <- which.min(abs(e0[time]-e0grid))
                this.mx[,time] <- apply(mlt.mat[,mltcols[c(ie01, ie02)]], 1, 
                                function(x) return(seq(x[1],x[2], length = 50)[iorde0]))
            }
        }
    }
    return(this.mx)
}

#' @export
#' @rdname mltgroup
#' @param e0m A time series of target male life expectancy.
#' @param e0f A time series of target female life expectancy.
#' @param \dots Additional arguments passed to the underlying function. 
#' 
mltj <- function(e0m, e0f, ...) {
    e0  <- list(female=e0f, male=e0m)
    res <- list()
    for(sex in names(e0)) {
        res[[sex]] <- list()
        if(!is.null(e0[[sex]]))
            res[[sex]]$mx <- mlt(e0 = e0[[sex]], sex = sex, ...)
    }
    return(res)
}

.apply.kannisto.if.needed <- function(mx, min.age.groups, ...) {
    if(nrow(mx$female$mx) <  min.age.groups || nrow(mx$male$mx) <  min.age.groups) {
    kan <- do.call("cokannisto", c(list(mx$male$mx, mx$female$mx), ...))
    mx$male$mx <- kan$male
    mx$female$mx <- kan$female
    }
    return(mx)
}

#' @title Mortality Prediction by Method Blending
#' @description Predict age-specific mortality rates using a blend of two different methods (Coherent Lee-Carter, 
#'     Pattern Mortality Decline, or Model Life Tables). Weights can be applied to fine-tune the blending mix.
#' @details The function allows to combine two different methods using given weights.
#'     The weights can change over time - by default they are interpolated from the starting weight 
#'     to the end weight. The projection is done for both sexes, so that coherent methods can be applied.
#' @param e0m A time series of future male life expectancy.
#' @param e0f A time series of future female life expectancy.
#' @param meth1 Character string giving the name of the first method to blend. It is one of 
#'     \dQuote{lc}, \dQuote{pmd} or \dQuote{mlt}, corresponding to Coherent Lee-Carter (function \code{\link{mortcast}}), 
#'      Pattern Mortality Decline (function \code{\link{copmd}}) and 
#'      Model Life Tables (function \code{\link{mltj}}), respectively.
#' @param meth2 Character string giving the name of the second method to blend. 
#'     One of the same choices as \code{meth1}.
#' @param weights Numeric vector with values between 0 and 1 giving the weight of \code{meth1}.
#'     If it is a single value, the same weight is applied for all time periods. 
#'     If it is a vector of size two, it is assumed these are weights for the first and the last
#'     time period. Remaining weights will be interpolated. Note that \code{meth2} is weighted 
#'     by \code{1 - weights}.
#' @param apply.kannisto Logical. If \code{TRUE} and if any of the methods results in less than 
#'     \code{min.age.groups} age categories, the coherent Kannisto method (\code{\link{cokannisto}}) 
#'     is applied to extend the age groups into old ages.
#' @param min.age.groups Minimum number of age groups. Triggers the application of Kannisto, see above. 
#' @param meth1.args List of arguments passed to the function that corresponds to \code{meth1}. 
#' @param meth2.args List of arguments passed to the function that corresponds to \code{meth2}. 
#' @param kannisto.args List of arguments passed to the \code{\link{cokannisto}} function if Kannisto is applied. 
#' @return List with elements \code{female} and \code{male}, each of which contains a matrix \code{mx}
#'     with the predicted mortality rates. In addition, it contains elements \code{meth1res} and \code{meth2res}
#'     which contain the results of the functions corresponding to the two methods. 
#'     Elements \code{meth1} and \code{meth2} contain the names of the methods. 
#'     A vector \code{weights} contains the final (possibly interpolated) weights.
#' @export
#' 
#' @seealso \code{\link{mortcast}}, \code{\link{copmd}}, \code{\link{mltj}}, 
#'     \code{\link{cokannisto}}
#'     
#' @examples
#' data(mxM, mxF, e0Fproj, e0Mproj, package = "wpp2017")
#' country <- "Brazil"
#' # estimate parameters from historical mortality data
#' mxm <- subset(mxM, name == country)[,4:16]
#' mxf <- subset(mxF, name == country)[,4:16]
#' rownames(mxm) <- rownames(mxf) <- c(0,1, seq(5, 100, by=5))
#' lcest <- lileecarter.estimate(mxm, mxf)
#' # project into future
#' e0f <- as.numeric(subset(e0Fproj, name == country)[-(1:2)])
#' e0m <- as.numeric(subset(e0Mproj, name == country)[-(1:2)])
#' # Blend LC and MLT
#' pred1 <- mortcast.blend(e0m, e0f, meth1 = "lc", meth2 = "mlt",
#'     meth1.args = list(lc.pars = lcest),
#'     meth2.args = list(type = "CD_North"),
#'     weights = c(1,0.25))
#' # Blend PMD and MLT
#' pred2 <- mortcast.blend(e0m, e0f, meth1 = "pmd", meth2 = "mlt",
#'     meth1.args = list(mxm0 = mxm[, "2010-2015"],
#'                       mxf0 = mxf[, "2010-2015"]))
#' # plot projection by time
#' plotmx <- function(pred, iage, main) 
#'     with(pred, {
#'         # blended projections 
#'         plot(female$mx[iage,], type="l", 
#'             ylim=range(meth1res$female$mx[iage,], 
#'                        meth2res$female$mx[iage,]),
#'             ylab="female mx", xlab="Time", main=main, col = "red")
#'         lines(meth1res$female$mx[iage,], lty = 2)
#'         lines(meth2res$female$mx[iage,], lty = 3)
#'         legend("topright", legend=c("blend", meth1, meth2),
#'                lty = 1:3, col = c("red", "black", "black"), bty = "n")
#'     })
#' age.group <- 3 # 5-9 years old
#' par(mfrow=c(1,2))
#' plotmx(pred1, age.group, "LC-MLT (age 5-9)")
#' plotmx(pred2, age.group, "PMD-MLT (age 5-9)")
#' 
mortcast.blend <- function(e0m, e0f, 
                          meth1 = "lc", meth2 = "mlt", weights = c(1, 0.5),
                          apply.kannisto = TRUE, min.age.groups = 28,
                          meth1.args = NULL, meth2.args = NULL, kannisto.args = NULL) {

    methods.allowed <- list(lc = "mortcast", mlt = "mltj", pmd = "copmd")
    meth1 <- match.arg(meth1, choices = names(methods.allowed))
    meth2 <- match.arg(meth2, choices = names(methods.allowed))
    
    npred <- length(e0m)
    w <- weights
    if(is.null(w)) w <- 0.5
    if(length(w) == 1) w <- rep(w, npred)
    if(!(length(w) == 2 || length(w) == npred))
        stop("Weights should be either of length 1 (constant weight), 2 (interpolated etween start and end) or the same size as e0.")
    if(any(w > 1 | w < 0)) stop("Weights must be between 0 and 1.")
    
    mx1 <- do.call(methods.allowed[[meth1]], c(list(e0m, e0f), meth1.args))
    mx2 <- do.call(methods.allowed[[meth2]], c(list(e0m, e0f), meth2.args))
    
    if(apply.kannisto) {
        mx1 <- .apply.kannisto.if.needed(mx1, min.age.groups, kannisto.args)
        mx2 <- .apply.kannisto.if.needed(mx2, min.age.groups, kannisto.args)
    }
    
    if(length(w) == 2 && npred > 2)  # interpolate weights
        w <- seq(w[1], w[2], length = npred)

    wmat <- NULL
    res <- list()
    for(sex in names(mx1)) {
        if(! sex %in% names(mx2)) next
        if(is.null(wmat))
           wmat <- matrix(w, ncol = npred, nrow = nrow(mx1[[sex]]$mx), byrow = TRUE)
        res[[sex]] <- list(mx = wmat * mx1[[sex]]$mx + (1-wmat) * mx2[[sex]]$mx)
    }
    return(c(res, list(meth1res = mx1, meth2res = mx2, 
                       meth1 = meth1, meth2 = meth2, weights = w)))
}