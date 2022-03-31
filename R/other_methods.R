#' @title Pattern of Mortality Decline Prediction
#' @description Predict age-specific mortality rates using the Pattern of mortality decline (PMD) method (Andreev et al. 2013).
#' @details These functions implements the PMD method introduced in Andreev et al. (2013) and its modifications. 
#'     It assumes that the future decline in age-specific mortality will follow a certain pattern 
#'     with the increase in life expectancy at birth (e0): 
#'     \deqn{log[mx(t)] = log[mx(t-1)] - k(t) \rho_x(t)}
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
#' @param mx0 A vector with starting age-specific mortality rates. In case of \code{modpmd} it can be 
#'     a matrix where rows correspond to age groups
#'     and columns correspond to time periods. Rownames define the starting ages of the age groups.
#' @param sex Either "male" or "female".
#' @param nx Size of age groups. Should be either 5 or 1.
#' @param interp.rho Logical controlling if the \eqn{\rho} coefficients should be interpolated 
#'     (\code{TRUE}) or if the raw (binned) version should be used (\code{FALSE}), as stored in 
#'     the dataset \code{\link{PMDrho}}.
#' @param kranges A vector of size two, giving the min and max of the \eqn{k} parameter which is 
#'     estimated to match the target \code{e0} using the bisection method.
#' @param keep.lt Logical. If \code{TRUE} additional life table columns are kept in the 
#'     resulting object.
#' @param keep.rho Logical. If \code{TRUE} the \eqn{\rho} coefficients are included in the resulting object.
#' @return Function \code{pmd} and \code{modpmd} return a list with the following elements: a matrix \code{mx}
#'     with the predicted mortality rates. If \code{keep.lt} is \code{TRUE}, it also 
#'     contains matrices \code{sr} (survival rates), and life table quantities \code{Lx} and \code{lx}.
#'     If \code{keep.rho} is \code{TRUE}, it contains a matrix \code{rho} where columns correpond 
#'     to the values in the \code{e0} vector and rows correspond to age groups.
#' @export
#' 
#' @seealso \code{\link{mortcast}}, \code{\link{mortcast.blend}}, \code{\link{PMDrho}}
#' @references
#' Andreev, K., Gu, D., Gerland, P. (2013). Age Patterns of Mortality Improvement by Level of Life Expectancy at Birth with Applications to Mortality Projections. Paper presented at the Annual Meeting
#' of the Population Association of America, New Orleans, LA. \url{https://paa2013.princeton.edu/papers/132554}.
#' 
#' Gu, D., Pelletier, F., Sawyer, C. (2017). Projecting Age-sex-specific Mortality: A Comparison of the Modified Lee-Carter and Pattern of Mortality Decline Methods, UN Population Division, 
#' Technical Paper No. 6. New York: United Nations. \url{https://population.un.org/wpp/Publications/Files/WPP2017_TechnicalPaperNo6.pdf}
#' 
#' @examples
#' data(mxF, e0Fproj, package = "wpp2017")
#' country <- "Hungary"
#' # get initial mortality for the current year
#' mxf <- subset(mxF, name == country)[,"2010-2015"]
#' names(mxf) <- c(0,1, seq(5, 100, by=5))
#' # get target e0
#' e0f <- subset(e0Fproj, name == country)[-(1:2)]
#' # project into future
#' pred <- pmd(e0f, mxf, sex = "female")
#' # plot first projection in black and the remaining ones in grey 
#' plot(pred$mx[,1], type = "l", log = "y", ylim = range(pred$mx),
#'     ylab = "female mx", xlab = "Age", main = country)
#' for(i in 2:ncol(pred$mx)) lines(pred$mx[,i], col = "grey")
#'
#' @rdname pmdgroup

pmd <- function(e0, mx0, sex = c("male", "female"), nx = 5, interp.rho = FALSE,
                kranges = c(0, 25), keep.lt = FALSE, keep.rho = FALSE, ...) {
    sex <- match.arg(sex)
    if(length(dim(e0)) > 0) e0 <- drop(as.matrix(e0)) # if it's a data.frame, it would not drop dimension without as.matrix
    if(length(dim(mx0)) > 0) mx0 <- drop(as.matrix(mx0)) 
    
    npred <- length(e0)
    nage <- length(mx0)
    # initialize results
    rho <- .find.pmd.rho(if(sex == "male") MortCast::RhoMales else MortCast::RhoFemales, 
                         e0, nage, npred, interp.rho = interp.rho, nx = nx)
    mx0l <- list(mx0)
    e0l <- list(e0)
    names(mx0l) <- names(e0l) <- sex
    res <- .do.copmd(e0l, mx0l, rho = rho, npred = npred, nx = nx, kranges = kranges,
                     keep.lt = keep.lt, ...)
    if(keep.rho)
        res[[sex]]$rho <- rho
    return(res[[sex]])
}

#' @export
#' @rdname pmdgroup
#' @details Function \code{modpmd} implements a modified version of \code{pmd} where the initial \eqn{log[mx(t_0)]}
#' is replaced by an \eqn{a_x} estimated as in \code{\link{leecarter.estimate}}, i.e. using possibly 
#' multiple years of historical \code{mx} and optionally smoothed. Arguments \code{ax.index}, \code{ax.smooth} and 
#' \code{ax.smooth.df} determine the estimation years and parameters of the smoothing. 
#' 
#' @param ax.index A vector of column indices of \code{mx} to be used to estimate the \eqn{a_x = E[log(mx(t_0))]} parameter.
#'     By default it is estimated as the average over all observed time periods, but this argument can restrict the time periods 
#'     to use.
#' @param ax.smooth Logical allowing to smooth the \eqn{a_x} over ages.
#' @param ax.smooth.df Degree of freedom for smoothing if \code{ax.smooth} is \code{TRUE}. 
#'     Default is half the length of \eqn{a_x}.
#' 


modpmd <- function(e0, mx0, sex = c("male", "female"), nx = 5, interp.rho = FALSE,
                kranges = c(0, 25), ax.index = NULL, ax.smooth = FALSE, 
                ax.smooth.df = NULL, keep.lt = FALSE, keep.rho = FALSE, ...) {
    sex <- match.arg(sex)
    mx0 <- .prepare.mx.for.modpmd(mx0, nx = nx)
    if(length(dim(e0)) > 0) e0 <- drop(as.matrix(e0)) # if it's a data.frame, it would not drop dimension without as.matrix
    npred <- length(e0)

    # compute ax
    mx0v <- .compute.mx0.for.modpmd(mx0, ax.index = ax.index, ax.smooth = ax.smooth, 
                                   ax.smooth.df = ax.smooth.df)

    nage <- length(mx0v) # mx0v is a vector
    # find rho and initialize results
    rho <- .find.pmd.rho(if(sex == "male") MortCast::RhoMales else MortCast::RhoFemales, 
                         e0, nage, npred, interp.rho = interp.rho, nx = nx)
    mx0l <- list(mx0v)
    e0l <- list(e0)
    names(mx0l) <- names(e0l) <- sex
    
    # compute pmd
    res <- .do.copmd(e0l, mx0l, rho = rho, npred = npred, nx = nx, kranges = kranges,
                     keep.lt = keep.lt, ...)
    if(keep.rho)
        res[[sex]]$rho <- rho
    return(res[[sex]])
}

.prepare.mx.for.modpmd <- function(mx, nx = 5){
    if(length(dim(mx))==0)
        mx <- as.matrix(mx)
    if(is.null(rownames(mx))) # default ages
        rownames(mx) <- if(nx > 1) c(0,1, seq(nx, by=nx, length=nrow(mx)-2)) else 0:(nrow(mx)-1)
    return(mx)
}

.compute.mx0.for.modpmd <- function(mx, ax.index = NULL, ax.smooth = FALSE, ax.smooth.df = NULL){
    # compute ax
    lmx <- log(mx)
    nest <- ncol(lmx)
    if(is.null(ax.index)) ax.index <- 1:nest
    ax <- apply(lmx[,ax.index, drop=FALSE], 1, sum, na.rm=TRUE) / length(ax.index)
    if(ax.smooth) {
        if (is.null(ax.smooth.df)) ax.smooth.df <- ceiling(length(ax)/2)
        ax.sm <- smooth.spline(ax, df = ax.smooth.df)$y
        ax[-1] <- ax.sm[-1] # keep value of the first age group
    }
    return(exp(ax))
}

.find.pmd.rho <- function(rho, e0, nage, npred, interp.rho = FALSE, nx = 5) {
    rhocols <- colnames(rho)
    rho.mids <- as.numeric(rhocols) # mid points
    brks <- c(0, rho.mids + 2.5, 200)
    # find rho for all e0
    this.rho <- matrix(0, nrow = nrow(rho), ncol = npred)
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
    # if this.rho has a different number of age groups than desired, interpolate between ages
    if(nx == 1) {
        interp <- function(x)
                      approx(as.integer(rownames(rho)), x, xout = seq(0, nage - 1), method = "linear")$y
        this.rho <- apply(this.rho, 2, interp)
    }
    return(this.rho)
}

.do.copmd <- function(e0l, mx0l, rho, npred, nx = 5, kranges = c(0, 25), keep.lt = FALSE, 
                      sexratio.adjust = FALSE, adjust.sr.if.needed = FALSE, 
                      adjust.with.mxf = FALSE, a0rule = c("ak", "cd")) {
    # e0l and mx0l should be named lists of e0 and mx0 arrays with names being male and/or female. 
    # PMD is performed on all elements of the list using the same rho
    sexes <- c("female", "male")
    if(! ((all(names(mx0l) %in% sexes)) | all(names(e0l) %in% sexes)))
        stop("Names of the e0l and mx0l lists must be 'male' and/or 'female'.")
    if(!nx %in% c(1, 5))
        stop("The nx argument must be either 5 or 1.")
    nage <- length(mx0l[[1]])
    if(nx == 5) {
        default.ages <- c(0, 1, seq(5, length = nage - 2, by = 5))
        resnage <-  nage-1 # number of age groups of the resulting matrices 
        age.groups <- default.ages[-2] # group 0-1 collapsed into 0-5
    } else {
        default.ages <- seq(0, nage - 1)
        resnage <-  nage # all ages
        age.groups <- default.ages
    }
    a0cat <- list(ak = 1, cd = 2)[[match.arg(a0rule)]]
    # initialize results
    zeromatsr <- matrix(0, nrow=resnage, ncol=npred)
    zeromatmx <- matrix(0, nrow=nage, ncol=npred)
    ressex <- list(mx=zeromatmx, lx=zeromatsr, sr=zeromatsr, Lx=zeromatsr)
    result <- list(female = ressex, male = ressex)
    sex.ratio.ini <- rep(0, nage)
    constraint <- -1
    nconstr <- 0
    # iterate over sexes - rho stays the same
    for(sex in sexes) { # important that female is processed first because of a possible sex constraint
        if(!sex %in% names(mx0l)) next
        #if(sex == "male") stop("")
        PMDres <- .C("PMD", as.integer(npred), as.integer(c(female=2, male=1)[sex]), 
                     as.integer(nage), as.integer(nx),
                 as.numeric(mx0l[[sex]]), as.numeric(rho), as.numeric(e0l[[sex]]), 
                 Kl=as.numeric(kranges[1]), Ku=as.numeric(kranges[2]), 
                 Constr = constraint, Nconstr = as.integer(nconstr), ConstrIfNeeded = as.integer(adjust.sr.if.needed == TRUE && sex == "male"),
                 FMx = as.numeric(result$female$mx), SRini = sex.ratio.ini, a0rule = as.integer(a0cat),
                 LLm = as.numeric(result[[sex]]$Lx), Sr=as.numeric(result[[sex]]$sr), 
                 lx=as.numeric(result[[sex]]$lx), Mx=as.numeric(result[[sex]]$mx), PACKAGE = "MortCast")
        ages <- names(mx0l[[sex]])
        if(is.null(ages)) ages <- default.ages
        result[[sex]]$mx <- matrix(PMDres$Mx, nrow=nage, 
                                   dimnames=list(ages, names(e0l[[sex]])))
        result[[sex]]$nx <- nx
        if(keep.lt) {
            result[[sex]]$sr <- matrix(PMDres$Sr, nrow=resnage, dimnames=list(age.groups, names(e0l[[sex]])))
            result[[sex]]$Lx <- matrix(PMDres$LLm, nrow=resnage, dimnames=list(age.groups, names(e0l[[sex]])))
            result[[sex]]$lx <- matrix(PMDres$lx, nrow=resnage, dimnames=list(age.groups, names(e0l[[sex]])))
        } else {
            result[[sex]]$sr <- NULL
            result[[sex]]$Lx <- NULL
            result[[sex]]$lx <- NULL
        }
        if(sex == "female" && (sexratio.adjust || adjust.sr.if.needed) && "male" %in% names(mx0l)) { # both sexes must be present if applying constraint
            # compute minimum male mx
            if(adjust.sr.if.needed) { # 
                sex.ratio.ini <- as.numeric(mx0l[["male"]]/mx0l[["female"]])
            } else {
                minmx <- NULL
                if(adjust.with.mxf) { # using female mx
                    minmx <- result$female$mx
                } else { # using Danan's regression
                    if(nx == 1) warning("No PMD regression adjustment for nx = 1")
                    else {
                        minmx <- matrix(-1, nrow = nrow(MortCast::PMDadjcoef), ncol = npred)
                        for(iage in 1:nrow(minmx)) {
                            coef <- MortCast::PMDadjcoef[iage, ]
                            minmx[iage,] <-  10^(coef[,"intercept"] + coef[,"lmxf"]*log10(result[[sex]]$mx[iage,]) + 
                                         coef[,"e0f"]*e0l$female + coef[,"e0f2"]*e0l$female^2 + coef[,"gap"]*(e0l$female - e0l$male))
                        }
                    }
                }
                if(!is.null(minmx)) {
                    minmx[,e0l$male > e0l$female] <- -1 # apply only if e0F >= e0M
                    constraint <- as.numeric(minmx)
                    nconstr <- nrow(minmx)
                }
            }
        }
    }
    return(result)
}

#' @export
#' @rdname pmdgroup
#' @param e0m A time series of target male life expectancy.
#' @param e0f A time series of target female life expectancy.
#' @param mxm0,mxf0 A vector with starting age-specific male/female mortality rates. If \code{use.modpmd} is \code{TRUE},
#'    this can be a matrix of historical mx (age x time) from which the starting values are estimated.
#' @param nx Size of age groups. Should be either 5 or 1.
#' @param use.modpmd Logical determining if the modified version of PMD (\code{modpmd}) should be used. 
#'    In such a case the starting values of mortality rates are estimated similarly to \eqn{a_x} in 
#'    \code{\link{leecarter.estimate}}, possibly from more than one time periods. In addition, a smoothing can be applied.
#' @param \dots Additional arguments passed to the underlying functions. For \code{copmd}, in addition to
#'      \code{kranges} and \code{keep.lt}, it can be \code{sexratio.adjust} which is 
#'      a logical controlling if a sex-ratio adjustment should be applied to prevent crossovers 
#'      between male and female mx. In such a case it uses coefficients from the \code{\link{PMDadjcoef}} dataset. 
#'      However, if the argument \code{adjust.with.mxf} is set to \code{TRUE} (in addition to \code{sexratio.adjust}),
#'      the adjustment is done using the 
#'      female mortality rates as the lower constraint for male mortality rates. 
#'      If the argument \code{adjust.sr.if.needed} is set to \code{TRUE}, a sex-ratio adjustment
#'      is performed dynamically, using the sex ratio in the previous time point. 
#'      In such a case, an adjustment in time t is applied only if there was a drop of sex ratio 
#'      below one at time t-1. Other arguments passed here in \code{copmd} can be \code{ax.index}, \code{ax.smooth} and
#'      \code{ax.smooth.df} which control the estimation of the initial mx if \code{use.modpmd} is \code{TRUE}.
#' @return Function \code{copmd} returns a list with one element for each sex 
#'     (\code{male} and \code{female}) where each of them is a list as described above.
#'     In addition if \code{keep.rho} is \code{TRUE}, element \code{rho.sex} 
#'     gives the sex-dependent (i.e. not averaged) \eqn{\rho_x} coefficient.
#' 
copmd <- function(e0m, e0f, mxm0, mxf0, nx = 5, interp.rho = FALSE, keep.rho = FALSE, 
                  use.modpmd = FALSE, ...) {
    e0  <- list(female=e0f, male=e0m)
    mx0  <- list(female=mxf0, male=mxm0)
    dotargs <- list(...)
    m0args <- list()
    if(use.modpmd)  # extract modpmd arguments
        m0args <- dotargs[names(dotargs) %in% names(as.list(args(.compute.mx0.for.modpmd)))]

    for(sex in names(e0)) {
        # convert to vectors or matrix as needed
        if(length(dim(e0[[sex]])) > 0) 
            e0[[sex]] <- drop(as.matrix(e0[[sex]])) # if it's a data.frame, it would not drop dimension without as.matrix
        if(use.modpmd){
            mx0[[sex]] <- .prepare.mx.for.modpmd(mx0[[sex]], nx = nx)
            mx0[[sex]] <- do.call(".compute.mx0.for.modpmd", c(list(mx0[[sex]]), m0args))
        } else {
            if(length(dim(mx0[[sex]])) > 0) 
                mx0[[sex]] <- drop(as.matrix(mx0[[sex]])) 
        }
    }
    npred <- length(e0$male)
    if(length(e0$female) != npred)
        stop("Mismatch in length of the e0 vectors.")
    nage <- length(mx0$male)
    if(length(mx0$female) != nage)
        stop("Mismatch in length of the mx0 objects.")

    # derive rho as an average over male and female
    if(nx == 5 && nage != nrow(MortCast::RhoMales)) {
        warning("Mismatch in length of mx0 and the coefficient dataset. mx truncated to ",
                nrow(MortCast::RhoMales), " age categories.")
        nage <- nrow(MortCast::RhoMales)
        for(sex in names(e0)) mx0[[sex]] <- mx0[[sex]][1:nage]
    }
    rho.male <- .find.pmd.rho(MortCast::RhoMales, e0$male, nage, npred, nx = nx, interp.rho = interp.rho)
    rho.female <- .find.pmd.rho(MortCast::RhoFemales, e0$female, nage, npred, nx = nx, interp.rho = interp.rho)
    rho <- (rho.male + rho.female)/2
    res <- do.call(".do.copmd", c(list(e0, mx0, rho = rho, npred = npred, nx = nx), 
                                  dotargs[!names(dotargs) %in% names(m0args)]))

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
#'     extracts the corresponding mortality pattern from \code{\link{MLTlookup}} (for abridged LT) 
#'     or \code{\link{MLT1Ylookup}} (for 1-year LT), 
#'     while interpolating between neighboring e0 groups.
#'     Function \code{mlt} is for one sex, while \code{mltj} can be used for both sexes.
#' @param e0 A time series of target life expectancy.
#' @param sex Either "male" or "female".
#' @param type Type of the model life table. Available options are \dQuote{CD_East}, \dQuote{CD_North}, 
#' \dQuote{CD_South}, \dQuote{CD_West}, \dQuote{UN_Chilean}, \dQuote{UN_Far_Eastern}, 
#' \dQuote{UN_General}, \dQuote{UN_Latin_American}, \dQuote{UN_South_Asian}. 
#' @param nx Size of age groups. Should be either 5 or 1.
#' @param \dots Not used.
#' @return Function \code{mlt} returns a matrix with the predicted mortality rates. Columns correspond 
#'     to the values in the \code{e0} vector and rows correspond to age groups. 
#'     Function \code{mltj} returns a list of such matrices, one for each sex.
#' @export
#' 
#' @seealso \code{\link{mortcast}}, \code{\link{mortcast.blend}}, \code{\link{pmd}}, \code{\link{MLTlookup}}
#' 
#' @references 
#' \url{https://www.un.org/development/desa/pd/data/extended-model-life-tables}
#' 
#' Coale, A., P. Demeny, and B. Vaughn. 1983. Regional model life tables and stable 
#' populations. 2nd ed. New York: Academic Press.
#' 
#' @examples
#' data(e0Fproj, package = "wpp2017")
#' country <- "Uganda"
#' # get target e0
#' e0f <- subset(e0Fproj, name == country)[-(1:2)]
#' # project into future using life table Cole-Demeny North
#' mx <- mlt(e0f, sex = "female", type = "CD_North")
#' # plot first projection in black and the remaining ones in grey 
#' plot(mx[,1], type = "l", log = "y", ylim = range(mx),
#'     ylab = "female mx", xlab = "Age", 
#'     main = paste(country, "5-year age groups"))
#' for(i in 2:ncol(mx)) lines(mx[,i], col = "grey")
#' 
#' # MLT for 1-year age groups
#' mx1y <- mlt(e0f, sex = "female", type = "CD_North", nx = 1)
#' plot(mx1y[,1], type = "l", log = "y", ylim = range(mx1y),
#'     ylab = "female mx", xlab = "Age", 
#'     main = paste(country, "1-year age groups"))
#' for(i in 2:ncol(mx1y)) lines(mx1y[,i], col = "grey")
#'     
#' @rdname mltgroup

mlt <- function(e0, sex = c("male", "female"), type = "CD_West", nx = 5, ...) {
    sex <- match.arg(sex)
    sexcode <- c(female=2, male=1)[sex]
    if(length(dim(e0)) > 0) e0 <- drop(as.matrix(e0)) # if it's a data.frame, it would not drop dimension without as.matrix
    if(!nx %in% c(1, 5))
        stop("The nx argument must be either 5 or 1.")
    lookup <- if(nx == 5) MortCast::MLTlookup else MortCast::MLT1Ylookup
    conv.type <- gsub(" |_", "", type)
    conv.mlttypes <- gsub(" |_", "", lookup$type)
    if(! conv.type %in% conv.mlttypes) {
        stop("Wrong MLT type. Available types:\n", paste(unique(lookup$type), collapse = ", "))
    }
    mlt <- lookup[conv.mlttypes == conv.type & 
                             lookup$sex == sexcode, c("age", "e0", "mx")]

    mltw <- reshape(mlt, direction = "wide", timevar = "e0", idvar = "age", sep = "_")
    mlt.mat <- mltw[, -1]
    colnames(mlt.mat) <- sapply(strsplit(colnames(mlt.mat), "_"), function(x) x[2])
    rownames(mlt.mat) <- mltw[, 1]
    npred <- length(e0)
    nage <- nrow(mlt.mat)

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
mltj <- function(e0m, e0f, ..., nx = 5) {
    e0  <- list(female=e0f, male=e0m)
    res <- list()
    for(sex in names(e0)) {
        res[[sex]] <- list()
        if(!is.null(e0[[sex]])) {
            res[[sex]]$mx <- mlt(e0 = e0[[sex]], sex = sex, ..., nx = nx)
            res[[sex]]$nx <- nx
        }
    }
    return(res)
}

#' @title Log-Quadratic Mortality Model
#' @description Predict age-specific mortality rates using the Log-Quadratic Mortality Model (Wilmoth et al. 2012).
#' @details The LogQuad method in this implementation projects mortality rates using the equation
#'      \deqn{\log(m_x) = a_x + b_x h + c_x h^2 + v_x k}
#'      where \eqn{a_x}, \eqn{b_x}, \eqn{c_x} and \eqn{v_x} are age-specific coefficients, \eqn{h = \log( 5q0 )} 
#'      (i.e. reflects child mortality), 
#'      and \eqn{k} should be chosen to match 45q15 (adult mortality) or set to 0 (default). The coefficients
#'      can be passed as inputs, or taken from the package default dataset \code{\link{LQcoef}} which 
#'      are taken from \url{https://u.demog.berkeley.edu/~jrw/LogQuad/}.
#'      
#'      For the given inputs and values of life expectancy e0, the function finds values of \eqn{h} that 
#'      best match e0, using life tables and the bisection method. It returns the corresponding mortality schedule
#'      for each value of e0.
#'      
#'      Function \code{logquad} is for one sex, while \code{logquadj} can be used for both sexes.
#' @param e0 Vector of target life expectancies.
#' @param sex Which sex does the give \code{e0} corresponds to.
#' @param my.coefs Data frame with columns \dQuote{sex}, \dQuote{age}, \dQuote{ax}, \dQuote{bx}, \dQuote{cx}, \dQuote{vx}.
#'      The \dQuote{sex} column should contain values \dQuote{female}, \dQuote{male} and/or \dQuote{total}.
#'      The \dQuote{age} column must be sorted so that it assures that rows correspond to ages in increasing order.
#'      Any \code{NA}s are internally converted to zeros. If not given, the dataset \code{\link{LQcoef}} is used.
#' @param q5ranges A vector of size two, giving the min and max of 5q0 used in the bisection method.
#' @param k Value of the \eqn{k} parameter.
#' @param keep.lt Logical. If \code{TRUE} additional life table columns are kept in the 
#'      resulting object.
#' @param \dots Additional life table arguments.
#' @return Function \code{logquad} returns a list with the following elements: a matrix \code{mx}
#'     with the predicted mortality rates. If \code{keep.lt} is \code{TRUE}, it also 
#'     contains matrices \code{sr} (survival rates), and life table quantities \code{Lx} and \code{lx}.
#'     Function \code{logquadj} returns a list of objects, one for each sex.
#' @export
#' @references 
#'      Wilmoth, J., Zureick, S., Canudas-Romo, V., Inoue, M., Sawyer, C. (2012). 
#'      A Flexible Two-Dimensional Mortality Model for Use in Indirect Estimation. 
#'      Population studies, 66(1), 1-28. \doi{doi:10.1080/00324728.2011.611411}
#' @seealso \code{\link{LQcoef}}, \code{\link{mortcast.blend}}, \code{\link{mortcast}}, \code{\link{pmd}}, \code{\link{mlt}}
#' 
#' @examples
#' data(e0Mproj, package = "wpp2017")
#' country <- "Brazil"
#' # get target e0
#' e0m <- as.numeric(subset(e0Mproj, name == country)[-(1:2)])
#' # project into future
#' pred <- logquad(e0m, sex = "male")
#' # plot first projection in black and the remaining ones in heat colors 
#' plot(pred$mx[,1], type = "l", log = "y", ylim = range(pred$mx),
#'     ylab = "male mx", xlab = "Age", main = country)
#' for(i in 2:ncol(pred$mx)) lines(pred$mx[,i], 
#'     col = heat.colors(20)[i])
#'     
#' @rdname lqgroup
#' @name logquad
logquad <- function(e0, sex = c("male", "female", "total"), my.coefs = NULL,
                    q5ranges = c(1e-4, 0.9), k = 0, keep.lt = FALSE, ...) {
    get.a0rule <- function(a0rule = c("ak", "cd"), ...)
        list(ak = 1, cd = 2)[[match.arg(a0rule)]]
    sex <- match.arg(sex)
    sex.code <- list(male=1, female=2, total=3)[[sex]]
    if(is.null(my.coefs)) {
        coefs <- MortCast::LQcoef
    } else coefs <- my.coefs
    colnames(coefs) <- tolower(colnames(coefs))
    if(!all(c("sex", "age", "ax", "bx", "cx", "vx") %in% colnames(coefs)))
        stop(paste("Missing columns in the coefficient dataset.\nRequired columns are",
                   paste(c("sex", "age", "ax", "bx", "cx", "vx"), collapse = ", "), 
                   "\nAvailable columns are", paste(colnames(coefs), collapse = ", ")))
    sex.coefs <- coefs[tolower(coefs$sex) == sex,]
    sex.coefs[is.na(sex.coefs)] <- 0 # replace NA with zero, otherwise it could not be passed to the C function
    nage <- nrow(sex.coefs)
    npred <- length(e0)
    ages <- c(0, 1, seq(5, length = nage - 2, by = 5))
    a0cat <- get.a0rule(...)
    
    # initialize results
    zeromatsr <- matrix(0, nrow=nage-1, ncol=npred)
    zeromatmx <- matrix(0, nrow=nage, ncol=npred)
    result <- list(mx=zeromatmx, lx=zeromatmx, sr=zeromatsr, Lx=zeromatsr, nx = 5)

    LQres <- .C("LQuad", as.integer(npred), as.integer(sex.code), 
                 as.integer(nage), as.numeric(e0), 
                 as.numeric(sex.coefs$ax), as.numeric(sex.coefs$bx),
                 as.numeric(sex.coefs$cx), as.numeric(sex.coefs$vx),
                 Q5l=as.numeric(q5ranges[1]), Q5u=as.numeric(q5ranges[2]), 
                 K = as.numeric(k), a0rule = as.integer(a0cat),
                 LLm = as.numeric(result$Lx), Sr=as.numeric(result$sr), 
                 lx=as.numeric(result$lx), Mx=as.numeric(result$mx), PACKAGE = "MortCast")
    
    result$mx <- matrix(LQres$Mx, nrow=nage, dimnames=list(ages, names(e0)))
    if(keep.lt) {
        result$sr <- matrix(LQres$Sr, nrow=nage-1, dimnames=list(ages[-2], names(e0)))
        result$Lx <- matrix(LQres$LLm, nrow=nage-1, dimnames=list(ages[-2], names(e0)))
        result$lx <- matrix(LQres$lx, nrow=nage, dimnames=list(ages, names(e0)))
    } else {
        result$sr <- NULL
        result$Lx <- NULL
        result$lx <- NULL
    }
    return(result)
}

#' @export
#' @rdname lqgroup
#' @param e0m A time series of target male life expectancy.
#' @param e0f A time series of target female life expectancy.
#' @param \dots Additional arguments passed to the underlying function. 
#' @name logquadj
#' 
logquadj <- function(e0m, e0f, ...) {
    e0  <- list(female=e0f, male=e0m)
    res <- list()
    for(sex in names(e0)) {
        res[[sex]] <- list()
        if(!is.null(e0[[sex]]))
            res[[sex]] <- logquad(e0 = e0[[sex]], sex = sex, ...)
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
#'     Coherent Pattern Mortality Decline, Log-Quadratic model, or Model Life Tables). Weights can be applied to fine-tune the blending mix.
#' @details The function allows to combine two different methods using given weights.
#'     The weights can change over time - by default they are interpolated from the starting weight 
#'     to the end weight. As the blended mortality rates do not necessarily match the target life expectancy, 
#'     scaling is applied to improve the match, controlled by the \code{match.e0} argument. 
#'     The projection is done for both sexes, so that coherent methods can be applied.
#' @param e0m A time series of future male life expectancy.
#' @param e0f A time series of future female life expectancy.
#' @param meth1 Character string giving the name of the first method to blend. It is one of 
#'     \dQuote{lc}, \dQuote{pmd}, \dQuote{mlt} or \dQuote{logquad}, corresponding to Coherent Lee-Carter (function \code{\link{mortcast}}), 
#'      Pattern Mortality Decline (function \code{\link{copmd}}), Log-Quadratic model (function \code{\link{logquadj}}), and 
#'      Model Life Tables (function \code{\link{mltj}}), respectively. The \dQuote{logquad} method can only be used 
#'      with 5-year age groups.
#' @param meth2 Character string giving the name of the second method to blend. 
#'     One of the same choices as \code{meth1}.
#' @param weights Numeric vector with values between 0 and 1 giving the weight of \code{meth1}.
#'     If it is a single value, the same weight is applied for all time periods. 
#'     If it is a vector of size two, it is assumed these are weights for the first and the last
#'     time period. Remaining weights will be interpolated. Note that \code{meth2} is weighted 
#'     by \code{1 - weights}.
#' @param nx Size of age groups. Should be either 5 or 1.
#' @param apply.kannisto Logical. If \code{TRUE} and if any of the methods results in less than 
#'     \code{min.age.groups} age categories, the coherent Kannisto method (\code{\link{cokannisto}}) 
#'     is applied to extend the age groups into old ages.
#' @param min.age.groups Minimum number of age groups. Triggers the application of Kannisto, see above. 
#'     Change the default value if 1-year age groups are used (see Example).
#' @param match.e0 Logical. If \code{TRUE} the blended mx is scaled so that it matches the input e0.
#' @param keep.lt Logical. If \code{TRUE} additional life table columns are kept in the 
#'     resulting object. Only used if \code{match.e0} is \code{TRUE}.
#' @param meth1.args List of arguments passed to the function that corresponds to \code{meth1}. 
#' @param meth2.args List of arguments passed to the function that corresponds to \code{meth2}. 
#' @param kannisto.args List of arguments passed to the \code{\link{cokannisto}} function if Kannisto is applied. 
#'     If 1-year age groups are used various defaults in the Kannisto function need to be changed (see Example).
#' @param \dots Additional life table arguments.
#' @return List with elements \code{female} and \code{male}, each of which contains a matrix \code{mx}
#'     with the predicted mortality rates. If the result has been scaled (\code{match.e0} is \code{TRUE}), the element 
#'     \code{mx.rawblend} contains the mx before scaling. Also in such a case, if \code{keep.lt} is \code{TRUE}, it also 
#'     contains matrices \code{sr} (survival rates), and life table quantities \code{Lx} and \code{lx}.
#'     In addition, the return object contains elements \code{meth1res} and \code{meth2res}
#'     which contain the results of the functions corresponding to the two methods. 
#'     Elements \code{meth1} and \code{meth2} contain the names of the methods. 
#'     A vector \code{weights} contains the final (possibly interpolated) weights.
#' @export
#' 
#' @seealso \code{\link{mortcast}}, \code{\link{copmd}}, \code{\link{mltj}}, \code{\link{logquad}},
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
#' e0f <- subset(e0Fproj, name == country)[-(1:2)]
#' e0m <- subset(e0Mproj, name == country)[-(1:2)]
#' 
#' # Blend LC and MLT
#' pred1 <- mortcast.blend(e0m, e0f, meth1 = "lc", meth2 = "mlt",
#'     meth1.args = list(lc.pars = lcest),
#'     meth2.args = list(type = "CD_North"),
#'     weights = c(1,0.25))
#'     
#' # Blend PMD and MLT
#' pred2 <- mortcast.blend(e0m, e0f, meth1 = "pmd", meth2 = "mlt",
#'     meth1.args = list(mxm0 = mxm[, "2010-2015"],
#'                       mxf0 = mxf[, "2010-2015"]))
#'                       
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
#' # Blend LC and MLT for 1-year age groups
#' #########################################
#' # First interpolate e0 to get 1-year life expectancies (for first five years)
#' e0m1y <- approx(as.double(e0m[,1:2]), n = 5)$y
#' e0f1y <- approx(as.double(e0f[,1:2]), n = 5)$y
#' # derive toy mx in order to get some LC parameters
#' mxm1y <- mlt(seq(70, 72, length = 4), sex = "male", nx = 1)
#' mxf1y <- mlt(seq(78, 79, length = 4), sex = "female", nx = 1)
#' lcest1y <- lileecarter.estimate(mxm1y, mxf1y, nx = 1)
#' 
#' # projections
#' pred3 <- mortcast.blend(e0m1y, e0f1y, meth1 = "lc", meth2 = "mlt",
#'     weights = c(1,0.25), min.age.groups = 131, nx = 1, 
#'     meth1.args = list(lc.pars = lcest1y),
#'     kannisto.args = list(est.ages = 90:99, proj.ages = 100:130))
#'     
#' # plot results
#' par(mfrow=c(1,1))
#' plot(0:130, pred3$female$mx[,5], log = "y", type = "l", col = "red")
#' lines(0:130, pred3$male$mx[,5], col = "blue")
#' 
#' @name mortcast.blend
mortcast.blend <- function(e0m, e0f, 
                          meth1 = "lc", meth2 = "mlt", weights = c(1, 0.5), nx = 5, 
                          apply.kannisto = TRUE, min.age.groups = 28, 
                          match.e0 = TRUE, keep.lt = FALSE, 
                          meth1.args = NULL, meth2.args = NULL, kannisto.args = NULL, ...) {
    get.a0rule <- function(a0rule = c("ak", "cd"), ...)
        match.arg(a0rule)
    
    methods.allowed <- list(lc = "mortcast", mlt = "mltj", pmd = "copmd", logquad = "logquadj")
    meth1 <- match.arg(meth1, choices = names(methods.allowed))
    meth2 <- match.arg(meth2, choices = names(methods.allowed))
    
    npred <- length(e0m)
    w <- weights
    if(is.null(w)) w <- 0.5
    if(length(w) == 1) w <- rep(w, npred)
    if(!(length(w) == 2 || length(w) == npred))
        stop("Weights should be either of length 1 (constant weight), 2 (interpolated between start and end) or the same size as e0.")
    if(any(w > 1 | w < 0)) stop("Weights must be between 0 and 1.")
    
    if(!nx %in% c(1, 5))
        stop("The nx argument must be either 5 or 1.")
    
    a0rule <- get.a0rule(...)
    a0cat <- list(ak = 1, cd = 2)[[a0rule]]
    meth1.args[["nx"]] <- nx
    meth2.args[["nx"]] <- nx
    meth1.args[["a0rule"]] <- a0rule
    meth2.args[["a0rule"]] <- a0rule
    
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
    if(match.e0) {
        e0l  <- list(female=e0f, male=e0m)
        nage <- nrow(res[[sex]]$mx)
        if(nx == 5) {
            default.ages <- c(0, 1, seq(5, length = nage - 2, by = 5))
            resnage <-  nage-1 # number of age groups of the resulting matrices 
            age.groups <- default.ages[-2] # group 0-1 collapsed into 0-5
        } else {
            default.ages <- seq(0, nage - 1)
            resnage <-  nage # all ages
            age.groups <- default.ages
        }
        zeromatsr <- matrix(0, nrow=resnage, ncol=npred)
        newmx <- res[[sex]]$mx
        for(sex in names(res)) {
            ressex <- list(lx=zeromatsr, sr=zeromatsr, Lx=zeromatsr)
            adjres <- .C("adjust_mx", as.integer(npred), as.integer(c(female=2, male=1)[sex]), 
                     as.integer(nage), as.integer(nx),
                     as.numeric(res[[sex]]$mx), as.numeric(e0l[[sex]]), 
                     a0rule = as.integer(a0cat),
                     LLm = as.numeric(ressex$Lx), Sr=as.numeric(ressex$sr), 
                     lx=as.numeric(ressex$lx), Mx=as.numeric(newmx), PACKAGE = "MortCast")
            res[[sex]]$mx.rawblend <- res[[sex]]$mx
            res[[sex]]$mx[] <- adjres$Mx
            if(keep.lt) {
                res[[sex]]$sr <- matrix(adjres$Sr, nrow=resnage,
                                           dimnames=list(age.groups, names(e0l[[sex]])))
                res[[sex]]$Lx <- matrix(adjres$LLm, nrow=resnage,
                                           dimnames=list(age.groups, names(e0l[[sex]])))
                res[[sex]]$lx <- matrix(adjres$lx, nrow=resnage, 
                                           dimnames=list(age.groups, names(e0l[[sex]])))
            } else {
                res[[sex]]$sr <- NULL
                res[[sex]]$Lx <- NULL
                res[[sex]]$lx <- NULL
            }
        }
    }
    return(c(res, list(meth1res = mx1, meth2res = mx2, 
                       meth1 = meth1, meth2 = meth2, weights = w)))
}
