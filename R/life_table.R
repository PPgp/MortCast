#' @title Life Table Function
#' @description Function for obtaining life table quantities from mortality rates.
#' @details Computes a life table corresponding to given mortality rates for either 5- or 1-year age groups. The implementation follows
#'    Preston et al. (2001), including the choice of ax (see Table 3.3 on page 48). 
#' @param mx Vector of age-specific mortality rates nmx. If \code{abridged} is \code{TRUE} (default), 
#'    the elements correspond to 1m0, 4m1, 5m5, 5m10, \dots. 
#'    If \code{abridged} is \code{FALSE}, they correspond to 1m0, 1m1, 1m2, 1m3, \dots.
#' @param sex Which sex the mortality rates correspond to. 
#' @param abridged Is it an abridged life table (\code{TRUE}, default) or not (\code{FALSE}). 
#'    In the former case, the \code{mx} vector is interpreted as corresponding to age groups 0, 1-4, 5-9, 10-14, \dots.
#'    If \code{FALSE}, the \code{mx} vector is interpreted as corresponding to one-year age groups, i.e. 0, 1, 2, 3, \dots.
#' @param a0rule Rule for approximation of a0. "ak" (default) uses the Andreev-Kingkade method (DR 2015), "cd" uses the 
#'    Coale-Demeany method.  
#' @param radix Base of the life table.
#' @param open.age Open age group. If smaller than the last age group of \code{mx}, the life table is truncated. 
#'    It does not have any effect if larger than the last age group.
#' @return Data frame with rows corresponding to age groups and the following columns:
#'    \describe{
#'       \item{age}{Starting year of the age group.}
#'       \item{mx}{Age-specific mortality rates as passed into the \code{mx} argument.}
#'       \item{qx}{Probability of dying between ages x and x+n.}
#'       \item{lx}{Number of survivors at age x.}
#'       \item{dx}{Number of deaths between ages x and x+n.}
#'       \item{Lx}{Person-years lived between ages x and x+n.}
#'       \item{sx}{Survival rate from age x to x+n. Note that in an abridged life table, sx always refers to 5-year intervals. 
#'                 Here, sx in the first row is the survival from births to the second age group, sx in the second row 
#'                 is the survival from age 0-4 to age 5-9, third row has the survival from 5-9 to 10-14 etc.}
#'       \item{Tx}{Person-years lived after age x.}
#'       \item{ex}{Life expectancy at age x.}
#'       \item{ax}{Average person-years lived in the interval by those dying in the interval. For young ages, it follows Preston et al. (2001), Table 3.3 on page 48.
#'                 For compatibility with computations done at the UN, we set ax for ages 5 and 10 in the abridged version
#'                 to 2.5. For an unabridged life table, ax is set to 0.5 for all but first and last age groups.}
#' }
#' @references 
#'    Preston, S.H., Heuveline, P., Guillot, M. (2001). Demography: Measuring and Modeling Population Processes. Oxford: Blackwell Publishers Ltd.
#' @export
#' @examples
#' data(mxF, e0Fproj, package = "wpp2017")
#' # get female mortality of Mexico for the current year
#' country <- "Mexico"
#' mxf <- subset(mxF, name == country)[,"2010-2015"]
#' life.table(mxf, sex = "female")
#' 
life.table <- function(mx, sex = c("male", "female", "total"), abridged = TRUE, a0rule = c("ak", "cd"), 
                       radix = 1, open.age = 130){
    # If abridged is TRUE, the first two elements of mx must correspond to 0-1 and 1-4. 
    sex <- match.arg(sex)
    sex <- list(male=1, female=2, total=3)[[sex]]
    a0cat <- list(ak = 1, cd = 2)[[match.arg(a0rule)]]
    if(abridged) {
        ages <- c(0, 1, seq(5, length = length(mx)-2, by = 5))
        LTfct <- "LifeTableAbridged"
    } else {
        ages <- seq(0, length = length(mx))
        LTfct <- "LifeTable"
    }
    nage <- length(ages)
    Lx <- lx <- qx <- Tx <- sx <- dx <- ax <- rep(0, nage)
    nagem1 <- nage-1
    resage <- as.integer(ages)
    resage <- resage[resage <= open.age]
    nresage <- length(resage)
    if(!all(c(0,1) %in% resage)) stop("Ages 0 and 1 must be included in the mortality age groups.")
    
    nas <- rep(NA,nresage)
    if(any(is.na(mx))) # there are NAs in mx
        return(data.frame(age=resage, mx=mx[1:nresage], qx=nas, lx=nas, dx=nas, Lx=nas, 
                          sx=nas, Tx=nas, ex=nas, ax=nas))
    
    LTC <- do.call(".C", list(LTfct, as.integer(sex), as.integer(nagem1), as.numeric(mx), as.integer(a0cat),
              Lx=Lx, lx=lx, qx=qx, ax=ax, Tx=Tx, sx=sx, dx=dx, PACKAGE = "MortCast"))
    LT <- data.frame(age=as.integer(ages), mx=mx, qx=LTC$qx, lx=LTC$lx, dx=LTC$dx, Lx=LTC$Lx,  sx=LTC$sx, Tx=LTC$Tx, 
                     ex=LTC$Tx/LTC$lx, ax=LTC$ax)
    LT$ax[nage] <- LT$ex[nage]

    if(radix != 1) 
        LT <- transform(LT, lx = lx * radix, 
                        dx = dx * radix,
                        Lx = Lx * radix,
                        Tx = Tx * radix
        )
    if(open.age < ages[nage]) {
        ## truncate life table with open age group < 130+ (Patrick Gerland's code)
        LT <- LT[1:nresage,]
        ## mx for open age /* reciprocal of 1/ex for open age group*/ 
        LT$mx[nresage] <- 1 / LT$ex[nresage]
        LT$qx[nresage] <- NA
        LT$dx[nresage] <- LT$lx[nresage] ## open age group dx = lx
        LT$Lx[nresage] <- LT$Tx[nresage] ## open age group Lx = Tx
        # Sx
        # for open age group, e.g., 85+ 
        ## penultimate age group -> Last entry of S(x,n) is S( 80+,5) = T( 85) / T( 80)
        ## for open age group itself: Sx cannot be computed due to trunaction Sx <- NA
        LT$sx[nresage-1] <- LT$Tx[nresage]/LT$Tx[nresage-1]
        LT$sx[nresage] <- NA
        LT$ax[nresage] <- LT$ex[nresage] ## for open age group ax = ex
    }
    if(abridged)
        rownames(LT) <- c(paste(LT$age[-nresage], pmax(LT$age[-1]-1,1), sep="-"), paste0(LT$age[nresage], "+"))
    else rownames(LT) <- LT$age
    return(LT)
}
