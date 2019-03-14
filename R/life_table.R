#' @title Life Table Function
#' @description Function for obtaining life table quantities from mortality rates.
#' @details Computes a life table corresponding to given mortality rates for 5-years age groups. 
#' @param mx Vector of age-specific mortality rates nmx. The elements correspond to 1m0, 4m1, 5m5, 5m10, \dots. 
#'    It can have no more than 28 elements which corresponds to age up to 130. 
#' @param sex For which sex is the life table.
#' @param ages Defines the age groups of \code{mx}. If \code{mx} has no names, the \code{ages} corresponds to 0, 1, 5, 10, \dots.
#' @param radix Base of the life table.
#' @param open.age Open age group. If smaller than the last age group of \code{mx}, the life table is truncated. 
#' @export
#' @examples
#' data(mxF, e0Fproj, package = "wpp2017")
#' # get female mortality of Mexico for the current year
#' country <- "Mexico"
#' mxf <- subset(mxF, name == country)[,"2010-2015"]
#' names(mxf) <- c(0, 1, seq(5, 100, by=5))
#' life.table(mxf, sex = "female")
#' 
life.table <- function(mx, sex = c("male", "female"), ages = names(mx), radix = 1, open.age = 130){
    # The first two elements of mx must correspond to 0-1 and 1-4. 
    # If include01 is FALSE, the first two age groups of the results are collapsed to 0-5
    sex <- match.arg(sex)
    sex <- list(male=1, female=2)[[sex]]
    if(is.null(ages))
        ages <- c(0, 1, seq(5, length(mx)-2))
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
    
    LTC <- .C("LifeTable", as.integer(sex), as.integer(nagem1), as.numeric(mx), 
              Lx=Lx, lx=lx, qx=qx, ax=ax, Tx=Tx, sx=sx, dx=dx)
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
    rownames(LT) <- c(paste(LT$age[-nresage], pmax(LT$age[-1]-1,1), sep="-"), paste0(LT$age[nresage], "+"))
    return(LT)
}
