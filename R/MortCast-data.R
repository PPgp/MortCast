#' @title Model Life Tables Lookup
#' @name MLTlookup
#' @aliases MLTlookup MLT1Ylookup
#' @docType data
#' @description Lookup tables containing values for various model life tables, including 
#' Coale-Demeny and UN life tables.
#' 
#' @usage data(MLTlookup)
#' data(MLT1Ylookup)
#'
#' @format Data frame with the following columns:
#' \describe{
#' \item{type}{Type of the model life table. Available options are \dQuote{CD_East}, \dQuote{CD_North}, 
#' \dQuote{CD_South}, \dQuote{CD_West}, \dQuote{UN_Chilean}, \dQuote{UN_Far_Eastern}, 
#' \dQuote{UN_General}, \dQuote{UN_Latin_American}, \dQuote{UN_South_Asian}. 
#' For the CD types, see Coale et al. (1983). For the UN types, see the link in References below.
#' }
#' \item{sex}{Code for distinguishing sexes. 1 is for male, 2 is for female.}
#' \item{age}{Starting age of an age group. In \code{MLTlookup} these are 0, 1, 5, 10, ... 130. 
#'            The \code{MLT1Ylookup} table contains 1-year ages ranging from 0 to 130.}
#' \item{e0}{Level of life expectancy, starting at 20 and going by steps of 2.5 up to 115.}
#' \item{mx}{Mortality rates.}
#' \item{lx, Lx, sx}{Other life table columns.}
#'}
#'
#' @seealso \code{\link{mlt}}
#' 
#' @references 
#' Coale, A., P. Demeny, and B. Vaughn. 1983. Regional model life tables and stable 
#' populations. 2nd ed. New York: Academic Press.
#' 
#' \url{https://www.un.org/development/desa/pd/data/extended-model-life-tables}
#' 
#' @source An updated version of these datasets were provided by Sara Hertog, United Nations Population Division, in October 2021 (package version >= 2.6-0).
#'      For previous version of the tables, install MortCast 2.5-0: `devtools::install_github("PPgp/MortCast@v2.5-0")`
#' 
#' @examples 
#' data(MLTlookup)
#' str(MLTlookup)
#' # CD West life table for male at e0 of 80
#' subset(MLTlookup, type == "CD_West" & sex == 1 & e0 == 80)
#' 
#' @keywords datasets
"MLTlookup"


#' @title Pattern Mortality Decline Lookup Tables
#' @name PMDrho
#' @aliases PMDrho RhoFemales RhoMales
#' @docType data
#' 
#' @description Data object containing two tables with \eqn{\rho} coefficients for the 
#' Pattern Mortality Decline method as implemented in the \code{\link{pmd}} function.
#' 
#' @usage data(PMDrho)
#'
#' @format Using \code{data(PMDrho)} loads two objects into memory: \code{RhoFemales} and
#'     \code{RhoMales}. They both are data frames with 22 rows corresponding to age groups, 
#'     and 17 columns corresponding to different levels of life expectancy in 5-years intervals 
#'     (from 50 to 135). The names of the columns reflect the middle of the respective interval.
#' 
#' 
#' @seealso \code{\link{pmd}}
#' 
#' @references 
#' Andreev, K. Gu, D., Gerland, P. (2013). Age Patterns of Mortality Improvement by Level of Life Expectancy at Birth with Applications to Mortality Projections. Paper presented at the Annual Meeting
#' of the Population Association of America, New Orleans, LA. \url{https://paa2013.princeton.edu/papers/132554}.
#' 
#' Gu, D., Pelletier, F. and Sawyer, C. (2017). Projecting Age-sex-specific Mortality: A Comparison of the Modified Lee-Carter and Pattern of Mortality Decline Methods, UN Population Division, 
#' Technical Paper No. 6. New York: United Nations. \url{https://population.un.org/wpp/Publications/Files/WPP2017_TechnicalPaperNo6.pdf}
#' 
#' @examples 
#' data(PMDrho)
#' head(RhoFemales)
#' head(RhoMales)
#' 
#' # plot a few male patterns
#' e0lev <- colnames(RhoMales)[c(1, 5, 9, 13, 17)]
#' plot(RhoMales[, e0lev[1]], type="l", log="y", ylim=range(RhoMales[,e0lev]),
#'     ylab="male rho", xlab="Age")
#' for(i in 2:length(e0lev)) lines(RhoMales[,e0lev[i]], lty = i)
#' legend("bottomleft", legend = e0lev, lty = 1:length(e0lev), bty= "n")
#' 
#' 
#' @keywords datasets
#' @rdname PMDrho
NULL

#' @title Coefficients for Sex Ratio Adjustments in the PMD Method
#' @name PMDadjcoef
#' @aliases PMDadjcoef
#' @docType data
#' 
#' @description Data object containing a table of coefficients to be used to adjust the sex ratio in the 
#' coherent Pattern Mortality Decline method as implemented in the \code{\link{copmd}} function. To invoke 
#' the adjustment, argument \code{sexratio.adjust} should be set to \code{TRUE}.
#' 
#' @usage data(PMDadjcoef)
#'
#' @format Data frame containing columns \dQuote{age},  \dQuote{intercept}, \dQuote{lmxf}, \dQuote{e0f}, \dQuote{e0f2}, and \dQuote{gap}. 
#'     Rows correspond to age groups. The values are estimates of the following  regression 
#'     \deqn{\log_{10} mx^M = \beta_0 + \beta_1\log_{10} mx^F + \beta_2 e_0^F + \beta_3 (e_0^F)^2 + \beta_4(e_0^F - e_0^M)}
#'     The order of the columns starting with intercept corresponds to the order of the coefficients in the above equation.
#'
#' @source The coefficients were estimated and provided by Danan Gu, UN Population Division.
#' 
#' @seealso \code{\link{copmd}}
#' 
#' @references 
#' Gu, D., Pelletier, F. and Sawyer, C. (2017). Projecting Age-sex-specific Mortality: A Comparison of the Modified Lee-Carter and Pattern of Mortality Decline Methods, UN Population Division, 
#' Technical Paper No. 6. New York: United Nations. \url{https://population.un.org/wpp/Publications/Files/WPP2017_TechnicalPaperNo6.pdf}
#' 
#' @examples 
#' data(PMDadjcoef)
#' PMDadjcoef
#' 
#' 
#' @keywords datasets
"PMDadjcoef"

#' @title Coefficients for the Log-Quadratic Mortality Model
#' @name LQcoef
#' @aliases LQcoef
#' @docType data
#' 
#' @description Data object containing a table of coefficients to be used in the Log-Quadratic Model 
#'     as implemented in the \code{\link{logquad}} function. 
#' 
#' @usage data(LQcoef)
#'
#' @format Data frame containing columns \dQuote{sex}, \dQuote{age}, \dQuote{ax}, \dQuote{bx}, \dQuote{cx}, \dQuote{vx}. 
#'     Rows correspond to sex and age groups. 
#'
#' @source \url{https://u.demog.berkeley.edu/~jrw/LogQuad/}
#' 
#' @seealso \code{\link{logquad}}
#' 
#' @references 
#'      Wilmoth, J., Zureick, S., Canudas-Romo, V., Inoue, M., Sawyer, C. (2012). 
#'      A Flexible Two-Dimensional Mortality Model for Use in Indirect Estimation. 
#'      Population studies, 66(1), 1-28. \doi{doi:10.1080/00324728.2011.611411}
#' 
#' @examples 
#' data(LQcoef)
#' head(LQcoef)
#' 
#' 
#' @keywords datasets
"LQcoef"

