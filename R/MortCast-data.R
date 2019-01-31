#' @title Model Life Tables Lookup
#' @docType data
#' @description Lookup table containing values for various model life tables, including 
#' Coale-Demeny and UN life tables.
#' 
#' @usage data(MLTlookup)
#'
#' @format Data frame with the following columns:
#' \describe{
#' \item{type}{Type of the model life table. Available options are \dQuote{CD East}, \dQuote{CD North}, 
#' \dQuote{CD South}, \dQuote{CD West}, \dQuote{UN Chilean}, \dQuote{UN Far_East_Asian}, 
#' \dQuote{UN General}, \dQuote{UN Latin}, \dQuote{UN South_Asian}. 
#' For the CD types, see Coale et al. (1983). For the UN types, see \url{https://population.un.org/wpp/Download/Other/MLT}.
#' }
#' \item{sex}{Code for distinguishing sexes. 1 is for male, 2 is for female.}
#' \item{age}{Starting age of an age group. These are 0, 1, 5, 10, ... 130.}
#' \item{e0}{Level of life expectancy, starting at 20 and going by steps of 2.5 up to 105.}
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
#' \url{https://population.un.org/wpp/Download/Other/MLT}
#' 
#' @examples 
#' data(MLTlookup)
#' str(MLTlookup)
#' # CD West life table for male at e0 of 80
#' subset(MLTlookup, type == "CD West" & sex == 1 & e0 == 80)
#' 
#' @keywords datasets
"MLTlookup"


#' @title Pattern Mortality Decline Lookup Tables
#' @name rhoPMD
#' @aliases rhoPMD RhoFemales RhoMales
#' @docType data
#' 
#' @description Data object containing two tables with \eqn{\rho} coefficients for the 
#' Pattern Mortality Decline method as implemented in the \code{\link{pmd}} function.
#' 
#' @usage data(rhoPMD)
#'
#' @format Using \code{data(rhoPMD)} loads two objects into memory: \code{RhoFemales} and
#'     \code{RhoMales}. They both are data frames with 22 rows corresponding to age groups, 
#'     and 17 columns corresponding to different levels of life expectancy in 5-years intervals 
#'     (from 50 to 135). The names of the columns reflect the middle of the respective interval.
#' 
#' 
#' @seealso \code{\link{pmd}}
#' 
#' @references 
#' Andreev, K. Gu, D., Gerland, P. (2013). Age Patterns of Mortality Improvement by Level of Life Expectancy at Birth with Applications to Mortality Projections. Paper presented at the Annual Meeting
#' of the Population Association of America, New Orleans, LA. \url{http://paa2013.princeton.edu/papers/132554}.
#' 
#' Gu, D., Pelletier, F. and Sawyer, C. (2017). Projecting Age-sex-specific Mortality: A Comparison of the Modified Lee-Carter and Pattern of Mortality Decline Methods, UN Population Division, 
#' Technical Paper No. 6. New York: United Nations. \url{https://population.un.org/wpp/Publications/Files/WPP2017_TechnicalPaperNo6.pdf}
#' 
#' @examples 
#' data(rhoPMD)
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
#' @rdname rhoPMD
NULL



