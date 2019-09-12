#' Test data set for Bayesian factorial HMNL logistic.
#'
#' Synthetic data
#'
#' @format  A list with four conponenrts  "srDtbAll" "thetaMatTrue" "bMatSjList"   "lambda_sList"
#' \describe{
#'   \item{srDtbAll} {  First list component a data table with following fields}
#'   \item{Sj}{Factor chosenb syllable (one of ba,da,bi,di,bu,du)}
#'   \item{Catsym} {character one of A B C}
#'   \item{onez}{ a column of ones}
#'   \item{x}{strandardized x values}
#'   \item{thetaMatTrue} { second element =a matrix [nVar x nCat]  of underlying population theta}
#'   \item{bMatSjList} { third element =a matrix [nSj x nVar]  of underlying beta+theta coefs for subjects}
#'   \item{lambda_sList} {fourth el. list [[ nSj]] with comps "lambda_s" "lambda" "s_u" of underlying lapse parameters}
#'   }
#' @source {Made up using vignette CreateTest1D3Cat_HMNL_data.Rmd}
# test1D3CatLapseDatList
