# Remove all but pars from stanfit object
# from https://github.com/stan-dev/rstan/issues/426#issue-239658968
# Modified by tmn to check for stanfit object
#' Trim a stanfit object to only wanted names
#'
#' @param object stanfit object to prune
#' @param pars a charvec of parameter names to retain (lp_ is always included)
#'
#' @return returns a trimmed stanfit object
#' @export
#'
#' @examples
#' print('No example')
# cleanObject <- function(object,pars){
 pruneStanfitObj <- function(object,pars){
    pars <- c(pars,'lp__')
    # tmn next line:
    if( !("stanfit" %in% class(fit1NC) )) stop("First arg should be stanfit")
    nn <- paste0('^',pars,'(\\[|$)',collapse="|")
    ids <-  grep(nn,  object@sim$fnames_oi)
    # browser()
     ids.2 <- which(names(object@par_dims) %in% pars)
    for(i in 1:4){
        a <- attributes(object@sim$samples[[i]])
        x <- object@sim$samples[[i]][ids]
        for(j in c('names','inits','mean_pars'))
            a[[j]] <- a[[j]][ids]
        attributes(x) <- a
        object@sim$samples[[i]] <- x
    }
    object@par_dims <- object@par_dims[ids.2]
    object@sim$dims_oi <-   object@sim$dims_oi[ids.2]
    object@sim$pars_oi<- object@sim$pars_oi[ids.2]
    object@sim$fnames_oi <-  object@sim$fnames_oi[ids]
    object@sim$n_flatnames <- length(object@sim$fnames_oi)
    object
}
