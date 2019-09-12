#----napp2logistic()----
#' Convert a (equi-covariance) NAPP model to logistic regression coefficients
#'
#' @param meanColVecs  [nVar x nGroup] matrix of  means
#' @param covMat [nVar x nVar] common covariance matrix
#' @param priors [nGroup x 1] vector of priors
#'
#' @return bMat [(nVar+1) x 1] vector of coefficients s.t. xmatAug %*% b yields logistic scores
#'  and xmatAug rbind(1, xmat)  where xmat is rbind(xvar1,xvar2,xvar3...)
#' @export
#'
#' @examples
#'   mu=matrix(c(-2, 0, 2,
#' -2, 0, -2),2,3, byrow = FALSE);
#' covMat=matrix(c(1, .5,
#'                .5, 1), 2,2, byrow=FALSE);
#' bMat=napp2logistic(meanColVecs,covMat))
#' print(bMat)
napp2logistic = function (meanColVecs,covMat,priors=NULL){
    # browser()
    # # Convert equicovariance NAPP model coefs to equivalent logistics
    # # function [b0,b]=napp2logistic(meanColVecs,covMat,priors);
    # # INPUT
    # # 1) meanColVecs  --   nVar-by-nGroup matrix of mean coefficients
    # # 2) covMat  --   nVar-by-nVar covariance matrix
    # # 3) priors  --   nGroup-by-1 vector of prior probabilities
    # # OUTPUT
    # # 1) bMat --(nVar+1) x nGroup matrix intercepts (first column is intercepts)
    # # Copyright (c) T M Nearey 2012
    # # Version 1.1 12-Jul-2012
    # ## *** SEE ALSO ***
    # # linClassifier2logisticCoefs
 if (!( is.matrix(meanColVecs) && is.matrix(covMat))){
     stop('meanColVecs. and covMat must both be matrices.')
 }
 nVar=nrow(meanColVecs)
 nGroup = ncol(meanColVecs)
 if (! all(dim(covMat)==nVar)){
     stop('covMat must be square matrix  with same number of rows as meanColVecs')
 }
    isempty=function(x) is.null(x) || length(x)==0
    if (isempty(covMat)){
        covMat = error('Must have covMat input'); #  # default
    }



    if (is.null(priors)){
        priors=1.0/nGroup*rep(1,nGroup)
    }

    if (length(priors) != nGroup){
        stop('Priors must be same length as nGroup')}
    priors=matrix(priors)


    #INPUT
    # meeanColVecs(nVar x nGroup)
    #      Given mean vectors of k groups (one per column of meanColVecs,of  length nVar)
    # covMat
    #	Given pooled covariance matrix,
    # calculate logistics of corresponding NAPP models
    # Output
    # b0  are 1 x nGroup intercepts
    # b  nVar x ngroups coefficients of variables
    #
    # Note if x is a row vector of stimuli,
    # APP(x,g)= exp(b0(g)+x*b(:,g))/sum(exp(b0+x*b))
    #  See BMDP 1981 p 680

    invcovmat=MASS::ginv(covMat)
    lnpriors=log(priors)
    #	beta1(:,gcr.useg)=invcovmat*gcr.mu;
    # 	alpha1(gcr.useg)=gcr.lnprior(gcr.useg)-.5*gcr.mu(:,gcr.useg)'*invcovmat*gcr.mu(:,gcr.useg);
    bMat=matrix(NA,nVar+1,nGroup)
    for( i in 1:nGroup){
        bMat[-1,i]=invcovmat %*% meanColVecs[ ,i]
        bMat[1,i]=lnpriors[i]-.5*t(meanColVecs[ ,i,drop=FALSE]) %*% invcovmat %*% meanColVecs[ ,i,drop=FALSE]
    }
 return(bMat)
}



