## ------------------------------------------------------------------------
#' Build a  main effects list for a family of coefficient names
#'
#' @param thisFamilysCoefs  coefficients, usually for a single family (or for entire set of families for a formula)
#'
#' @return  mainEffectsList  named list of main effects with levels and attributes nfac gfacType and trueFactorNames
#' @export
#'
#' @examples
#' expandedCoefs=c("C|b",  "C|d",  "V|a",  "V|i",
#'   "V|u",  "C|b:V|a",   "C|d:V|a",  "C|b:V|i",  "C|d:V|i",  "C|b:V|u",
#'    "C|d:V|u",  "C|b:F1", "C|d:F1",  "V|a:F2",  "V|i:F2",  "V|u:F2",
#'     "C|b:cond|1",  "C|d:cond|1","C|b:cond|2",  "C|d:cond|2",
#'      "C|b:cond|3",  "C|d:cond|3")
#' thisMainEffectsList=buildFamilyMainEffectsList(expandedCoefs)
#' print(str(thisMainEffectsList))

buildFamilyMainEffectsList=function(thisFamilysCoefs){
    gfacLev=getGFacAndLevNames(thisFamilysCoefs)
    # These shold be handled ONE FAMILY at a time...
    # Construct complete main effects list
    updateMainEffectsList=function(tgfacLev,mainEffectsList){
         # browser()
        for (igfac in 1:length(tgfacLev$fac)){
            tfName=tgfacLev$fac[igfac]
            tfLev=tgfacLev$lev[igfac]
            if (length(mainEffectsList)==0|| (!tfName %in% names(mainEffectsList))){
                # New factorname encountered encode its first level
                mainEffectsList[[tfName]]=tfLev
            }else if (!(tfLev%in% mainEffectsList[[tfName]])){
                # Add another level
                mainEffectsList[[tfName]]=append(mainEffectsList[[tfName]],tfLev)
            }
        }
        return(mainEffectsList)
    } # end updateMainEffectsList
    mainEffectsList=list()
    for (icoef in 1:length(gfacLev)){
        mainEffectsList=updateMainEffectsList(gfacLev[[icoef]],mainEffectsList)
        # print(mainEffectsList)
    }
    # Decorate with attributes
    mainEffectNames=names(mainEffectsList)
    nfac=length(mainEffectNames)
    # Fix covariates (one level "factors")
    nlevsVec=sapply(mainEffectsList,length)
    gfacType=rep('unknown',nfac)
    names(gfacType)=mainEffectNames
    for (i in 1:nfac){
        # browser()
        if (nlevsVec[i]==1){
        if(!isTRUE(mainEffectsList[[i]][1] =="" )){
            warning('Expecting empty main effects list')
            browser()
        }
            #  Note we're replacing the whole vector because it was a character vector
            mainEffectsList[[i]]=1.0
            gfacType[i]='covariate'
        }else{
            gfacType[i]='factor'
        }
    }
    # print(mainEffectsList)
    # print(gfacType)
    attr(mainEffectsList,'nfac') =nfac
    attr(mainEffectsList,'gfacType')=gfacType
    attr(mainEffectsList, 'trueFactorNames') = names(gfacType[gfacType=='factor'])
    # str(mainEffectsList)
    return(mainEffectsList)
}


## ------------------------------------------------------------------------
# So now need to create a list of all combinations of all factor levels
# nlevsVec=sapply(mainEffectsList,length)
#' Construct a  factorial data table (facDtb)
#' A factorial data table (facDtb) is a   data table of main effect factor levels and filler covariates (1.0)  sufficient to cover the highest level  of factors and covariate interactions. "True factors" are identified by name in an attribute. Non-true factors are filler covariates. This is often called one family at a time.
#'
#' @param mainEffectsList  a main effects list as created by buildFamilyMainEffectsList() or a charvec of expanded coefficient names that can be made into one
#'
#' @export
#'
#' @examples
#' expandedCoefs=c("C|b",  "C|d",  "V|a",  "V|i",
#'   "V|u",  "C|b:V|a",   "C|d:V|a",  "C|b:V|i",  "C|d:V|i",  "C|b:V|u",
#'    "C|d:V|u",  "C|b:F1", "C|d:F1",  "V|a:F2",  "V|i:F2",  "V|u:F2",
#'     "C|b:cond|1",  "C|d:cond|1","C|b:cond|2",  "C|d:cond|2",
#'      "C|b:cond|3",  "C|d:cond|3")
#' famwiseDtb=constructFactorialDtb(expandedCoefs)
#' print(str(famwiseDtb))

constructFactorialDtb=function(mainEffectsList){
    if(is.character(mainEffectsList)){
        # Convert to main effects list
        mainEffectsList=buildFamilyMainEffectsList(mainEffectsList)
    }
    nlevsVec=sapply(mainEffectsList,length)
    ntot=prod(nlevsVec)
    # browser()
    factorialDataList=list()
    mainEffectNames=names(mainEffectsList)
    for (fac in mainEffectNames){
        if (length(mainEffectsList[[fac]])==1){
            # covarate
            factorialDataList[[fac]]=rep(1.0,ntot)
        }else{
            # a truefactor to be
            factorialDataList[[fac]]=rep('',ntot)
        }
    }
    # print(str(factorialDataList))
    # browser()
    nfac=length(mainEffectNames)
    for (i in 1:ntot){
        inxvec=arrayInd(i,nlevsVec)
        for (j in 1:length(inxvec)){
            # set(facDtb,i,j,mainEffectsList[[facName s[j]]][inxvec[j]])
            thisFac=mainEffectNames[j]
            # This screws up covariates again you should fix them ealier
            factorialDataList[[thisFac]][i]=mainEffectsList[[mainEffectNames[j]]][inxvec[j]]
        }
    }
    # print(str(factorialDataList))
    # browser()

    # Convert true factors  factors with proper levels
    facDtb=as.data.table(factorialDataList)
    # print(str(facDtb))
    # This is avaialble as attribute of mainEffectsList
    trueFacNames=attr(mainEffectsList,'trueFactorNames');

    for (fac in names(facDtb))
        if( fac %in% trueFacNames){

            facDtb[ , (fac):= factor(facDtb[[fac]], levels=mainEffectsList[[fac]])
                    ]
        }
    # copy attributes
    attr(facDtb,'trueFactorNames')=trueFacNames
    # print(str(facDtb))
    # browser()
    return(facDtb)
}


## ------------------------------------------------------------------------
#' Make family-wise data matrix from a facDtb
#'
#' @param facDtb
#'
#' @return  a  stimulus-response MNL-compatible model.matrix
#' @export
#'
#' @examples
#' print('no example')
makeFamilyWiseDmxFromFactorialDtb=function(facDtb) {
    #  Note this AUTOMATICALLY determines the formula
    #  from the list of true factors.
    #  The formula is the "*" combiation of the true factors.
    #  However we eliminate all but the highest order terms
    #  as they have presumably been dealt with elsewhere
    # we should check this somewhere
    # Note teh column names will not include any covariate terms only true factors.
    # This will have to be dealt with in creating the bigMatrix
    mmContrastList=list()
    trueFacNames=attr(facDtb,'trueFactorNames')
    formStr=paste0("~", paste0(trueFacNames,collapse='*'))
    # browser()
    stopifnot(!is.null(trueFacNames))
    for (sv in trueFacNames){
        # Compile the contrast matrix
        mmContrastList[[sv]]="mySumContrFun"
    }
    mmFac=model.matrix(as.formula(formStr),data=facDtb,contrasts.arg = mmContrastList)
    # str(mmFac)
    # remove everyghing but the highest order interaction cols.. that's determiable simply from number of ":" in colnamesru
    # it's one less than number of true factors
    ncolon=length(trueFacNames)-1;
    inxHiOrd=which(grepl(':',colnames(mmFac))==ncolon)
    # Still need to remove intercept for main effects (no colons)
    inxHiOrd=inxHiOrd[inxHiOrd>1]
    return(mmFac[,inxHiOrd,drop=FALSE])
}


## ------------------------------------------------------------------------
#' Create a coefficient expander list
# ' Create objects to map sum-to-zero  contrasts->  redundant sum-to-zero  indicator coefs
#'
#' @param expandedCoefNames - character vector of expanded coef names (available from e.g. make_rsSumDmxFrom_srDtb object) a e.g. c("C|d", "V|a",... "C_d:V_a: ...., C|b:V|a:y")

#' @param shouldReturnBigMatrix  default FALSE. If true, RM_bigMat component returned

#' @return  result: a list with components familyNames, XM_list and optionally XM_bigMat

#' @details
#'  Given an character vector of expanded (indicator) coefficient names for a HMNL, this function creates a result list consisting of
#     result$familyNames the family names of groups of coefficients. EG
#    if coefNames == c("C|d", "V|a",... "C_d:V_a: ...., C|b:V|a:y"), the returned
#'    familyNames would include "C", "V", C:V", "C:V:y"
#'    result$XM_list a list of expansion matrices for expanding sum-to-zero k-1
#' levels  per factor contrast coefficients  to sum-to-zero indicator coefficients
#' (all k levels per factor)
#' and  optionally result$XM_bigMat  a large (ncoefExpanded x ncoefContrast) matrix with non-zero blocks constisting of
#  result$RM_list[[ifam]] ifam=1:nfam.

#' @param shouldMakeBigMatrix  - should a big  exoander matrix be created or just the list of component matrices
#'
#' @return result : a list with conponents:
#'  $familyNames the family names of groups of coefficients. EG
#'  $XM_list a list of expansion matrices for expanding sum-to-zero k-1
#'          levels  per factor contrast coefficients  to sum-to-zero indicator coefficients
#'          (all k levels per factor) Note the column names will contain only factors not covaraites!
#'    $XM_bigMat  (optionally) a large (ncoefExpanded x ncoefContrast) matrix with non-zero blocks constisting of
#        result$RM_list[[ifam]] ifam=1:nfam. The "bigMat" will have the correct coefficient names of the sum2zero coefs
#        The column names will be empty as they can't be reconstructed easily from the coef names only.
#        They should correspond to a set of make_rsSumDmx column names compatible with the genreating formula.
#        This is done in make_rsSumDmxFrom_srDtb()
#        This is a fragile implementation.
#' @export
#'
#' @examples
#' tlist=createCoefficientExpanderList(c("C|b",  "C|d",  "V|a",  "V|i",
#'   "V|u",  "C|b:V|a",   "C|d:V|a",  "C|b:V|i",  "C|d:V|i",  "C|b:V|u",
#'    "C|d:V|u",  "C|b:F1", "C|d:F1",  "V|a:F2",  "V|i:F2",  "V|u:F2",
#'     "C|b:cond|1",  "C|d:cond|1","C|b:cond|2",  "C|d:cond|2",
#'      "C|b:cond|3",  "C|d:cond|3"), shouldMakeBigMatrix=TRUE)
#' print(str(tlist))
createCoefficientExpanderList=function(expandedCoefNames,shouldMakeBigMatrix=TRUE){
    famWiseExpandedCoefList=  getFamilyWiseCoefList(expandedCoefNames)
    familyNames=names(famWiseExpandedCoefList)
    # browser()
    # Quick test:
    # Pre multiply by big ResidualMaker matrix should yield equivalent resul
    # Expand the coefficients onemfamily at a time
    coefficientExpanderList=list()
    nExpandedCoef=length(expandedCoefNames)
    # Can/must we assume that famWiseCoefList will be in order?
    # Go thru families
    if (shouldMakeBigMatrix){
        # Make this too wide, we'll narrow it later
        XM_bigMat=matrix(0,nExpandedCoef,nExpandedCoef)
        rownames(XM_bigMat)=expandedCoefNames
        colNamesBigMat=c();
    } else {
        XM_bigMat=NULL
    }
    maxColBigMat=0;
    inxSumContrOffset=0
    XM_list=list()
    for (fam in familyNames){
        # catln(fam)
        # famTestCoefNames=famWiseTestCoefList[[fam]]
        # inxTest=which(famTestCoefNames%in%testCoefNames)
        #
        # browser()
        famExpandedCoefNames=famWiseExpandedCoefList[[fam]]
        thisMainEffectsList=buildFamilyMainEffectsList(famExpandedCoefNames)
        # print(str(thisMainEffectsList))
        thisFamDtb= constructFactorialDtb(thisMainEffectsList)
        # print(str(thisFamDtb))
        famExpanderDmx= makeFamilyWiseDmxFromFactorialDtb(thisFamDtb)

        rownames(famExpanderDmx)=famExpandedCoefNames
        XM_list[[fam]]=thisFamDtb
        # print(famExpanderDmx)
        # print(str(famExpanderDmx))
        if (shouldMakeBigMatrix){
            #  Order is crucial expanded must be first which(1:5 %in% c(3:5) )
            inxExpanded=which(expandedCoefNames %in% famExpandedCoefNames)
            nFamSumContrCoefs=ncol(famExpanderDmx)
            inxSumContrCoefs=inxSumContrOffset+(1:nFamSumContrCoefs)
            inxSumContrOffset=inxSumContrOffset+nFamSumContrCoefs
            maxColBigMat=max(inxSumContrCoefs)
            XM_bigMat[inxExpanded,inxSumContrCoefs]=famExpanderDmx[, , drop=FALSE]
            # colNamesBigMat=c(colNamesBigMat,colnames(famExpanderDmx))
            # catln('colNamesBigMat')
            # print(colNamesBigMat)
            # browser()
            # print(fam)
            # catln('inxSumContrCoefs (cols)',inxSumContrCoefs)
            # catln('inxExpanded (rows)',inxExpanded)
            # catln("BigMat so far")
            # print(XM_bigMat[ 1:max(inxExpanded), 1:max(inxSumContrCoefs),drop=FALSE])
            # print(famExpanderDmx)
            # browser()
        } # endif should makeBigMat
    }
    if(shouldMakeBigMatrix){
        # prune and name columns
        XM_bigMat=XM_bigMat[,1:maxColBigMat]
        colnames(XM_bigMat)=NULL
    }
    return(list(familyNames=familyNames, XM_list=XM_list, XM_bigMat=XM_bigMat))
}


## ------------------------------------------------------------------------
#' Check marginal sum-to-zero constraints of expanded coefficients
#'
#' @param expandedCoefs  a  [K x 1] matrix of named {fac|lev:} coefficients (expanded) that obey marginal sum-to-zero
#'
#' @return logical (TRUE if all sum to zero constraints are met)
#' @export
#'
#' @examples
#' print ('No examples see vignett TestExpandSumToZeroCoefs.Rmd')
checkMarginalsSumToZero=function(expandedCoefs){
    # Check that each family subset sums to zero globally and marginally
    # The logic is simple, the programming tricky.
    # A whole family of coefficients should sum to zero.
    # For each factor k (column k of inxMat), every set of coefficients sharing all the same levels on all other (not-k) factors should add up to zero (along factor k)
    #
    checkMargMatSumToZero = function(){
        sumNearZero=function(x){
            abs(mean(x))<sqrt(.Machine$double.eps)}
        # browser()
        nCol=ncol(inxMat)
        chk= sumNearZero(famCoefs)
        if (!chk){
            message('Not global sumtozero')
            browser()
            return(FALSE)
        }
        # if a single main effect factor, only the global sum matters
        if (nCol==1){
            return(chk)
        }
        # Otherwise, check each dimension (column at a time)
        inxDim=1:nCol
        for (k in 1:nCol){
            #  inxMatM1 "index matrix minus 1 column"
            #   indices of all factors except k
            inxMatM1=inxMat[ , inxDim[-k], drop=FALSE]
            # Unique combinations of all non-k factors
            # (marginal to column k) patterns
            uInxMatM1=unique(inxMatM1, drop=FALSE)
            # Check that coefs sharing each unique marginal pattern sum to zero
            for (iur in 1:nrow(uInxMatM1)){
                binxMarg=apply(inxMatM1, 1, function(x) all(x==uInxMatM1[iur,]))
                # catln('k', k, 'uInxMatM1[iur,]', uInxMatM1[iur,])
                # catln('binxMarg',binxMarg, 'which binxMarg', which(binxMarg))
                tsum = sum(famCoefs[binxMarg,1])
                tchk = sumNearZero(famCoefs[binxMarg,1])
                if (!tchk){
                    message('chk failed on column', k, 'tsum', tsum)
                    browser()
                }
                chk=chk&tchk
                # catln('k,tchk, chk, tsum', k,tchk,chk,tsum)
            } # end for iur
        } # end for k
        # catln('Check all', chk)
        return(chk)
    }# end  checkMargMatSumToZero()

    if (!is.matrix(expandedCoefs) || dim(expandedCoefs)[2]!=1){
        stop('Arg expandedCoers should be a column vector')
    }
    famWiseCoefs=getFamilyWiseCoefList(rownames(expandedCoefs))
    checkAllFams=TRUE
    for (fam in names(famWiseCoefs)){
        famCoefNames=famWiseCoefs[[fam]]
        # 1 needed expandedCoefs is col vector
        famCoefs=expandedCoefs[famCoefNames, 1, drop=FALSE]
        # browser()
        facDtb=constructFactorialDtb(famCoefNames)
        # print(str(facDtb))
        trueFactors=attr(facDtb,'trueFactorNames')
        facDtbFacOnly=facDtb[,trueFactors,with=FALSE]
        # print(str(facDtbFacOnly))
        # Convert to numeric matrix sapply works fine
        inxMat=sapply(facDtbFacOnly,as.numeric)
        # print(inxMat)
        # print(str(inxMat))
        checkThisFam=checkMargMatSumToZero()
        if (!checkThisFam){
            catln('Fam', fam, 'did not sum to zero on all margins')
        }
        checkAllFams=checkThisFam & checkAllFams
    }
    return(checkAllFams)
}



