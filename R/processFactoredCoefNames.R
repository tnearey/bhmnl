## ------------------------------------------------------------------------
# library(stringr)
# library(nutils)
# Value
# For str_match, a character matrix. First column is the complete match, followed by one column for each capture group. For match_all, a list of character matrices.

catln <- function(...) {cat(...,'\n')}


## ------------------------------------------------------------------------
#' getTestCoefNames
#'
#' @return list of test coefficient names
#' @export
#'
#' @examples
#' print(' no examples')
getTestCoefNames=function(){
  coefNames=c(
    "C|b",
    "C|d",
    "V|a",
    "V|i",
    "V|u",
    "C|b:V|a",
    "C|d:V|a",
    "C|b:V|i",
    "C|d:V|i",
    "C|b:V|u",
    "C|d:V|u",
    "C|b:F1",
    "C|d:F1",
    "V|a:F2",
    "V|i:F2",
    "V|u:F2",
    "C|b:cond|1",
    "C|d:cond|1",
    "C|b:cond|2",
    "C|d:cond|2",
    "C|b:cond|3",
    "C|d:cond|3"
  )
}


## ------------------------------------------------------------------------
#' getGFacAndLevNames
#' Get general factor and factor level names
#`  @details
#`  gfactors are "general factors" that include categorical factors  (cfactors)  and
#`  continuous covariates. The former have levels, the latter don;t
#`
#' @param theseCoefs  a character vector of coefficient names of form <gfacName>("|"<levelName>)
#' separated by ":"
#' @return  returns a gfacLevList with components $fac and $lev
#'  $levels is NA for  covariate
#' @export
#' @examples
#' print('noexamples')

getGFacAndLevNames=function(theseCoefs){

  # Extract factor names and levels from coefficient names
  #  But ignore factors without levels (covariates)
  # Looking
  # Regex 101 '([^|:]+)(?:[|])?([^:]*)(?:[:])?
  # match 1 whole pattern.
  # match 2 ([^|:]+) factor;  string of chars not inluding | or :
  # non matching group maybe "|" ; maybe a trailing | after the factor
  # match 3 ([^:]*) level; a  possibly empty string of non-: chracters
  # non matching group maybe ":"
  rexPat='([^|:]+)(?:[|])?([^:]*)(?:[:])?'
  theMatchList=str_match_all(theseCoefs,rexPat)
  gfacLevList=list()
  ii=0
  for (i in 1:length(theMatchList)){
    # if there is no level, it's a covariate so don't
    # accumulate it in factors
    tlev=theMatchList[[i]][ ,3 ]
    # First element of tlev will tell us if its a covariate
    if (!is.na(tlev[1]))
      ii=ii+1
    gfacLevList[[ii]]=list()
    gfacLevList[[ii]]$fac=theMatchList[[i]][,2 ]
    gfacLevList[[ii]]$lev=theMatchList[[i]][,3 ]
  }
  #-dbg     catln('Returning from getGFacAndLevNames')
  # showvars(gfacLevList)
  return(gfacLevList)
}
# getGFacAndLevNames(coefNames)
# unlist(getGFacAndLevNames(coefNames)) This will return a long character vector



## ------------------------------------------------------------------------
#' getOneFamilyDfrAndCFacStr Get a data frame and a categorical factor string of a coefficient family name
#'
#' @param thisFam -- a coefficient family name  (string) egs C:V:F1, C,V
#' @param coefNames -- character vector of coefficient names C|b:F1 etc.
#' @param gcoefFamCv -- a character vector of coefFam names
#'
#' @return a list with components dfr - data frame reflecting true factor  and cFacStr a string with the "true" multi level factors involved (covariates  removed)
#' @details
#'  dataframe  and  cFacStr are useful useful for creating model matrix to
#'  remove the marginal values of factors mm = model.matrix.
#' @export
#'
#' @examples
#' print('noexamples')
getOneFamilyDfrAndCFacStr=function(thisFam,coefNames,gcoefFamCv){
  # browser()
  thiscoefFamCFacCv=c()
  #-dbg     catln('thisFam', thisFam)
  # browser()
  inxTheseFamCoefs=which(gcoefFamCv==thisFam)
  nTheseCoefs=length(inxTheseFamCoefs)
  theseNames=coefNames [inxTheseFamCoefs]
  # Calculate the number of factors 1 + number of colons
  nTheseFacs=1+str_count(theseNames[1],':')
  theseFacsAndLevsList=getGFacAndLevNames(theseNames)
  # Example:
  # "gfacLevList": list
  # [[1]]
  # [[1]]$fac
  # [1] "C" "F1"
  #
  # [[1]]$lev
  # [1] "b" NA
  #
  # [[2]]
  # [[2]]$fac
  # [1] "C" "F1"
  #
  # [[2]]$lev
  # [1] "d" NA
  #
  # We can get the family factor names from first level
  thisFamiyFacNames=theseFacsAndLevsList[[1]]$fac
  #-dbg     catln('thisFamiyFacNames',thisFamiyFacNames)
  # # Even ones are the levels for each factor
  # Only need first coef levels to determine if it's a factor or covariate
  firstCoefLevNames=theseFacsAndLevsList[[1]]$lev
  # browser()
  # dfList will later be converted to data.frame
  dfList=list()
  iTrueFac=0
  for (ifac in 1:nTheseFacs){
    # If  level of firstCoefLevNames[ifac]  is NA then we have a
    # covariate so we don't need it in the design: this skips it.
    # it also doesn't accumulate cocariate names in    thiscoefFamCFacCv
    if (!is.na(firstCoefLevNames[ifac])){
      iTrueFac=iTrueFac+1
      thiscoefFamCFacCv[iTrueFac]=thisFamiyFacNames[ifac]
      tfacLevs=rep('',nTheseCoefs)
      for (jcoef in 1:nTheseCoefs ){
        tfacLevs[jcoef]=theseFacsAndLevsList[[jcoef]]$lev[ifac]
      }
      dfList[[ thisFamiyFacNames[ifac] ]]=factor(tfacLevs)
      #   catln('dfList in factor loop')
      # print(dfList)
      # browser()
    }
    #
  } # end for ifac

  thisFamCatFacStr=paste0(unlist(   thiscoefFamCFacCv),collapse=':')

  # make the
  # We have enough info to construct the empty data frame    as.ch
  #-dbg     catln('thisFam', thisFam)
  #-dbg     catln('dfList for family')
  #-dbg         print(str(dfList))
  # need as.data.frame not data.frame (???this won't work if ther is moe than one list element???)
  dfr=as.data.frame(dfList)
  #  catln('dfr for this family isnull', is.null(dfr))
  #    print(str(dfr))
  # browser()
  return(list(dfr=dfr, cFacStr=thisFamCatFacStr))
}


## ------------------------------------------------------------------------
makeFormStrFromCFacStr=function(cFacStr, yVarName=''){
  formStr=paste0(yVarName,'~',str_replace(cFacStr,':','*'),
                 '-',cFacStr)
}

## ------------------------------------------------------------------------
#' Get family names from coefficient names (several coefNames may have same family name)
#'
#' @param coefNames a character vector of coefficient names
#'
#' @return a character vector of family names (of same length as input)
#' @seealso getFamilyWiseCoefList()
#' @export
#'
#' @examples
#' print(getFamNamesFromCoefNames(c('a|1:b_2:x','a_1','b_2', 'x')))
getFamNamesFromCoefNames=function(coefNames){
  # " Remove any <<|>><level> parts of
  # <factor>(<|><level>)?(<<:>><factor>(<|><level>)?)* patterns
  gcoefFamCv=str_replace_all( coefNames,"\\|.*?(:|$)","\\1")
  #-dbg print(gcoefFamCv)
  # famNames=unique(gcoefFamCv)
}
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#' makeModelMatFromFamilyDfrAndCFacStr
#'
#' @param dfr  data frame as created in getOneFamilyDfrAndCFacStr
#' @param cFacStr string with true facrors (getOneFamilyDfrAndCFacStr)
#'
#' @return  model matrix appropriate for marginal-centering of the factorial coeefficient family
#' @export
#'
#' @examples
#' print('no example  makeModelMatFromFamilyDfrAndCFacStr')
makeModelMatFromFamilyDfrAndCFacStr=function(dfr,cFacStr){
  formStr=makeFormStrFromCFacStr(cFacStr)
  mm= model.matrix(as.formula(formStr),data=dfr)
  return(mm)
}


## ------------------------------------------------------------------------
#' Create a residual maker matrix from coefficient names.
#'
#' @param coefNames  coefficient names a e.g. c("C|d", "V|a",... "C_d:V_a: ...., C|b:V|a:y")
#' @param shouldReturnBigMatrix  default FALSE. If true, RM_bigMat component returned
#' @param showLMTest  default FALSE. f True an example test is run for random coefficients.
#'
#' @return  result: a list with components familyNames, RM_list and optionally RM_bigMat

#' @details
#'  Given an character vector of coefficient names for a HMNL, this function creates
#   a result consisting of
#     result$familyNames the family names of groups of coefficients. EG
#    if coefNames == c("C|d", "V|a",... "C_d:V_a: ...., C|b:V|a:y"), the returned
#'    familyNames would include "C", "V", C:V", "C:V:y"
#'    result$RM_list a list of residual maker matrices for anihilating all marginal sums of
#'   factorial coefficients.
#' and result$RM_bigMat  a large (ncoef x ncoef) matrix with non-zero blocks constisting of
#  result$RM_list[[ifam]] ifam=1:nfam.

#'
#'
#'
#'
#'
#' @export
#' @examples
#' coefNames=getTestCoefNames()
#' result=createResidualMaker(coefNames, shouldReturnBigMatrix=TRUE, showLMTest=TRUE)
#' if (requireNamespace("nutils", quietly=TRUE)) {
#'    nutils::showvars(result)
#'  }else{
#'    print(str(result))
#'  }
createResidualMaker=function (coefNames, shouldReturnBigMatrix=FALSE, showLMTest=FALSE){

  fHatMat=function(X){
    hat= X %*% solve(t(X) %*% X) %*% t(X)
  }
  ## Check via lm that the sum to zero constraints actually work

  checkVia_lm= function(thisFamName,thisFamCFacName,tmpDfr,btmpRes){
    lmForm= makeFormStrFromCFacStr(coefFamCFacCv[thisFam],"btmpRes")
    tmpdfr$btmpRes=btmpRes
    catln('grand mean resids', mean(btmpRes))
    catln(lmForm)
    print(summary(lm(as.formula(lmForm),data=thisFamDfrAndTrueFacList$dfr)))
  }
  # ----------
  ncoefAll=length(coefNames)
  #-dbg catln('gcoefFamCv')
  # Strip levels to get gust General Factor families
  # So  C|b:V|i:x  --> C:V:x
  gcoefFamCv=str_replace_all( coefNames,"\\|.*?(:|$)","\\1")
  #-dbg print(gcoefFamCv)
  uFamStr=unique(gcoefFamCv)
  # Get the coeficient names for stamping on the rows/columns of Residual Maker sub-matrices
  coefNamesPerFamList=list()
  for (fam in uFamStr){
    coefNamesPerFamList[[fam]] = coefNames[gcoefFamCv==fam]
  }
  coefFamCFacCv=c()
  RM_list=list()
  for (thisFam in uFamStr){
    thisFamDfrAndTrueFacList= getOneFamilyDfrAndCFacStr(thisFam,coefNames,gcoefFamCv)
    # browser()
    thisModelMat= makeModelMatFromFamilyDfrAndCFacStr(thisFamDfrAndTrueFacList$dfr,thisFamDfrAndTrueFacList$cFacStr)
    coefFamCFacCv[thisFam] = thisFamDfrAndTrueFacList$cFacStr
    # coefFamTrueFacList[[thisFam]]=paste0(unlist(   thiscoefFamCFacCv),collapse=':')
    # showvars(thisModelMat)
    # https://en.wikipedia.org/wiki/Projection_matrix
    # Projection operator P (hat matrix by another name)
    Hat = fHatMat(thisModelMat)
    # Residual maker  operator (I-P)
    RM_list[[thisFam]]= diag(nrow(Hat))-Hat
    # Stamp the coef names on this:
    dimnames(RM_list[[thisFam]])=list(coefNamesPerFamList[[thisFam]],coefNamesPerFamList[[thisFam]])

    # browser()
    #
    if (showLMTest){
      ncoef=nrow(Hat)
      btmp=as.matrix(20*rnorm(ncoef))
      btmpRes= RM_list[[thisFam]] %*% btmp
      tmpdfr=thisFamDfrAndTrueFacList$dfr
      checkVia_lm(thisFamName,thisFamCFacName=coefFamCFacCv[thisFam],tmpDfr,btmpRes)
      # options(contrast=saveContrOpt)
      # browser()
    }
    #  Now pack this into the grand
  } #end for thisFam

  #  Create 0 matrix to hold the RM_bigMat -- the residual operator for
  #  ALL the family wise residual operator submatries (block diagonal)
  if (shouldReturnBigMatrix){
    RM_bigMat = matrix(0,ncoefAll,ncoefAll, dimnames=list(coefNames,coefNames))
    for (thisFam in uFamStr){
      inxTheseFamCoefs=which(gcoefFamCv==thisFam)
      RM_bigMat[inxTheseFamCoefs,inxTheseFamCoefs]=RM_list[[thisFam]]
    }
  }else{
    RM_bigMat=NULL
  }

  return( list(famNames=uFamStr,RM_list=RM_list,RM_bigMat=RM_bigMat))
  # showvars(coefFamCFacCv)

}

## ----getFamilyWiseCoefList()--------------------------------------------------------------------
#' Get the familynames for each coefficient and organize into list
#'
#'
#' @return  a named list (by family names) of coefficients
#' for each family
#' @seealso getFamNamesFromCoefNames()
#' @export
#'
#' @examples
#' fwCoefNames=getFamilyWiseCoefList(c("A|a, A|b, B|a, B|b, A|a:B|b"))
#' print(str(fwCoefNames))
getFamilyWiseCoefList=function(coefNames){
  famNames=getFamNamesFromCoefNames(coefNames)
  familyWiseCoefList=list()
  for (fam in unique(famNames)){
    familyWiseCoefList[[fam]] = coefNames[famNames==fam]
  }
return(familyWiseCoefList)
}

## ------------------------------------------------------------------------


# So we need to try to run ../inst/BayesHMNL.stan on some test data using
# bigmatrix... Then we can work on separate covariances per family of bcoefficients and stacking
#
# # data {
#   int<lower=2> nRespCat; // Number of response category alternatives (choices) in each scenario
#   int<lower=1> K; // Number of columns in each design matrix list X[t]
#   int<lower=1> nP; // Number of participants (respondants/ agents/ persons/perceivers)
#   int <lower=1> T; // total number of observation trials  (grand trials) //// int<lower=1> S; // Number of scenarios per respondent
#   int<lower=1,upper=nRespCat> Y[T]; // YB[nP, S]; // best choices
#   // int<lower=1,upper=nRespCat> YW[nP, S]; // worst choices
#   matrix[nRespCat, K] X[T]; // was   matrix[nRespCat, K] X[nP, S]; // matrix of attributes for each obs
#   int<lower=1, upper=nP> PID[T]; // Added tmn. Serial identifier for participant on each observation trial
# }

