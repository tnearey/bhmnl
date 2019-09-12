## ------------------------------------------------------------------------
# library(bhmnl)
# ' @importFrom stats as.formula model.matrix
# ' @importFrom utils  head
#'
## ---catln()---------------------------------------
#' Call `cat()` with trailing newline
#'
#' @param ... args to pass to cat()
#'
#' @export
#'
#' @examples
#' catln('Boo',27,'hoo')
catln <- function(...) {
  cat(...,'\n')
}
# ---normvec()--------------------------------------------------------

#' Normalize  elements of(hopefully all positive ) vector (max forced to 1.0)
#'
#' @param x vector of non-negatives
#'
#' @return normalized vector (sums to 1s)
#' @export
#'
#' @examples
#' normvec(c(1,2,3))
normvec=function(x){
  if (any(x<0)) stop(' x must be all non-negative')
  y=x/max(x)
  y=y/sum(y)
}
## from http://tr.im/hH5A thru https://gist.github.com/aufrank/83572
# ---logsumexp()--------------------------------------------------------

#' Log of sum of exp of vector (numerically stable.)
#'
#' @param x vector
#'
#' @return log(sum(exp(x))) but better
#' @export
#'
#' @examples
#'  logsumexp(c(1,2,3))
#
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}
# ---softmax()--------------------------------------------------------
#' Softmax function exp(x)/logsumexp(x)
#'
#' @param x vector (of eg log likelihoods)
#'
#' @return exp(x)/logsumexp(x)
#' @export
#'
#' @examples
#' softmax(log(c(1,2,3)))
softmax <- function (x) {
  exp(x - logsumexp(x))
}

## --showRespFacDtbLevs-----------------------------------------------------------
## Debugging
showRespFacDtbLevs=function(respFacDtb){
  for (n in names(respFacDtb)){
    catln(n, levels(respFacDtb[[n]]))
  }
}


## --showColInfo()----------------------------------------------------------------------
showColInfo=function(xdtb){
  catln('class of input', class(xdtb))
  for (n in names(xdtb)){
    t=xdtb[[n]]
    catln(n, class(t))
    nna=sum(is.na(t))
    nuniq=length(unique(t))
    if (is.factor(t)){
      catln(levels(t))
    }else if (is.numeric(t)){
      catln('min', min(t), 'max', max(t), 'nuniq', nuniq)
    }else{
      catln('nuniq', nuniq)
    }
    if (nna>0) {catln('*** n na ',nna, '***')}
  }
}


## -----checkAllAreFactors()-------------------------------------------------------------------
checkAllAreFactors=function(facNames,dtb){
  checkit=TRUE
  for (f in facNames){
    if (!is.factor(dtb[[f]])){
      catln('Not factor: ', f)
      checkit=FALSE
    }
  }
  if (!checkit) stop('Must all be factors')
  return(checkit)
}



## --make_respFacDtb()----------------------------------------------------------------------
#'  Make a response factor data table from an srDtb data.table
#'
#' @param srDtb  -- a stimulus-response data.table (one trial per row)
#'  or a data table compatible with respFacDtb except for masterCatLevels attributes
#' @param masterCatColName -- column name for master category
#' @param desiredFactorColNames -- column names to be used in formula
#' @param includeMasterCatInDesired -- boolean (default TRUE) whether to include masterCatColName among desiredFactorColNames

#' @details  Collect response factors from an srDtb and collect into a compact respFacDtb
#' Note if you want to exclude  masterCatColName from output of columns,
#' you must set includeMasterCatInDesired=FALSE
# .
#' @return - returns a respFacDtb data.table with  masterCatLevel attribute
#'    (and class == data.frame;data.table; respFacDtb). Note only has nlevels(masterCat) rows
#'    and columns are masterCat and any other factors that decompose mastercat
#'
#' @export
#'

make_respFacDtb=function(srDtb,masterCatColName,desiredFactorColNames, includeMasterCatInDesired=TRUE){

  facnames=unique(c(masterCatColName,desiredFactorColNames))
  if (!checkAllAreFactors(facnames,srDtb)) stop(' Programming error')
  masterCatLevels=levels(srDtb[[masterCatColName]]);
  inxFirstAppearanceMasterCat=c()
  if (includeMasterCatInDesired && !(masterCatColName %in% desiredFactorColNames )){
    desiredFactorColNames=c(masterCatColName,desiredFactorColNames)
  }
  masterCatLevs=levels(srDtb[[masterCatColName]])
  for (i in seq_along(masterCatLevs)){
    s=masterCatLevs[i];
    # Get first occurrence of each syl, cn and vow corresponding to each level of syl
    inxFirstAppearanceMasterCat[i]=which(srDtb[[masterCatColName]]==s)[1]
  }
  # dotdot https://stackoverflow.com/a/41052245/1795127
  # dotdot didnt work here use with=FALSE
  respFacDtb=srDtb[inxFirstAppearanceMasterCat, desiredFactorColNames, with=FALSE]
  # Make sure masterCat is legitimate (must be max number of levels)
  checkLevs=sapply(respFacDtb, function(x) length(levels(x)))
  if  (length(masterCatLevels)< max(checkLevs)){
    stop('MasterCat must have at least as many levels as any other response factor')
  }
  attr(respFacDtb,'masterCatLevels')= masterCatLevels;
  class(respFacDtb) = c(class(respFacDtb),'respFacDtb')
  return(respFacDtb)
}


## --padLevels()----------------------------------------------------------------------
padLevels=function(x){
  if ( '_' != levels(x)[1])
    y=factor(x,levels=c('_',levels(x)))
  else
    y=x;
}


## --padLevels_facRespDtb-------------------------------------------------------------------
padLevels_facRespDtb= function(respFacDtb){
  # Note may does NOT work "in place" as indexing via []
  # Expand the levels with a dummy blank one at the front
  # This is to trick model.matrix into specifying a full matrix
  # This appears to work:
  # Get a simple name for car::contr.Treatment and set the decorator
  # This could all be done inside data.table by referene world,
  # using .SD, .SDcols and .L apply,
  # I think we are triggering a copy and maybe a conversion to data.frame
  for (fac in names(respFacDtb)){
    respFacDtb[[fac]]=padLevels(respFacDtb[[fac]])
    # catln(fac, 'still data .table', is.data.table(respFacDtb))
  }
  # catln(' before return padLevels_facRespDtb')
  # showRespFacDtbLevs(respFacDtb)
  return(as.data.table(respFacDtb))
  if (!("respFacDtb" %in%  class(respFacDtb)) ) {
    message('Stamping respFacDtb  class')
    class(respFacDtb)= c(class(respFacDtb),'respFacDtb')
  }
}


## ---padLevels_facStim_srDtb()---------------------------------------------------------------------
#' Pad the levels of all factors of an srDtb data table.
#'
#' @param srDtb  a (trial-wise) data table with stimulus factors that
#'  may need "_" level padding
#' @param stimVarColNames names of columns including stimuli
#'
#' @return the srDtb with padded stimulus factor levels
#' @export
#'
#' @examples
#' print("no example")
padLevels_facStim_srDtb= function(srDtb,stimVarColNames){
  # Expand any factorial stimulus factors
  for (fac in stimVarColNames){
    # browser()
    # Need to address srdtb[[fac]] cause we don't want to return data.table
    if (is.factor(srDtb[[fac]])){
      # double dot here before := only because we want to extract the vector on rhs
      # Double dot construction ..fac := not working
      # srDtb[,..fac := padLevels(srDtb[[fac]])]
      # But parentesis does work in assignment
      srDtb[,(fac) := padLevels(srDtb[[fac]])]
      # srDtb[,fac := padLevels(srDtb[[fac]]), with=FALSE]
      # catln(fac, levels(srDtb[[fac]]))
      # browser()
    }

    #
    # if (is.factor (srDtb[[fac]])){
    # Version 1.12.2 should work with character vector and ..
    # srDtb[[fac]]=padLevels(srDtb[[fac]])
    # catln(fac, 'still data .table', is.data.table(srDtb))
    # catln('levs stimfac', fac)
    # print(levels(srDtb[[fac]]))
    # }
  }
  # catln(' before return padLevels_facRespDtb')

  return(as.data.table(srDtb))
}

## --make_stimRespFacDtb()----------------------------------------------------------------------
#' Make a stimulus-response-factor data.table from respFacDeb and srDtb
#'
#' @param respFacDtb response factor Data table
#' @param srDtb stimulus response data table
#' @param stimVarColNames names for stimulus variables in srDtb to be retained
#'
#' @return returns a stimRespFacDtb
#' @export
#'
#' @examples
#' cat('no examples \n')
make_stimRespFacDtb= function(respFacDtb,srDtb,stimVarColNames){
  # This is slowish but it's pretty small . Could be done with .SD .SDcols and lapply, but this is clearer.

  ### Set up the car:Treatment contrasts for clearer labeling. Then build the srDmx design matrix
  #  CBIND the columns of expanded  srDtb and respFacDtb factors, and make it a model.matrix
  # Then extract just the design part
  # paste the columns of row expanded respCat factors and row expanded stimulus properties
  # Data table is useful because we can copy its rows efficiently.

  nObs=nrow(srDtb)
  # Get this before expansion of factors
  nrespCat=nrow(respFacDtb)
  ## Need to keep track of stimulus properties extracted from srDtb
  inxSrDtbCopy=rep(1:nObs,times=nrespCat)
  inxRespFacCopy=rep(1:nrespCat,each=nObs)
  #  'dot dot' convention  https://stackoverflow.com/a/41052245/1795127
  stimRespFacStimDtb = cbind(respFacDtb[inxRespFacCopy, ],srDtb[inxSrDtbCopy, ..(stimVarColNames)])
}



## make_respFacStimDtb()------------------------------------------------------------------------
make_respFacStimDtb= function(respFacDtb,srDtb,stimVarColNames){
  # This is slowish but it's pretty small . Could be done with .SD .SDcols and lapply, but this is clearer.

  ### Set up the car:Treatment contrasts for clearer labeling. Then build the srDmx design matrix
  #  CBIND the columns of expanded  srDtb and respFacDtb factors, and make it a model.matrix
  # Then extract just the design part
  # paste the columns of row expanded respCat factors and row expanded stimulus properties
  # Data table is useful because we can copy its rows efficiently.

  nObs=nrow(srDtb)
  # Get this before expansion of factors
  nrespCat=nrow(respFacDtb)
  ## Need to keep track of stimulus properties extracted from srDtb
  inxSrDtbCopy=rep(1:nObs,each=nrespCat)
  inxRespFacCopy=rep(1:nrespCat,times=nObs)

  # catln('respFacDtbCopy dim', dim(respFacDtb[inxRespFacCopy, ]))
  # catln('srDtb dim', dim(srDtb[inxSrDtbCopy, ]))
  # Use dotdot notation to refer to column vars by char vec [[ not working]]
  # catln('Browser before ..')
  respFacStimDtb = cbind(respFacDtb[inxRespFacCopy, ],
                         srDtb[inxSrDtbCopy,stimVarColNames,with=FALSE])
}


## ---make_indicatorDmx()---------------------------------------------------------------------
#' make_indicatorDmx Make a Design Matrix (sr- or rs-)
#'
#' @param form a formula or formula string
#' @param respFacStimDtbPadded a padded factorial StimDtb
#' @param respVarColNames column names for response variavles
#' @param stimVarColNames column names for stimulus variables
#'
#' @return returns an rsDmx or an srDmx
#' @export
#'
#' @examples
#' cat('no examples \n')
make_indicatorDmx=function(form, respFacStimDtbPadded,respVarColNames,stimVarColNames){
  # Define contr.treat and decorator here.
  # Make the model matrix don't forget the 1.
  # Add also bare stimulus variables involved in formula
  # The intercept and stimvars will get stripped at end.
  # Note if respFacStimDtbPadded then you get an rsDnx with response cat
  # varying fastest down rows, if stimRespFacDtbPadded, then you get
  # an srDmx with stinulus varying fastest and response cats stacked.
  if (class(form)=='formula'){
    # https://stackoverflow.com/a/14671346/1795127
    # Eliminate extra spaces too
    formStr= gsub('  +',' ',Reduce(paste, deparse(form)))
  }else if (is.character(form)){
    formStr=form
  }else{
    stop(' Input `form` should be  of class character string or formula')
  }
  # Get rid of anything up to ~ in formStr
  formStr=sub('.*\\~','',formStr)
  varsInForm=all.vars(as.formula(paste0('~',formStr)))
  # Find the stim varsvariables that are included in the formula
  stimVarsInForm = intersect(stimVarColNames, varsInForm)
  respVarsInForm = intersect(respVarColNames,varsInForm)
  # catln('stimVarsInForm',stimVarsInForm)
  # catln('respVarsInForm',respVarsInForm)
  # Also get factored stimuli
  # factoredStimInForm=stimvarsInForm(sapply(is.factor,)
  # sapply(is.factor )
  #
  # Get names of stim factors in form to add then in to mmCotrastList
  binxStimFactors=sapply(stimVarsInForm, function(x) is.factor(respFacStimDtbPadded[[x]]))
  # browser()
  stimFactorsInForm=stimVarsInForm[binxStimFactors]

  mmContrastList=list()
  for (sv in c(respVarsInForm,stimFactorsInForm)){
    # Compile the contrast matrix
    mmContrastList[[sv]]="myIndicatorContrFun"
  }
  # catln('mmContrastList')
  # print(mmContrastList)

  # ???We have to add any stim vars to near the begining of the formula to keep from triggering "dummy" coding
  # in model.matrix. This includes factored stimulus variables (I think)

  nstimVarsAug= length(stimVarsInForm)
  # paste0( "1+ ", paste0(f,collapse=' + '), ' + ', 'a + b + c')
  # Paste the stimulus variables at the END?
  augFormStr=paste0(" ~ 1 + ",formStr  ,
                    " + " ,paste0(stimVarsInForm, collapse = " + ") )
  # catln('formStr' ,formStr)
  # catln('augFormStr',augFormStr )
  mmx=model.matrix(as.formula(augFormStr) ,respFacStimDtbPadded, contrasts.arg = mmContrastList)

  # catln('mmx')
  # summary(mmx)
  # catln(' COlnames mmxRef')
  #
  mmxNames=colnames(mmx)
  # Remove the columns  exra columns 1 + number of stimvars in formula
  # Good terms hould have contrast "|l<level> in in them either at begining or after :
  # binxKeep=grepl('(^|:).*?\\|',mmxNames)
  # Above doesn't work.
  #  We need to delete any column that doesn't have a stimulus variable
  #  This will do it.
  effectsList= lapply(strsplit(mmxNames,':'), function(x) gsub('\\|.*$','',x))
  binxKeep=sapply(effectsList, function(x) any(x %in% respVarsInForm))
  # catln('broswer')
  # browser()
  # catln('passing mmx names', mmxNames[binxKeep])
  srDmx=as.matrix(mmx)[,binxKeep]
  return(srDmx)
}


## ---make_sumDmx()---------------------------------------------------------------------
#'  Make a sum to zero Design Matrix (sr- or rs-)
#'
#' @param form a formula or formula string
#' @param respFacStimDtbPadded a padded factorial StimDtb
#' @param respVarColNames column names for response variavles
#' @param stimVarColNames column names for stimulus variables
#'
#' @return returns an rsDmx or an srDmx
#' @export
#'
#' @examples
#' cat('no examples \n')
make_sumDmx=function(form, respFacStimDtb,respVarColNames,stimVarColNames){
  # Define contr.treat and decorator here.
  # Make the model matrix don't forget the 1.
  # Add also bare stimulus variables involved in formula
  # The intercept and stimvars will get stripped at end.
  # Note if respFacStimDtb then you get an rsDnx with response cat
  # varying fastest down rows, if stimRespFacDtbPadded, then you get
  # an srDmx with stinulus varying fastest and response cats stacked.
  if (class(form)=='formula'){
    # https://stackoverflow.com/a/14671346/1795127
    # Eliminate extra spaces too
    formStr= gsub('  +',' ',Reduce(paste, deparse(form)))
  }else if (is.character(form)){
    formStr=form
  }else{
    stop(' Input `form` should be  of class character string or formula')
  }
  # Get rid of anything up to ~ in formStr
  formStr=sub('.*\\~','',formStr)
  varsInForm=all.vars(as.formula(paste0('~',formStr)))
  # Find the stim varsvariables that are included in the formula
  stimVarsInForm = intersect(stimVarColNames, varsInForm)
  respVarsInForm = intersect(respVarColNames,varsInForm)
  # catln('stimVarsInForm',stimVarsInForm)
  # catln('respVarsInForm',respVarsInForm)
  # Also get factored stimuli
  # factoredStimInForm=stimvarsInForm(sapply(is.factor,)
  # sapply(is.factor )
  #
  # Get names of stim factors in form to add then in to mmCotrastList
  binxStimFactors=sapply(stimVarsInForm, function(x) is.factor(respFacStimDtb[[x]]))
  # browser()
  stimFactorsInForm=stimVarsInForm[binxStimFactors]

  mmContrastList=list()
  for (sv in c(respVarsInForm,stimFactorsInForm)){
    # Compile the contrast matrix
    mmContrastList[[sv]]="mySumContrFun"
  }
  #-dbg   catln('mmContrastList')
  #-dbg   print(mmContrastList)

  # ???We have to add any stim vars to near the begining of the formula to keep from triggering "dummy" coding
  # in model.matrix. This includes factored stimulus variables (I think)

  nstimVarsAug= length(stimVarsInForm)
  # paste0( "1+ ", paste0(f,collapse=' + '), ' + ', 'a + b + c')
  # Paste the stimulus variables at the END?
  augFormStr=paste0(" ~ 1 + ",formStr  ,
                    " + " ,paste0(stimVarsInForm, collapse = " + ") )
  # catln('formStr' ,formStr)
  # catln('augFormStr',augFormStr )
  mmx=model.matrix(as.formula(augFormStr) ,respFacStimDtb, contrasts.arg = mmContrastList)

  # catln('mmx')
  # summary(mmx)
  # catln(' COlnames mmxRef')
  #
  mmxNames=colnames(mmx)
  # Remove the columns  exra columns 1 + number of stimvars in formula
  # Good terms hould have contrast "|l<level> in in them either at begining or after :
  # binxKeep=grepl('(^|:).*?\\|',mmxNames)
  # Above doesn't work.
  #  We need to delete any column that doesn't have a stimulus variable
  #  This will do it.
  effectsList= lapply(strsplit(mmxNames,':'), function(x) gsub('\\|.*$','',x))
  binxKeep=sapply(effectsList, function(x) any(x %in% respVarsInForm))
  # catln('broswer')
  # browser()
  # catln('passing mmx names', mmxNames[binxKeep])
  srDmx=as.matrix(mmx)[,binxKeep]
  return(srDmx)
}


## --myIndicatorContrFun()----------------------------------------------------------------------
#'   A decorated contrast function always forcing INDICATOR levels
#'
#' @param levLabs  charvec of level labels
#' @param base base category
#' @param contrasts  contrast=TRUE (ignored)
#'
#' @return  contr returns contrast matrix
#' @export
#'
#' @examples
#' cat('no examples\n')
myIndicatorContrFun = function (levLabs, base = 1, contrasts = TRUE){
  # Specialized function for labeling contrsts witn <facname>|<levelLab>
  # I.e. pipe symbol as prefix for level,
  # inspired by car:contr.Treatment
  #-dbg   catln('class (levLabs) in myIndicatorContrFun', class(levLabs))
  #-dbg   catln('base', 'contrasts', base, contrasts)
  # readline('myIndicatorContrFun :')
  if (!is.character(levLabs)) {
    stop('levLabs (first arg) needs to be a character vector.')
  }
  n <- length(levLabs)

  if (!(base==1 && contrasts)){
    stop( 'myIndicatorContrFun must have base=1 and contrasts=TRUE')
  }
  contrastNames <- paste0('|', levLabs)
  contr <- array(0, c(n, n), list(levLabs, contrastNames))
  diag(contr) <- 1
  if (contrasts) {
    if (n < 2)
      stop("Contrast requires at least 2 d.f.")
    if (base < 1 | base > n)
      stop("Base level must be >= 1 and <= n")
    contr <- contr[, -base, drop = FALSE]
  }
  return(contr)
}
## --mySumContrFun()----------------------------------------------------------------------
#'   A decorated contrast function always forcing sum-to-zero  levels (last cat is ref)
#'
#' @param levLabs  charvec of level labels
#' @param base ( not used should be null)
#' @param contrasts  contrast=TRUE (ignored)
#'
#' @return  contr returns contrast matrix
#' @export
#'
#' @examples
#' cat('no examples\n')
mySumContrFun = function (levLabs, base = NULL, contrasts = TRUE){
  # Specialized function for labeling contrsts witn <facname>|<levelLab>
  # I.e. pipe symbol as prefix for level,
  # inspired by car:contr.Treatment
  #-dbg   catln('class (levLabs) in mySumContrFun', class(levLabs))
  #-dbg   catln('base', 'contrasts', base, contrasts)
  # readline('mySumContrFun :')
  if (!is.character(levLabs)) {
    stop('levLabs (first arg) needs to be a character vector.')
  }
  n <- length(levLabs)

  if (! (is.null(base) && contrasts)){
    stop( 'mySumContrFun must have base=1 and contrasts=TRUE')
  }
  contrastNames <- paste0('|', levLabs)
  contr <- array(0, c(n, n), list(levLabs, contrastNames))
  diag(contr) <- 1
  if (contrasts) {
    if (n < 2)
      stop("Contrast requires at least 2 d.f.")
    contr[n,]=-1
    contr <- contr[, -n, drop = FALSE]
  }
  return(contr)
}

# ====== make_rsIndicatorDmxFrom_srDtb() ==============
#' make_rsIndicatorDmxFrom_srDtb
#'
#' @param form  a formula or formula string with response factors and interactions of
#'   response factors with stimulus parameters   e.g. ~ C*V+ C:F1+ V:F2

#' @param srDtb  a stimulus(-response0 data table
#' (Note: trial-by-trial responses are required only for model fitting and
#' for compling respDtbs from full srDtb's .
#' stimulus only tables can be used in e.g. test data generation)
#' @param stimVarColNames  names of stimulus factors
#' @param respFacDescr  either a list of response factor names  (first of which must be the "master category")
#'       or a respFacDtb (of which all response col names will be used )
#' @param includeMasterCatInDesired  (defalt TRUE) include the master response category in
#'  the rsDmx response factors
#'
#' @return a rsDmx: a response-stimulus data matrix with attributes 'masterCatLevels'
#'
#' @export
#'
#' @examples
#' cat('no examples \n')
make_rsIndicatorDmxFrom_srDtb=function(form=NULL, srDtb=NULL,stimVarColNames=NULL,
                                       respFacDescr=NULL,
                                       includeMasterCatInDesired=TRUE){
  # respFacDescr - can be an respFacDtb or a character vector of response factor names in srDtb
  #   if the latter, then it must have the master category first (e.g. is syl )
  # srDtb can be an original srDtb or just corresponding stimulus colunns
  # browser()

  # Convert respFacDescr to respFacDtb type if necessary
  if (is.character(respFacDescr)){
    respVarColNames=respFacDescr;
    respFacDescr=make_respFacDtb(srDtb,masterCatColName = respFacDescr[1],
                                 desiredFactorColNames =respFacDescr[-1],
                                 includeMasterCatInDesired = includeMasterCatInDesired )
  }
  if (!('respFacDtb' %in% class(respFacDescr))){
    message("respFacDescr must be char vec or a respFacDtb")
    browser()
  }
  #  by now respFacDescr is a respFacDtb, so we can grab all its column names
  respVarColNames=names(respFacDescr)

  masterCatLevels=attr(respFacDescr,'masterCatLevels')

  ## Padd out all factors to indicators
  respFacDtbPadded=padLevels_facRespDtb(respFacDescr)

  #-dbg   catln('after respFacDtb padding')
  #-dbg  showRespFacDtbLevs(respFacDtbPadded)
  # Pad any factorial stimulus levels
  srDtbPadded= padLevels_facStim_srDtb(srDtb, stimVarColNames)
  #-dbg  catln('srDtbPadded')
  #-dbg  showColInfo(srDtbPadded)
  respFacStimDtbPadded= make_respFacStimDtb(respFacDtbPadded,srDtbPadded,stimVarColNames)
  # stimRespFacDtbPadded= make_stimRespFacDtb(respFacDtbPadded,srDtbPadded,stimVarColNames)
  #-dbg  catln('respFacStimDtbPadded')
  #-dbg  showColInfo(respFacStimDtbPadded)
  # showColInfo(stimRespFacDtbPadded)
  # Costruct a stimulus-response design matrix
  #-dbg  catln('respVarColNames',respVarColNames)
  #-dbg  catln('stimVarColNames',stimVarColNames)
  rsDmx= make_indicatorDmx(form, respFacStimDtbPadded,respVarColNames,stimVarColNames)
  # catln('rsDmx')
  # print(head(rsDmx,n=20))
  attr(rsDmx,'masterCatLevels')=masterCatLevels
  # browser()
  return(rsDmx)
}
# ---- make_catRespMMFormStr-----
#' make a "safe" categorial response model.matrix formula string, forcing all : to *
#'
#' @param formSpec  a formula or formula
#' @param respVarColNames charvec of response variable column names (all should be factors)
#' @param stimVarColNames charvec of response variable column names (factors or covariates)
#'
#' @return returns a character string
#' @export
#'
#' @details
#' Currently this is not used.
#' @examples
#' respVarColnames=c('syl','C','V')
#' stimVarColNames=c('F1','F2','cond')
#'formStr=make_catRespMMFormStr(formSpec="syl~C*V+C:F1+C:F2+C:cond")
#' print(formStr)
#'
make_catRespMMFormStr=function(formSpec,respVarColNames,stimVarColNames){
  # Prunes off any  LHS and converts all ':' to '*' to ensure
  # proper handling of factor interactions
  # users of this may have to remove any stimulus only columns after
  formStr=paste0(as.character(formSpec),collapse='')
  # Prune off lhs (not sure why \\~ is necessary)
  formStr=sub('.*\\~','~',formStr)
  # Just convert all : to *
  # This ensures all lower order interactions are there
  # Model matrix will sort them out later if there are redundancies
  formStr=gsub(':', '*', formStr,fixed = TRUE)
  # print(formStr)
  # browser()
  return(formStr)
}

# ---- make_rsSumDmxFrom_srDtb -----
#' Make a response-varying-first Sum-to-Zero Dmx
#'
#'  Make a design matrix with attributes masterCatLevels among others
#'
#' @param form  a formula or formula string with response factors and interactions of
#'   response factors with stimulus parameters   e.g. ~ C*V+ C:F1+ V:F2
#'   Note that stimulus only strings will be ignored
#'   and that  : interactions between stimulus and responses will frist be
#'   converted to  so C:F1:V:F2 would be pre converted to C*F1:V*F2
#'   This is necessary to create correct design matrices for things like V:cond
#'   where cond is a STIMULUS FACTOR
#'

#' @param srDtb  a stimulus(-response0 data table
#' (Note: trial-by-trial responses are required only for model fitting and
#' for compling respDtbs from full srDtb's .
#' stimulus only tables can be used in e.g. test data generation)
#' @param stimVarColNames  names of stimulus factors
#' @param respFacDescr  either a list of response factor names  (first of which must be the "master category")
#'       or a respFacDtb (of which all response col names will be used )
#' @param includeMasterCatInDesired  (defalt TRUE) include the master response category in
#'  the rsDmx response factors
#'
#' @return a rsDmx: a response-stimulus data matrix with attributes 'masterCatLevels' 'expandedCoefNames'
#'     and  coefExpanderList with components
#'    "masterCatLevels" , "expandedCoefNames", "expanderMatrix, expandedRowFamNames"
#'    "sumColFamNames",  "expandedCoefFamIndicesList", "sumColRowCoefIndicesList"
#'    See Details.
#'
#' @details
#'  The output rsDmx matrix has rows in groups nCat (= attr(rsDmx, 'masterCatLevels))
#'  corresponding to the  number of categories in  the input column identified in
#'  masterCatColName. It also has the following attributes:
#'  \itemize{
#'   \item expandedCoefNames corresponds to the column names or output rsDmx;
#'   \item expandedCoefNames corresponds to the column names or output rsDmx;
#'   \item expanderMatrix to a matrix such that:
#'   \item \code{expandedCoefs = coefs \%*\% t(expanderMatrix)}    applied to the coefficients
#'    of the columnns of the sum to zero coded rsDmx returns the expanded coefficient matrix
#'   (with the last reference category of each factor appropriately computed)
#'   \item expandedCoefFamIndicesList - named (by family of coefficients) list of integer indices into
#'  the rows of the expanderMatrix. Can be use to access gropus of the the  expandedCoefs above
#' }
#
#' @export
#'
#' @examples
#' n=50
#' srDtb= data.table(resp=factor(sample(c("A","B","C"),n,replace=TRUE) 1),
#' x=nrand(n),y=nrand(n)))
#' stimVarColNames=c('x','y')
#' respFacDescr=list("resp")
#' make_rsSumDmxFrom_srDtb(~resp*x+resp*y, srDtb, stimVarColnames, respFacDescr)
make_rsSumDmxFrom_srDtb=function(form=NULL, srDtb=NULL,stimVarColNames=NULL,
                                 respFacDescr=NULL,
                                 includeMasterCatInDesired=TRUE){
  # Convert respFacDescr to respFacDtb type if necessary
  if (is.character(respFacDescr)){
    respVarColNames=respFacDescr;
    respFacDescr=make_respFacDtb(srDtb,masterCatColName = respFacDescr[1],
                                 desiredFactorColNames =respFacDescr[-1],
                                 includeMasterCatInDesired = includeMasterCatInDesired )
  }
  if (!('respFacDtb' %in% class(respFacDescr))){
    message("respFacDescr must be char vec or a respFacDtb")
    browser()}
  #  by now respFacDescr is a respFacDtb, so we can grab all its column names
  respVarColNames=names(respFacDescr)

  masterCatLevels=attr(respFacDescr,'masterCatLevels')
  ## ONE ROW indicator processing to get column names for full indicstor dmx
  ## **** Padding to create the row names of a full-indicator data.matrix.**** ----
  ## We're really only going after data.matrix created column names. Sholdn't need
  ## ore than one row... as long as factor levels have been padded.
  respFacDtbPadded1Row=padLevels_facRespDtb(respFacDescr[1,])
  #-dbg   catln('after respFacDtbPadded1Row padding')
  #-dbg  showRespFacDtbLevs(respFacDtbPadded1Row)
  # Pad any factorial stimulus levels
  # Should only have to pad first row of srDtb
  srDtbPadded1Row= padLevels_facStim_srDtb(srDtb[1,], stimVarColNames)
  # catln('srDtbPadded1Row')
  # showColInfo(srDtbPadded1Row)
  respFacStimDtbPadded1Row= make_respFacStimDtb(respFacDtbPadded1Row,srDtbPadded1Row,stimVarColNames)
  # stimRespFacDtbPadded= make_stimRespFacDtb(respFacDtbPadded,srDtbPadded,stimVarColNames)
  # catln('respFacStimDtbPadded1Row')
  # showColInfo(respFacStimDtbPadded1Row)
  # showColInfo(stimRespFacDtbPadded)
  # Costruct a stimulus-response design matrix
  # catln('respVarColNames',respVarColNames)
  # catln('stimVarColNames',stimVarColNames)
  rsIndicatorDmx1Row= make_indicatorDmx(form, respFacStimDtbPadded1Row,respVarColNames,stimVarColNames)
  # catln('rsIndicatorDmx1Row')
  # print(head(rsDmx,n=20))
  # We access this degenerate (vectorish) structure with "names" not col names
  expandedCoefNames=names(rsIndicatorDmx1Row)
  attr(rsIndicatorDmx1Row,'masterCatLevels')=masterCatLevels
  ## ****** NEXT create the sum-to-zero version -- using complete unpaded srDtb and respFacDtb****
  respFacStimDtb= make_respFacStimDtb(respFacDescr,srDtb,stimVarColNames)
  rsSumDmx= make_sumDmx(form, respFacStimDtb,respVarColNames,stimVarColNames)

  # SO NOW WE NEED TO BUILD THE coefficient Expander MATRIX  Should be pseudo inverse of
  attr(rsSumDmx,'masterCatLevels')=masterCatLevels
  attr(rsSumDmx,'expandedCoefNames')=expandedCoefNames

  expanderMatrix=createCoefficientExpanderList(expandedCoefNames)$XM_bigMat
  expandedFamNames=getFamNamesFromCoefNames(expandedCoefNames)
  sumCoefFamBNames=getFamNamesFromCoefNames(colnames(rsSumDmx))
  sumCoefNames=colnames(rsSumDmx)
  # $ Assign the correct  column names
  colnames(expanderMatrix)= sumCoefNames
  assertthat::assert_that(all(expandedCoefNames==rownames(expanderMatrix)),
                          msg="expandedCoefNames not all equal rownames")
  #      $expandedRowIndicesList[[famName]]

  # print(str(rsSumDmx))
  # browser()
  # coefExpanderList with components
  #      $expanderMatrix [nIndicator x nSumToZero] with named rows and columns
  #'      $expandedRowFamNames = getFamNamesFromCoefNames(rownames(expanderMatrix))
  #'      $sumColFamNames=getFamNamesFromCoefNames(colnames(expanderMatrix))
  #'      # for famName in unique(expandedRowFamNames)
  #'
  #'      $expandedCoefFamIndicesList[[famName]]=which(expandedRowFamNames==famName)
  #'      $sumColRowCoefIndicesList[[famName]]=which(expandedRowFamNames==colname)
  #'

  #      expanderMatrix [nIndicator x nSumToZero] with named rows and columns
  expandedRowFamNames = getFamNamesFromCoefNames(rownames(expanderMatrix))
  sumColFamNames=getFamNamesFromCoefNames(colnames(expanderMatrix))
  expandedCoefFamIndicesList=list()
  sumColRowCoefIndicesList=list()
  for (famName in unique(expandedRowFamNames)){
    expandedCoefFamIndicesList[[famName]]=which(expandedRowFamNames==famName)
    sumColRowCoefIndicesList[[famName]]=which(expandedRowFamNames==famName)
  }
  attr(rsSumDmx, 'expanderMatrix') <- expanderMatrix
  attr(rsSumDmx, "expandedRowFamNames") <-expandedRowFamNames
  attr(rsSumDmx, "sumColFamNames") <-sumColFamNames
  attr(rsSumDmx, "expandedCoefFamIndicesList") <-expandedCoefFamIndicesList
  attr(rsSumDmx, "sumColRowCoefIndicesList") <-sumColRowCoefIndicesList

# print(str(coefExpanderList))


return(rsSumDmx)

}
