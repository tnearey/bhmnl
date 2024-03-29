% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factorialResponseLogisticExpandedDesign.R
\name{make_rsSumDmxFrom_srDtb}
\alias{make_rsSumDmxFrom_srDtb}
\title{Make a response-varying-first Sum-to-Zero Dmx}
\usage{
make_rsSumDmxFrom_srDtb(form = NULL, srDtb = NULL,
  stimVarColNames = NULL, respFacDescr = NULL,
  includeMasterCatInDesired = TRUE)
}
\arguments{
\item{form}{a formula or formula string with response factors and interactions of
response factors with stimulus parameters   e.g. ~ C*V+ C:F1+ V:F2
Note that stimulus only strings will be ignored
and that  : interactions between stimulus and responses will frist be
converted to  so C:F1:V:F2 would be pre converted to C*F1:V*F2
This is necessary to create correct design matrices for things like V:cond
where cond is a STIMULUS FACTOR}

\item{srDtb}{a stimulus(-response0 data table
(Note: trial-by-trial responses are required only for model fitting and
for compling respDtbs from full srDtb's .
stimulus only tables can be used in e.g. test data generation)}

\item{stimVarColNames}{names of stimulus factors}

\item{respFacDescr}{either a list of response factor names  (first of which must be the "master category")
or a respFacDtb (of which all response col names will be used )}

\item{includeMasterCatInDesired}{(defalt TRUE) include the master response category in
the rsDmx response factors}
}
\value{
a rsDmx: a response-stimulus data matrix with attributes 'masterCatLevels' 'expandedCoefNames'
    and  coefExpanderList with components
   "masterCatLevels" , "expandedCoefNames", "expanderMatrix, expandedRowFamNames"
   "sumColFamNames",  "expandedCoefFamIndicesList", "sumColRowCoefIndicesList"
   See Details.
}
\description{
Make a design matrix with attributes masterCatLevels among others
}
\details{
The output rsDmx matrix has rows in groups nCat (= attr(rsDmx, 'masterCatLevels))
 corresponding to the  number of categories in  the input column identified in
 masterCatColName. It also has the following attributes:
 \itemize{
  \item expandedCoefNames corresponds to the column names or output rsDmx;
  \item expandedCoefNames corresponds to the column names or output rsDmx;
  \item expanderMatrix to a matrix such that:
  \item \code{expandedCoefs = coefs \%*\% t(expanderMatrix)}    applied to the coefficients
   of the columnns of the sum to zero coded rsDmx returns the expanded coefficient matrix
  (with the last reference category of each factor appropriately computed)
  \item expandedCoefFamIndicesList - named (by family of coefficients) list of integer indices into
 the rows of the expanderMatrix. Can be use to access gropus of the the  expandedCoefs above
}
}
\examples{
n=50
srDtb= data.table(resp=factor(sample(c("A","B","C"),n,replace=TRUE) 1),
x=nrand(n),y=nrand(n)))
stimVarColNames=c('x','y')
respFacDescr=list("resp")
make_rsSumDmxFrom_srDtb(~resp*x+resp*y, srDtb, stimVarColnames, respFacDescr)
     $expandedRowFamNames = getFamNamesFromCoefNames(rownames(expanderMatrix))
     $sumColFamNames=getFamNamesFromCoefNames(colnames(expanderMatrix))
     # for famName in unique(expandedRowFamNames)

     $expandedCoefFamIndicesList[[famName]]=which(expandedRowFamNames==famName)
     $sumColRowCoefIndicesList[[famName]]=which(expandedRowFamNames==colname)

}
