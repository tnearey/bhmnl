% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factorialResponseLogisticExpandedDesign.R
\name{make_respFacDtb}
\alias{make_respFacDtb}
\title{Make a response factor data table from an srDtb data.table}
\usage{
make_respFacDtb(srDtb, masterCatColName, desiredFactorColNames,
  includeMasterCatInDesired = TRUE)
}
\arguments{
\item{srDtb}{-- a stimulus-response data.table (one trial per row)
or a data table compatible with respFacDtb except for masterCatLevels attributes}

\item{masterCatColName}{-- column name for master category}

\item{desiredFactorColNames}{-- column names to be used in formula}

\item{includeMasterCatInDesired}{-- boolean (default TRUE) whether to include masterCatColName among desiredFactorColNames}
}
\value{
- returns a respFacDtb data.table with  masterCatLevel attribute
   (and class == data.frame;data.table; respFacDtb). Note only has nlevels(masterCat) rows
   and columns are masterCat and any other factors that decompose mastercat
}
\description{
Make a response factor data table from an srDtb data.table
}
\details{
Collect response factors from an srDtb and collect into a compact respFacDtb
Note if you want to exclude  masterCatColName from output of columns,
you must set includeMasterCatInDesired=FALSE
}
