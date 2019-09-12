# sourceRFiles.R
#
# https://stackoverflow.com/questions/6979917/how-to-unload-a-package-without-restarting-r
# detach("package:vegan", unload=TRUE)
doit=readline('Want me to unload bhmnl and source files y/n? ')
if(doit=='y'){
    if("package:bhmnl" %in% search()) detach("package:bhmnl", unload=TRUE)
    flist=c(
        "/Users/tnearey/tmnRPkgs/bhmnl/R/bhmnl-package.R",
        "/Users/tnearey/tmnRPkgs/bhmnl/R/factorialResponseLogisticExpandedDesign.R",
        "/Users/tnearey/tmnRPkgs/bhmnl/R/ParseFactoredCoefNames.R"
    )
    for (i in 1:length(flist))
        source(flist[i],echo=FALSE)
}
