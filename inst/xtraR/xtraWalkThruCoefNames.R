
# This code is surplus
catln('coefFamCv')
coefFamCv=str_replace_all( coefNames,"\\|.*?(:|$)","\\1")
uFamStr=unique(coefFamCv)
# famCoefList$inxCoefs, famCoefList$names, famCoefList$coefFacList
famCoefList=list()
for (thisFam in uFamStr){
    inxTheseFamCoefs=which(coefFamCv==thisFam)
    theseNames=coefNames [inxTheseFamCoefs]
    names(inxTheseFamCoefs)=theseNames
    # I think we can walk throu this and create a data.frame

    famCoefList[[thisFam]]=list(inxCoefs=inxTheseFamCoefs,
                                names=theseNames,
                                coefFacList=str_split(theseNames,':'))
}
catln('famCoefList')
print(str(famCoefList))

