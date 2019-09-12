### >>>  FactorialHMNL_CreateData
#' Title
#'
#'
#' @param R  number of respondants
#' @param ntot number of total observations
#' @param seed rando seed
#'
#' @return
#' @export
#'
# myFactorialHMNL_CreateData=function(R,ntot=100,seed=123){
# Setup nu ----
stop('Depricated')
    R=5



    if (!exists(".Random.seed")) runif(30)
    savedSeed <- .Random.seed
    on.exit({.Random.seed <- savedSeed; cat('Random seed restored.\n')})




    set.seed(seed)

    # C= 6
    #
    sylLabs=c('ba','da','bi','di','bu','du')
    sylDfr= data.frame( syl= factor(sylLabs, levels = sylLabs))
    sylDfr$cn= as.factor(as.character(substr(sylDfr$syl,1,1)))
    sylDfr$vow=as.factor(as.character(substr(sylDfr$syl,2,2)))

    # met x be F1, y be F2.  Q=1 be stressed, Q=2 be unstressed Q=3 be reduced.
    # # for pars = c('One','F1','F2','Q1','Q2','Q3')
    nCn=2
    nVw=3
    nStimCols= 6
    # Standardized F1 and F2
    # BetaMu Arrays
    BetaMu=array(0,c(nCn,nVow,nStimCols),
                 dimnames=list(
                     cn=levels(sylDfr$cn),
                     vow=levels(sylDfr$vow),
                     stimPar=c('One','F1','F2','Q1','Q2','Q3')
                 )
    )

    # Leave  BetaMu[ , , 'One']  at zero.

  BetaMu['b', ,'F1'] = BetaMu['b', ,'F1'] -.1
  BetaMu['d', ,'F1'] = BetaMu['d', ,'F1'] +.1
  BetaMu[, 'i','F1'] = BetaMu[, 'i','F1'] -.5
  BetaMu[, 'a','F1'] = BetaMu[, 'a','F1'] +1
  BetaMu[, 'u','F1'] = BetaMu[, 'u','F1'] -.5

  BetaMu['b', ,'F2'] = BetaMu['b', ,'F2'] -.2
  BetaMu['d', ,'F2'] = BetaMu['d', ,'F2'] +.2
  BetaMu[, 'i','F2'] = BetaMu[, 'i','F2'] +1
  BetaMu[, 'a','F2'] = BetaMu[, 'a','F2'] +0
  BetaMu[, 'u','F2'] = BetaMu[, 'u','F2'] -1


    # the Q values
  BetaMu['b','i' ,'Q1'] =BetaMu['b','i' ,'Q1'] +.05
  BetaMu['d','a' ,'Q2'] = BetaMu['d','a' ,'Q2'] +0
  BetaMu['b', 'u','Q3'] =  BetaMu['b', 'u','Q3'] -.05
     cat('BetaMu filled \n')
    print(BetaMu)
     # R
    # Random effects per subject per parameter -----
    BetaSj=list()

    for (r in 1:R){
        BetaSj[[r]]=BetaMu
        for (dname in dimnames(BetaMu)$stimPar){
             if (dname=='One'){
                 rng=.1
             }else{
                 rng=range(as.vector(BetaMu[,,dname]))
             }
             BetaSj[[r]][,,dname]=BetaMu[,,dname] + .1*rng*rnorm(nCn*nVw)
        }

    }
    BetaSj
    BetaSjCoefMat= t(sapply(BetaSj,as.vector))

    }
    #
    # dfrWide=data.frame(syl=sample(c('ba','da','bi','di','bu','du'),
#     #                               ntot,replace=TRUE),x=1:ntot,y=sample(1:ntot))
#     #dbg set x and y to make it easier to see contrasts
#     #  Added to test factorial stimuli the factored stimulus variable m has no effect on responses.
#     dfrWide=data.frame(syl=sample(c('ba','da','bi','di','bu','du'),ntot,replace=TRUE),
#                        x=rnorm(ntot), y=rnorm(ntot),
#                        m =as.factor(sample(1:3,ntot,replace=TRUE))
#     )
#     #  build a chosen as clone of syl to use as a choice variable
#     dfrWide$chosen=dfrWide$syl
#     dfrWide$cn=as.factor(as.character(substr(dfrWide$syl,1,1)))
#     dfrWide$vow=as.factor(as.character(substr(dfrWide$syl,2,2)))
#     # Build in some correlation with a model...
#     # Very sloppy but should be good enough for test
#     if (add.random){
#         dfrWide$id=sample(N,ntot,replace=TRUE)
#         # Note integer sample must have integer N , not (1:N)
#         rand_bx=sapply(N,rnorm)
#         rand_iy=sapply(N,rnorm)
#     }
#
#     inx_b= dfrWide$cn=='b'
#     inx_i=dfrWide$vow=='i'
#     # The random perturbations in any
#     ubx={ if(add.random) rand_bx[dfrWide$id[inx_b]] else 0}
#     uiy={ if(add.random) rand_iy[dfrWide$id[inx_i]] else 0}
#
#     dfrWide$x[inx_b]=dfrWide$x[inx_b]-(2+ ubx )
#     dfrWide$y[inx_i]=dfrWide$y[inx_i]+(2+  uiy )
#
#     #   print(dfrWide)
#     # Restore seed to a random value (based onclock)
#     # Example from ?set.seed
#     # Should be taken care of by on exit
#     #   rm(.Random.seed)
#     #   runif(1)
#     #  .Random.seed[1:6]
# #     return(dfrWide)
# # }
