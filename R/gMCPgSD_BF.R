#' Performs a graph based multiple test procedure using Bonferroni test for a given graph and
#' unadjusted p-values
#'
#' @param g Matrix of local alpha to be passed. Each element defined how much of the local alpha reserved for the hypothesis corresponding to its row index is passed on to the hypothesis corresponding to its column index.
#' @param w.start Vector of initial weights
#' @param t Stage index
#' @param h Enpoint index
#' @param p Vector of unadjusted p-values

#' @param alpha A numeric specifying the maximal allowed type one error rate
#' @param timing Interim analysis timing
#' @param sfpar Spending function parameter
#' @param debug Debug indicator
#'
#' @return a dataframe  to indicate which hypotheses are rejected. 1=rejected, 0=not rejected
#' @import gMCP
#' @import gsDesign
#' @export
#'
#' @examples
#'
gMCPgSD_BF <- function(g,w.start,t,h,p,alpha=0.025,timing,sfpar,debug){
  w.all<-gMCP::generateWeights(g,w.start)
  w.tmp<-w.start
  jd<-NULL
  cont <- TRUE
  ss<-0
  if(length(t)==1) ss<-1
  ## Step 0: initialize IA ##
  t.tmp<-1
  gonext=0
  if(debug==1)print(paste('t=',t.tmp,sep=''))
  while(cont){
    if(gonext==1 & debug==1) print(paste('t=',t.tmp,sep=''))
    gonext=0
    p.tmp <- p[t.tmp,]

    spend.tmp <- alpha*w.tmp
    alpha.tmp <- rep(0,length(h))
    ## Step 1:Compute nominal significance level for each hypothesis test ##
    for(j in 1:length(h)){
      if(spend.tmp[j]>0){
        if(ss==0){
          tmp.time <- timing[,j]
          tmp.id <- which(!is.na(tmp.time))
          if(sum(tmp.id==t.tmp)==1){
            t.tmp2 <- which(tmp.id==t.tmp)
            tmp.time <- tmp.time[!is.na(tmp.time)]
            ct<-gsDesign(k = length(tmp.time), test.type = 1, alpha = spend.tmp[j], sfu = sfHSD,sfupar=sfpar[j],timing=tmp.time)$upper$bound[t.tmp2]
            alpha.tmp[j]<- pnorm(ct, mean = 0, sd = 1, lower.tail = F, log.p = FALSE)
          } else {
            alpha.tmp[j]<-0
          }

          #alpha.tmp[j]<- sfLDOF(spend.tmp[j],timing,sf)$spend[t.tmp]
        } else {
          alpha.tmp[j] <- spend.tmp[j]
        }

      }
    }
    if(debug==1)print(alpha.tmp)
    ## Step 2:Check if p-value <= nominal significance ##
    id <- which( p.tmp <= alpha.tmp)
    jd <- c(jd,id)
    # If there is p-value passed,update w.tmp #
    if(length(id)>0){
      if(length(jd)<length(h)){
        l.jd <- length(jd)
        l.jd2 <- length(h)-l.jd
        w.tmp <- w.all[rowSums(as.matrix(w.all[,jd])==0)==l.jd &rowSums(as.matrix(w.all[,1:length(h)][,-jd]==1))==l.jd2 ,(length(h)+1):(2*length(h))]
      } else {
        if(debug==1)print('All h are rejected, stop')
        cont <- FALSE
      }

    } else { # If there is no p-value passed,go to next step #
      t.tmp <- t.tmp + 1
      gonext=1
      if(t.tmp>length(t)){
        cont <- FALSE
      }
    }

  }

  res <- matrix(0, nr=1, nc = length(h))
  res.name <- rep(NA,length(h))
  for(jj in 1:length(h)){
    res[1,jj] <- (jj %in% jd)*1
    res.name[jj] <- paste0('H',jj,sep='')
  }

  res <- data.frame(res)
  colnames(res) <- res.name

  if(length(jd)==0) print("No hypothesis is rejected")
  return(res)

}
