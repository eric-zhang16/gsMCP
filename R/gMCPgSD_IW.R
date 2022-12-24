#' Performs a graph based multiple test procedure for a given graph and
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
#' @return a dataframe with column 'R' indicating whether each endpoint can be rejected, in the same order as input vector h
#' @import gMCP
#' @import gsDesign
#' @export
#'
#' @examples
#'
gMCPgSD_MS <- function(g,w.start,t,h,p,alpha=0.025,timing,sfpar,debug){
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
        if(ss==1){
          alpha.tmp[j] <- spend.tmp[j]
        } else {
          tmp.time <- timing[,j]
          tmp.id <- which(!is.na(tmp.time))
          if(sum(tmp.id==t.tmp)==1){
            t.tmp2 <- which(tmp.id==t.tmp)
            tmp.time <- tmp.time[!is.na(tmp.time)]
            ct <- gsDesign(k = length(tmp.time), test.type = 1, alpha = spend.tmp[j], sfu = sfHSD,sfupar=sfpar[j],timing=tmp.time)$upper$bound[t.tmp2]
            alpha.tmp[j] <- pnorm(ct, mean = 0, sd = 1, lower.tail = F, log.p = FALSE)
          } else {
            alpha.tmp[j]<-0
          }

        }

      }
    }
    if(debug==1)print(paste(c('Bonferroni threholds=',alpha.tmp),collapse=" "))
    ## Step 2:Check if p-value <= nominal significance ##
    id <- which( p.tmp <= alpha.tmp) # Hypothesis that are rejected #
    jd <- c(jd,id)
    if(length(jd)==0){
      id.c <- h    # id.c for un-rejected #
    } else {
      id.c <- h[!h%in%jd]
    }
    if(length(id)>0){
      if(length(id.c)>0){
        if(debug==1)print(paste(c(id,' are rejected by Bonferroni,',id.c,' are remaining, update graph and continue Bonferroni test'),collapse=" "))
        l.jd <- length(jd)
        l.jd2 <- length(h)-l.jd
        w.tmp <- w.all[rowSums(as.matrix(w.all[,jd])==0)==l.jd &rowSums(as.matrix(w.all[,1:length(h)][,-jd]==1))==l.jd2 ,(length(h)+1):(2*length(h))]
      } else {
        if(debug==1)print(paste(c(id,' are rejected by Bonferroni,no h is remaining, stop and all h are rejected'),collapse=" "))
        cont <- FALSE
      }

    } else if(length(id.c)>=2){ # if  >=2 un-rejected  #
      # Update w to extract weights under id.c #
      if(debug==1)print(paste(c('No one rejected by Bonferroni, remaining hypothesis=',length(id.c),',check unrejected using Simes'),collapse=" "))
      l.c <- length(id.c)
      l.c2 <- length(h)-l.c
      w.c <- w.all[rowSums(as.matrix(w.all[,1:length(h)][,-id.c]==0))==l.c2 ,]
      if(debug==1)print(w.c)
      simes.rej <- 0
      for(ii in id.c){ # Check each hypothesis (ii) in id.c #
        w.ii <- w.c[as.matrix(w.c[,1:length(h)][,ii]==1) ,] # For each hypothesis (ii) in id.c, extract weights for combo involves that hypothesis #
        rej.ii <- rep(0,dim(w.ii)[1])
        cnt.jj <-0
        for(jj in 1:dim(w.ii)[1]){ # For each combo (jj) involves ii, check if any hypothesis (jjj) can be rejected #
          cnt.jj <- cnt.jj + 1
          combo.j <- w.ii[,1:length(h)][jj,]
          w.j <- w.ii[,(length(h)+1):(2*length(h))][jj,]
          id.jj <- which(combo.j==1) # hopthesis in combo jj, id no. refer to h #
          w.jj <- w.j[id.jj]
          for(id.jjj in id.jj){ # To see if any jjj in combo jj can be rejected #
            ## id.jjj: id no. refer to h
            id2.jjj <- which(p.tmp[id.jj] <= p.tmp[id.jjj]) # Find hpothesis with p-value <= p.jjj #
            alpha.spend.jjj <- sum(w.jj[id2.jjj])*alpha
            if(alpha.spend.jjj>0){
              if(ss==1){
                alpha.jjj <- alpha.spend.jjj
              } else {
                tmp.time <- timing[,id.jjj]
                tmp.id <- which(!is.na(tmp.time))
                if(sum(tmp.id==t.tmp)==1){
                  t.tmp2 <- which(tmp.id==t.tmp)
                  tmp.time <- tmp.time[!is.na(tmp.time)]
                  ct <-gsDesign(k = length(tmp.time), test.type = 1, alpha = alpha.spend.jjj, sfu = sfHSD,sfupar=sfpar[id.jjj],timing=tmp.time)$upper$bound[t.tmp2]
                  alpha.jjj <- pnorm(ct, mean = 0, sd = 1, lower.tail = F, log.p = FALSE)
                } else {
                  alpha.jjj<-0
                }

              }

              if(p.tmp[id.jjj]<=alpha.jjj){
                rej.ii[cnt.jj] <- 1
                if(debug==1)print(paste(c('For',ii,':',id.jjj,'in',id.jj,'rejected by',alpha.jjj,'with weights',w.jj[id2.jjj]),collapse=" "))
                break
              } else {
                if(debug==1)print(paste(c('For',ii,':',id.jjj,'in',id.jj,'NOT rejected by',alpha.jjj,'with weights',w.jj[id2.jjj]),collapse=" "))
              }
            }

          } # end of combo jj #
          if(rej.ii[cnt.jj]==0)break

        } # end of all combo involves ii
        if(sum(rej.ii)==dim(w.ii)[1]){ # ii is rejected #
          jd <- c(ii,jd)
          simes.rej <- 1
          if(debug==1)print(paste(c(ii,'is rejected by Simes'),collapse=" "))
          break
        }

      } # end of check all ii #

      if(simes.rej==0){
        if(debug==1)print(paste('No one rejected by Simes'))
        t.tmp <- t.tmp + 1
        gonext=1
        if(t.tmp>length(t)){
          cont <- FALSE
        }
      } else { # Someone rejected by Simes #
        if(length(jd)==length(h)){
          cont <- FALSE
        } else {
          l.jd <- length(jd)
          l.jd2 <- length(h)-l.jd
          w.tmp <- w.all[rowSums(as.matrix(w.all[,jd])==0)==l.jd & rowSums(as.matrix(w.all[,1:length(h)][,-jd]==1))==l.jd2 ,(length(h)+1):(2*length(h))]
        }
      }
    } else {  # if unrejected < 2, go to next IA #
      if(debug==1)print(paste('No one rejected by Bonferroni and unrejected from Bonferroni test < 2, go to next step'))
      t.tmp <- t.tmp + 1
      gonext=1
      if(t.tmp>length(t)){
        cont <- FALSE
      }

    }

  }
  
  res <- data.frame(H=h,R=jd)
  return(res)
}
