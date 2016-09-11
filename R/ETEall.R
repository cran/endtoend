####################################################################################
#THEORETICAL
####################################################################################

ETE <- function(p1,p2,L,N)
{
  if(p1 >= 1 | p1 <= 0) stop("p1 must be a real number in (0,1)")
  if(p2 >= 1 | p2 <= 0) stop("p2 must be a real number in (0,1)")
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")

  cat("   ","\n")
  cat("###########################################","\n")
  cat("END TO END - THEORETICAL RESULTS","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Data success probability         p1 = ",p1),"\n")
  cat(paste("ACK success probability          p2 = ",p2),"\n")
  cat(paste("Maximum number of transmissions  L  = ",L),"\n")
  cat(paste("Number of Hops                   N  = ",N),"\n")
  cat("   ","\n")

  pp = p1*p2

  if(L==Inf)
  {
    ETData = ((1 - p1^N)/(1-p1)) * ((1)/(p1*p2)^N)
    ETACK  = ((p1^N)*((1 - p2^N)/(1-p2))) * ((1)/(p1*p2)^N)
    ETS    = ETData + ETACK

    # ETDataV = diff(c(0,((1 - p1^seq(1,N))/(1-p1)) * ((1)/(p1*p2)^seq(1,N))))
    # ETACKV  = diff(c(0,((p1^seq(1,N))*((1 - p2^seq(1,N))/(1-p2))) * ((1)/(p1*p2)^seq(1,N))))
    # ETSV    = ETDataT + ETACKT

    REC_ETData = p1*((1 - p1^N)/(1-p1)) * ((1)/(p1*p2)^N)
    REC_ETACK  = p2 * ((p1^N)*((1 - p2^N)/(1-p2))) * ((1)/(p1*p2)^N)
    REC_ETS    = REC_ETData + REC_ETACK

    # REC_ETDataV = diff(c(0,p1*((1 - p1^seq(1,N))/(1-p1)) * ((1)/(p1*p2)^seq(1,N))))
    # REC_ETACKV  = diff(c(0,p2 * ((p1^seq(1,N))*((1 - p2^seq(1,N))/(1-p2))) * ((1)/(p1*p2)^seq(1,N))))
    # REC_ETSV    = REC_ETDataV + REC_ETACKV
  }else
  {
    if((p1*p2)<0.05)
    {
      ETData = L/(1-p1)
      ETACK  = 0
      ETS    = ETData + ETACK

      REC_ETData = (p1*L)/(1-p1)
      REC_ETACK  = 0
      REC_ETS    = (p1*L)/(1-p1)
    }else{
    ETData = ((1 - p1^N)/(1-p1)) * ((1 - (1-(p1*p2)^N)^L)/(p1*p2)^N)
    ETACK  = ((p1^N)*((1 - p2^N)/(1-p2))) * ((1 - (1-(p1*p2)^N)^L)/(p1*p2)^N)
    ETS    = ETData + ETACK

    # ETDataV = diff(c(0,((1 - p1^seq(1,N))/(1-p1)) * ((1 - (1-(p1*p2)^seq(1,N))^L)/(p1*p2)^seq(1,N))))
    # ETACKV  = diff(c(0,((p1^seq(1,N))*((1 - p2^seq(1,N))/(1-p2))) * ((1 - (1-(p1*p2)^seq(1,N))^L)/(p1*p2)^seq(1,N))))
    # ETSV    = ETDataT + ETACKT

    REC_ETData = p1*((1 - p1^N)/(1-p1)) * ((1 - (1-(p1*p2)^N)^L)/(p1*p2)^N)
    REC_ETACK  = p2 * ((p1^N)*((1 - p2^N)/(1-p2))) * ((1 - (1-(p1*p2)^N)^L)/(p1*p2)^N)
    REC_ETS    = REC_ETData + REC_ETACK

    # REC_ETDataV = diff(c(0,p1*((1 - p1^seq(1,N))/(1-p1)) * ((1 - (1-(p1*p2)^seq(1,N))^L)/(p1*p2)^seq(1,N))))
    # REC_ETACKV  = diff(c(0,p2 * ((p1^seq(1,N))*((1 - p2^seq(1,N))/(1-p2))) * ((1 - (1-(p1*p2)^seq(1,N))^L)/(p1*p2)^seq(1,N))))
    # REC_ETSV    = REC_ETDataV + REC_ETACKV
    }
  }

  # Pr1    = 1-(1-p1)^L
  PrS    = 1-(1-p1^N)^L
  # PrSV   = 1-(1-p1^seq(1,N))^L

  res    = round(matrix(data = c(c(PrS),c(ETData),c(ETACK),c(ETS),c(REC_ETData),c(REC_ETACK),c(REC_ETS)),nrow = 7,ncol = 1,byrow = T),4)
  rownames(res)=c("Success Probability",
                  "Expected Data Transmissions",
                  "Expected ACK Transmissions",
                  "Expected Total Transmissions",
                  "Expected Data Receptions",
                  "Expected ACK Receptions","Expected Total Receptions")
  colnames(res)= c("Total")
  return(res)
}

# ETE(p1 = 0.3,p2 = 0.4,L = Inf,N = 3)
# ETE(p1 = 0.3,p2 = 0.4,L = 4,N = 3)
####################################################################################
#MONTE CARLO
####################################################################################

MCETE = function(p1,p2,L,N,M=5000)
{
  if(p1 >= 1 | p1 <= 0) stop("p1 must be a real number in (0,1)")
  if(p2 >= 1 | p2 <= 0) stop("p2 must be a real number in (0,1)")
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")
  if(M%%1!=0 | M<1) stop("N must be a positive integer")

  cat("   ","\n")
  cat("###########################################","\n")
  cat("END TO END - MONTE CARLO SIMULATION RESULTS","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Data success probability         p1 = ",p1),"\n")
  cat(paste("ACK success probability          p2 = ",p2),"\n")
  cat(paste("Maximum number of transmissions  L  = ",L),"\n")
  cat(paste("Number of Hops                   N  = ",N),"\n")
  cat(paste("Monte Carlo Simulations          M  = ",M),"\n")
  cat("   ","\n")

  prog2 = function(p1,p2,L,N)
  {
    pos     = 1
    trans   = 0
    ack     = 0
    maxpos  = 1
    ok      = FALSE
    okACK   = FALSE
    arr     = 0
    abj     = 0
    failUP  = 0
    failACK = 0
    failT   = 0
    chegou  = 0

    #mientras hayan intentos y mientas no haya llegado
    while(failT < L && okACK == FALSE)
    {
      #mientas np se llegue al final y hayan intentos
      back   = FALSE

      while(pos<=N && failT < L)
      {
        maxpos = max(maxpos,pos)
        u1 = runif(1)
        trans = trans +1
        if(u1<p1){arr = arr + 1;pos = pos + 1}else{pos = 1;failUP = failUP + 1;failT = failT + 1}
      }

      if(pos == N+1)
      {
        maxpos = N+1
        posACK = pos
        ok = TRUE
        chegou=1
        while(posACK >1 && back == FALSE)
        {
          u2 = runif(1)
          ack = ack +1
          if(u2<p2){abj = abj + 1;posACK = posACK - 1}else{back = TRUE;pos = 1;failACK = failACK + 1;failT = failT + 1}
          if(posACK == 1){okACK = TRUE}
        }
      }
    }
    return(list(pos=pos,tDATA = trans,tACK = ack,
                trans = trans+ack,RDATA = arr,
                RACK = abj,Rtrans = arr+abj,ok=chegou))
    cat("   ","\n")
  }

  resul = matrix(data = NA,nrow = M,ncol = 8)
  pb <- txtProgressBar(min = 0, max = M, style = 3)

  for(k in 1:M)
  {
    run = prog2(p1,p2,L,N)
    resul[k,1]=run$pos
    resul[k,2]=run$tDATA
    resul[k,3]=run$tACK
    resul[k,4]=run$trans
    resul[k,5]=run$RDATA
    resul[k,6]=run$RACK
    resul[k,7]=run$Rtrans
    resul[k,8]=run$ok
    setTxtProgressBar(pb, k)
  }
  close(pb)

  success = sum(resul[,8])/M
  tarr = mean(resul[,2])
  tabj = mean(resul[,3])
  tt = tarr + tabj
  rarr = mean(resul[,5])
  rabj = mean(resul[,6])
  rt = rarr + rabj

  table = round(rbind(success,tarr,tabj,tt,rarr,rabj,rt),4)
  rownames(table)=c("MC Success Probability","MC Mean Data Transmissions","MC Mean ACK Transmissions","MC Mean Total Transmissions","MC Mean Data Receptions","MC Mean ACK Receptions","MC Mean Total Receptions")
  colnames(table)= c("Total")
  return(table)
  cat("   ","\n")
}

#RUN

#Theoretical
# ETE(p1=0.10,p2=0.10,L=17,N=10)
# #MonteCarlo Simulations
# MCETE(p1=0.10,p2=0.10,L=17,N=10,M=5000)

# 1-(0.15^2)^10
# 17/(1-0.15)

# N=10
# pseq = seq(from = 0.1,to = 0.9,by = 0.1)
# Lseq = c(27,13,8,6,5,4,3,2,2)
#
# SIM = matrix(data = NA,nrow = 9,ncol = 13,byrow = TRUE)
#
# for(i in 1:9)
# {
#   p1 = p2 = pseq[i]
#   L = Lseq[i]
#   SIM[i,] = c(ETE(p1=p1,p2=p2,L=L,N=10),MCETE(p1=p1,p2=p2,L,N=10,M=10^6)[2:7])
# }
#
# SIM[,2:7]-SIM[,8:13]
#
# TABLE = cbind(pseq,Lseq,SIM)
# library(xtable)
# xtable(TABLE)


# N=5
# pseq = seq(from = 0.1,to = 0.9,by = 0.1)
# Lseq = c(27,13,8,6,5,4,3,2,2)
#
# SIM2 = matrix(data = NA,nrow = 9,ncol = 13,byrow = TRUE)
#
# for(i in 1:9)
# {
#   p1 = p2 = pseq[i]
#   L = Lseq[i]
#   SIM2[i,] = c(ETE(p1=p1,p2=p2,L=L,N=5),MCETE(p1=p1,p2=p2,L,N=5,M=10^6)[2:7])
# }
#
# SIM2[,2:7]-SIM2[,8:13]
#
# TABLE2 = cbind(pseq,Lseq,SIM2)
# library(xtable)
# xtable(TABLE2)
#
#
#
#
# N=7
# pseq = seq(from = 0.1,to = 0.9,by = 0.1)
# Lseq = c(27,13,8,6,5,4,3,2,2)
#
# SIM3 = matrix(data = NA,nrow = 9,ncol = 13,byrow = TRUE)
#
# for(i in 1:9)
# {
#   p1 = p2 = pseq[i]
#   L = Lseq[i]
#   SIM3[i,] = c(ETE(p1=p1,p2=p2,L=L,N=7),MCETE(p1=p1,p2=p2,L,N=7,M=10^6)[2:7])
# }
#
# SIM3[,2:7]-SIM3[,8:13]
#
# TABLE3 = cbind(pseq,Lseq,SIM3)
# library(xtable)
# xtable(TABLE3)







# N=7
# pseq = seq(from = 0.1,to = 0.9,by = 0.1)
# Lseq = c(27,13,8,6,5,4,3,2,2)
#
# SIM4 = matrix(data = NA,nrow = 9,ncol = 13,byrow = TRUE)
#
# for(i in 1:9)
# {
#   p1 = p2 = pseq[i]
#   L = Lseq[i]
#   SIM4[i,] = c(ETE(p1=p1,p2=p2,L=L,N=7),MCETE(p1=p1,p2=p2,L,N=7,M=10^6)[2:7])
# }
# SIM4[,1]
#
# ERROR4 = (abs(SIM4[,2:7]-SIM4[,8:13]))
#
#
# boxplot(abs(SIM4[,2:7]-SIM4[,8:13]))
#
#
# TABLE4 = cbind(pseq,Lseq,SIM4)
# library(xtable)
# xtable(TABLE4,digits = 4)
#
