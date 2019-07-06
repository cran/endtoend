globalVariables(c("..count..","..density..","ind", "values",
                  "Success Probability", "Expected Data Transmissions",
                  "Expected ACK Transmissions", "Expected Total Transmissions",
                  "Expected Data Receptions", "Expected ACK Receptions",
                  "Expected Total Receptions"))

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

####################################################################################
#no output
####################################################################################


ETE0 = function (p1, p2, L, N){
  pp = p1 * p2
  if (L == Inf) {
    ETData = ((1 - p1^N)/(1 - p1)) * ((1)/(p1 * p2)^N)
    ETACK = ((p1^N) * ((1 - p2^N)/(1 - p2))) * ((1)/(p1 *
                                                       p2)^N)
    ETS = ETData + ETACK
    REC_ETData = p1 * ((1 - p1^N)/(1 - p1)) * ((1)/(p1 *
                                                      p2)^N)
    REC_ETACK = p2 * ((p1^N) * ((1 - p2^N)/(1 - p2))) * ((1)/(p1 *
                                                                p2)^N)
    REC_ETS = REC_ETData + REC_ETACK
  }else {
    if ((p1 * p2) < 0.05) {
      ETData = L/(1 - p1)
      ETACK = 0
      ETS = ETData + ETACK
      REC_ETData = (p1 * L)/(1 - p1)
      REC_ETACK = 0
      REC_ETS = (p1 * L)/(1 - p1)
    }else {
      ETData = ((1 - p1^N)/(1 - p1)) * ((1 - (1 - (p1 *
                                                     p2)^N)^L)/(p1 * p2)^N)
      ETACK = ((p1^N) * ((1 - p2^N)/(1 - p2))) * ((1 -
                                                     (1 - (p1 * p2)^N)^L)/(p1 * p2)^N)
      ETS = ETData + ETACK
      REC_ETData = p1 * ((1 - p1^N)/(1 - p1)) * ((1 - (1 -
                                                         (p1 * p2)^N)^L)/(p1 * p2)^N)
      REC_ETACK = p2 * ((p1^N) * ((1 - p2^N)/(1 - p2))) *
        ((1 - (1 - (p1 * p2)^N)^L)/(p1 * p2)^N)
      REC_ETS = REC_ETData + REC_ETACK
    }
  }
  PrS = 1 - (1 - p1^N)^L
  res = round(matrix(data = c(c(PrS), c(ETData), c(ETACK),
                              c(ETS), c(REC_ETData), c(REC_ETACK), c(REC_ETS)), nrow = 7,
                     ncol = 1, byrow = T), 4)
  # rownames(res) = c("Success Probability", "Expected Data Transmissions",
  #                   "Expected ACK Transmissions", "Expected Total Transmissions",
  #                   "Expected Data Receptions", "Expected ACK Receptions",
  #                   "Expected Total Receptions")
  # colnames(res) = c("Total")
  return(res)
}


stochastic_ETE = function(dist1,p11,p12,dist2,p21,p22,L,N,M=10^5,printout=TRUE,plotspdf=TRUE){
  if(L != Inf && (L%%1!=0 | L<0)) stop("L must be a positive integer")
  if(N%%1!=0 | N<1) stop("N must be a positive integer")
  if(dist1 == "uniform"){
    p1 = runif(M,p11,p12)
  }else{
    if(dist1 == "beta"){
      p1 = rbeta(M,p11,p12)
    }else{
      stop("p1 distribution must be either 'uniform' or 'beta'.")
    }
  }
  if(dist2 == "uniform"){
    p2 = runif(M,p21,p22)
  }else{
    if(dist2 == "beta"){
      p2 = rbeta(M,p21,p22)
    }else{
      stop("p1 distribution must be either 'uniform' or 'beta'.")
    }
  }

  out = apply(X = cbind(p1,p2),MARGIN = 1, function(x) ETE0(x[1], x[2], L, N))
  outsum = matrix(apply(out,1,mean),7,1)
  rownames(outsum) = c("Success Probability", "Expected Data Transmissions",
                       "Expected ACK Transmissions", "Expected Total Transmissions",
                       "Expected Data Receptions", "Expected ACK Receptions",
                       "Expected Total Receptions")
  colnames(outsum) = c("Total")

  df = data.frame(p1,p2,t(out))
  colnames(df) = c("p1","p2",rownames(outsum))

  #print(summary(df))
  #install.packages("pastecs")
  #library(pastecs)

  stats = as.data.frame(t(stat.desc(df))[,-c(1:3)])

  p1 = ggplot(df,aes(p1)) +
    geom_histogram(aes(y=..density.., fill = "p1"),alpha = 0.4,color="gray40",breaks = seq(min(p1),max(p1),length.out = round(diff(range(p1))/0.0625)+1)) +
    geom_histogram(aes(x = p2, y = ..density..,fill = "p2"),alpha = 0.4,color="gray40",breaks = seq(min(p2),max(p2),length.out = round(diff(range(p2))/0.0625)+1))+
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(fill = "")+
    xlab("probability") +
    ylab("density") +
    ggtitle("Data and ACK Success Probabilities") +
    xlim(0,1) +
    theme_classic()

  p2 = ggplot(df,aes(x=`Success Probability`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 1,
                   breaks = seq(min(df$`Success Probability`),max(df$`Success Probability`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_vline(xintercept = outsum[1],color="gray40",lwd=1.2)  +
    ggtitle(paste("Success Probability (mean = ",round(outsum[1],3),")",sep = "")) +
    ylab("relative frequency")+ theme_classic()


  p3 = ggplot(df,aes(x=`Expected Data Transmissions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 2,
                   breaks = seq(min(df$`Expected Data Transmissions`),max(df$`Expected Data Transmissions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[2],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Data Transmissions (mean = ",round(outsum[2],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p4 = ggplot(df,aes(x=`Expected ACK Transmissions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 3,
                   breaks = seq(min(df$`Expected ACK Transmissions`),max(df$`Expected ACK Transmissions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[3],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected ACK Transmissions (mean = ",round(outsum[3],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p5 = ggplot(df,aes(x=`Expected Total Transmissions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 4,
                   breaks = seq(min(df$`Expected Total Transmissions`),max(df$`Expected Total Transmissions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[4],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Total Transmissions (mean = ",round(outsum[4],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p6 = ggplot(df,aes(x=`Expected Data Receptions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 5,
                   breaks = seq(min(df$`Expected Data Receptions`),max(df$`Expected Data Receptions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[5],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Data Receptions (mean = ",round(outsum[5],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  p7 = ggplot(df,aes(x=`Expected ACK Receptions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 6,
                   breaks = seq(min(df$`Expected ACK Receptions`),max(df$`Expected ACK Receptions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[6],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected ACK Receptions (mean = ",round(outsum[6],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()


  p8 = ggplot(df,aes(x=`Expected Total Receptions`)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)),alpha = 0.2,color="gray40",fill = 7,
                   breaks = seq(min(df$`Expected Total Receptions`),max(df$`Expected Total Receptions`),length.out = 17)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_vline(xintercept = outsum[7],color="gray40",lwd=1.2)  +
    ggtitle(paste("Expected Total Receptions (mean = ",round(outsum[7],3),")",sep = "")) +
    ylab("relative frequency") +
    theme_classic()

  df2 = stack(df[,c(4:9)])

  p9 = ggplot(df2, aes(x = ind, y = values)) +
    geom_boxplot(fill = rev(2:7),alpha = 0.2,color="gray40") +
    coord_flip() +
    #scale_y_continuous(trans='sqrt') +
    xlab("") +
    scale_x_discrete(limits = rev(levels(df2$ind))) +
    theme_classic()

  if(isTRUE(printout)){
    cat(paste("Monte Carlo simulations          M  = ", M), "\n")
    cat(paste("Maximum number of transmissions  L  = ", L), "\n")
    cat(paste("Number of Hops                   N  = ", N), "\n")
    print(stats)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    print(p9)
  }

  if(isTRUE(plotspdf)){
    plotsPath = paste("ETE",format(Sys.time(),"%d%m%y_%H%M%S"),".pdf",sep="")
    pdf(file=plotsPath,width = 8.27,height = 5.83)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    print(p9)
    dev.off()
    print(paste("Plots file ",plotsPath," saved in working directory ",getwd(),".",sep = ""))
  }
  return(list(data=df,stats = stats))
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
