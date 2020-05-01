# 
# if (which(rownames(installed.packages())=="vegan")==0) install.packages("vegan")
# library(vegan)
# 
# ### Example 1 with two communities
# dataset=read.table("ExampleCom.txt", h=T, row.names=1)  ### example with two communities
# dataset=as.data.frame(t(dataset))
# dim(dataset)
# 
# # Example 2 with three communities
# dataset=read.table("ExampleCom3.txt", h=T, row.names=1)
# dataset=as.data.frame(t(dataset))
# dim(dataset)

###############################################################################################
###############################################################################################
############# Function for computing Nielsen estimator of Hurlbert ENS (k=2) See Dauby & Hardy    ##############
###############################################################################################
###############################################################################################


#' Nielsen index
#'
#' Compute Nielsen index
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param x vector of taxa abundances
#' @importFrom vegan diversity
#' 
#' 
#' @example 
#' \dontrun{
#' test <- sample(seq(1, 1000, 1), 30)
#' 
#' nielsen2(test)
#' }
#' 
#' @export
nielsen2 <- function(x) {
  Q <- vegan::diversity(x, index = "simpson")
  
  niel <- ((sum(x)-1)^2)/((1-Q)*(sum(x)+1)*(sum(x)-2)+3-sum(x))
  
  return(niel)
}

#' Nielsen indices for community matrix
#'
#' Compute Nielsen index for a matrix taxa x sample units
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param data matrix or dataframe with species as columns and sample units (e.g. plots) as rows
#' 
#' @importFrom vegan diversity
#' 
#' @export
nielsen <- function(data) {   #### the 'dataset' must be arranged as in the examples.
  
  niel_res <- apply(data, MARGIN = 1, function(x) nielsen2(x[x>0]))
  
  return(niel_res)
}

###############################################################################################
###############################################################################################
############# Function for computing Hurlbert ENS (k)   ##############
############### Need 'vegan' package ########################################
###############################################################################################

#' ENS(k) index
#'
#' Compute Equivalent Number of Species for a given parameter k
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param x vector of taxa abundances
#' @param k_param integer higher than 2 and lower or equal than the total number of individual in the plot
#' 
#' @importFrom vegan diversity
#' 
#' 
#' @example 
#' \dontrun{
#' test <- sample(seq(1, 100, 1), 30)
#' 
#' ensk2(test, k_param = 10)
#' 
#' ensk2(test, k_param = 100000)
#' }
#' 
#' @export
ensk2 <-
  function(x, k_param, verbose = TRUE) {
    
    if (k_param <= sum(x)) {
      Sk <- vegan::rarefy(x, k_param)
      Nt <- as.numeric(Sk)
      for (j in 1:10000) {
        Nt1 <- Sk / (1 - (1 - 1 / Nt) ^ k_param)
        if (Nt1 - Nt < 0.0001) {
          break
        }
        else{
          Nt <- Nt1
        }
      }
      return(as.numeric(Nt))
    }else{
      
      if(verbose) message("k parameter cannot be higher than the total number of individual in the sample")
      return(NA)
    }
  }

#' ENS(k) indices for community matrix
#'
#' Compute ENS(k) indices for a matrix taxa x sample units and a vector of ka values
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param data matrix or dataframe with species as columns and sample units (e.g. plots) as rows
#' @param k vector of integer
#' 
#' @importFrom vegan diversity
#' 
#' @export
ensk <-  function(data, k, verbose = TRUE) {      #### fontion pour calculer le ENS de E(k)
  
  ens_list <- list()
  for (k_param in k) {
    cat(" ", k_param)
    
    ens_res <-
      apply(data, MARGIN = 1, function(x)
        ensk2(
          x = x[x > 0],
          k_param = k_param,
          verbose = verbose
        ))
    
    ens_list[[length(ens_list) + 1]] <- ens_res
    
  }
  
  ENS <- do.call("cbind", ens_list)
    
  colnames(ENS) <- paste0("ens_", k)
  
  if (!is.null(rownames(data)))
    rownames(ENS) <- rownames(data)
  return(ENS)
  
  }




###############################################################################################
####################  Function for estimating the ratio of diversity sensu stricto between two communities  using Hurlbert ENS for a given value of parameter K 
####################  See Dauby & Hardy 2012, Figure 4.
####################  Use the function ENSK (see above) 
###########################################################################
###############################################################################################

ensk_ratio <- function(data, k) {
  TableEk <-
    matrix(NA, length(data[, 1]), max(rowSums(data)) - 1)
  for (i in 1:(max(rowSums(data)) - 1))
    TableEk[, i] <- ensk(data, k = i + 1, verbose = FALSE)
  TrueRatio <- matrix(NA, nrow = nrow(data) * (nrow(data) - 1) / 2, ncol =
                        1)
  Names = c()
  w = 1
  for (i in 1:(nrow(data) - 1)) {
    for (j in (i + 1):nrow(data)) {
      A <- i
      B <- j
      test <- c()
      if (k >= min(rowSums(data)))
        stop("k must be at least equal to the minimum sample size")
      
      ratio1 <- TableEk[A, k] / TableEk[B, k]
      
      if (is.na(ratio1) != TRUE) {
        if (ratio1 < 1) {
          A <- j
          B <- i
          ratio1 <- TableEk[A, k] / TableEk[B, k]
        }
        if (round(ratio1 * (k + 1), 0) <= max(rowSums(data))) {
          if (is.na(TableEk[A, ratio1 * (k + 1)]) != TRUE)  {
            for (z in 1:100000) {
              if (round(ratio1 * (k + 1), 0) <= max(rowSums(data)))
              {
                ratio2 = TableEk[A, ratio1 * (k + 1)] / TableEk[B, k + 1]
              } else
              {
                ratio2 <- 0
                break
              }
              test <- cbind(test, ratio2)
              if (ratio2 - ratio1 < 0.0001) {
                break
              } else
              {
                ratio1 <- ratio2
              }
            }
            
            TrueRatio[w, 1] <-
              ratio2
          }
        }
        
        if (!is.null(rownames(data))) {
          com1 <- rownames(data)[i]
          com2 <- rownames(data)[j]
        } else{
          com1 <- "com1"
          com2 <- "com2"
        }
        
        Names <-
          c(Names, paste(com1, ifelse(A > B, A, B), "/", com2, ifelse(A < B, A, B)))
        w <- w + 1
      }
    }
  }
  
  rownames(TrueRatio) <- Names
  
  return(TrueRatio)
}   ### end function


# 
# ENSkRatio(dataset, 40) #### Gives the ratio for each comparison of two communities. Ratio is always higher or equal to 1
# ENSkRatio(dataset, 50) #### If a value of 0 or NA is given, it means the ratio cannot be estimated because the community with the highest diversity has not been sufficiently sampled (see Dauby & Hardy 2012), page 669 
# 

###############################################################################################
###############################################################################################
############# Function for computing Gamma diversity for TWO communities'samples for a serie of k value - NECESSARY for estimating PAIRWISE beta diversity (see below)   ##############
############# When samples vary in size (see page 9-10 of Dauby & Hardy 2011):
############# Rarefaction procedure is used to weight samples equally when computing the gamma component. Here, pooled samples are obtained after rarefying the community sample to the lowest sample size in the dataset.
############## This function is rather slow if samples vary in size and if Kmax is high #################################################################################
###############################################################################################


GammaPairwise <-  function(data=dataset, Kmax=kmax){
       if (dim(data)[1]!=2)
        stop("data must comprise two communities'samples")
       if (Kmax>min(rowSums(data)))
        stop("Kmax must be lower of equal to the lowest sample size observed in all communities")
       library(vegan)
        tabK=c()
       if (rowSums(data[1,])==rowSums(data[2,])) {
          TabSk=c()
                  comGamma=table(c(rep(names(data[2,]),data[2,]),rep(names(data[2,]),data[2,])))
                  for (i in 2:Kmax) {TabSk=c(TabSk,rarefy(comGamma,i))}
                 for (k in 2:Kmax) {
                     Sk=TabSk[k-1]
                     Nt=TabSk[k-1]
                             for (i in 1:10000) {
                                 Nt1=Sk/(1-(1-1/Nt)^(k))
                                 if (Nt1-Nt<0.0001) {break}
                                 else{Nt=Nt1}
                                                       }# end loop i
                            tabK=c(tabK,Nt)
                                  } # end loop k
                                  }
        if (length(tabK)==0) {
          com1 = data[1,]
          subcom1 = com1[,which((com1)>0)]
          com2 = data[2,]
          subcom2= com2[,which((com2)>0)]
          min.size=min(sum(com1),sum(com2))
          com3=t(rbind(com1,com2))
          comGamma=table(c(rep(names(subcom1),subcom1),rep(names(subcom2),subcom2)))


       tabK<-c()
       for (k in 2:Kmax) {
         probaNoAll=c()
             for (f in 1:length(comGamma))          {
                 sp.nb1 <- subset(com3[,1], names(com3[,1])==names(comGamma[f]))
                 sp.nb2 <- subset(com3[,2], names(com3[,2])==names(comGamma[f]))

                 tab2=c()
                 for (l in 0:k) {
                 probaNo=1

                                            prod1 <- dbinom(l,k,1/2)
                                             prod2 <- choose(sum(com1)-sp.nb1,l)/choose(sum(com1),l)
                                             prod3 <- choose(sum(com2)-sp.nb2,k-l)/choose(sum(com2),k-l)

                                   probaNo=probaNo*prod1*prod2*prod3
                                   tab2=c(tab2,probaNo)

                                  }        ## end loop k
                             probaNoAll=c(probaNoAll,sum(tab2))
                                               } ## end f loop - species

               Sk=sum(1-probaNoAll)
                                               ### number of expected species among k individuals
                    Nt=Sk
                            for (i in 1:10000) {
                                Nt1=Sk/(1-(1-1/Nt)^(k))
                                if (Nt1-Nt<0.0001) {break}
                                else{Nt=Nt1}
                                               }
                    tabK=c(tabK,Nt)
                        } #end K loop
                        } #end 'else' loop
                return(tabK)

                                              }  ### end function
# 
# GammaPairwise(dataset, 2)     #### Example for computing Gamma diversity for k comprised between 2 and 10

###############################################################################################
###############################################################################################
############# Function for computing Beta Diversity for a serie of K value   ##############
#############
############# A limitation of Hurlbert's ENS is that k cannot be higher than the total size of the sample constraining the possible relative weight that can be attributed to rare species
###############################################################################################
###############################################################################################


BetaPairwise <-  function(data=dataset, Gamma=tabK, Kmax=kmax){
       if (dim(data)[1]!=2)
        stop("data must comprise two communities'samples")
       if (Kmax>min(rowSums(data)))
        stop("Kmax must be lower of equal to the lowest sample size observed in all communities")
       if (Kmax>(length(Gamma)+1))
        stop("Kmax must be lower of equal to length of Gamma vector")
       library(vegan)
        com1 = data[1,]
        subcom1 = com1[,which((com1)>0)]
        com2 = data[2,]
        subcom2= com2[,which((com2)>0)]
        tab1=c()
         for (k in 2:Kmax)     {
            rar21=rarefy(t(as.matrix(subcom1)), k)
            rar22=rarefy(t(as.matrix(subcom2)), k)
            meanal=(rar21+rar22)/2
            Sk=meanal
            Nt=Sk
            for (j in 1:10000) {
                Nt1=Sk/(1-(1-1/Nt)^k)
                if (Nt1-Nt<0.0001) {break}
                else{Nt=Nt1}
            }
            meanal=Nt
            tab1=cbind(tab1,meanal)
          }

         Tabbeta=c()
          for (g in 1:(Kmax-1)) {
              ratio1=Gamma[g]/tab1[,g]
              if (ratio1<1) {ratio1=1}
              else{ratio1=ratio1}
              for (k in 1:100000){
                      if  (ratio1*g>length(Gamma)) {ratio2=NA}
                          else{
                  ratio2=Gamma[g*ratio1]/tab1[,g]
                  if (ratio2-ratio1<0.0001) {break}
                  else{ratio1=ratio2}
                               }
                                }
              Tabbeta=cbind(Tabbeta, ratio2)
              Tabbeta=as.data.frame(Tabbeta)
                           }
              for (i in 1:length(Tabbeta)) {names(Tabbeta)[i]=paste("K",i+1)}
              return(Tabbeta)
                                                } ### end function

                                
#### Examples
# 
# tabK <- GammaPairwise(dataset,20)   #### Ideally the Kmax value for computing Gamma should the highest (the minimum sample size)
# BetaPairwise(dataset,tabK,10) #### Example for computing ENS beta diversity for k value comprised between 2 and 10
# 
# ### Example for computing pairwise beta diversity
# w=1
# Kmax=20
# Betatable=matrix(data=NA,nrow=nrow(dataset)*(nrow(dataset)-1)/2,ncol=(Kmax-1))
# for (p in 1:(nrow(dataset)-1)) {
#       for (z in (p+1):nrow(dataset)){
#       comFull<-rbind(dataset[p,],dataset[z,])
#       tabK <- GammaPairwise(comFull,20)
#       BETA <- BetaPairwise(comFull,tabK,20)
#       Betatable[w,] <- as.matrix(BETA)
#       w=w+1
#       colnames(Betatable)=colnames(BETA)
#                                 }
#                                 }
# Betatable
# 


      