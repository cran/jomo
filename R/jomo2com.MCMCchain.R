jomo2com.MCMCchain <-
  function(Y.con=NULL, Y.cat=NULL, Y.numcat=NULL, Y2.con=NULL, Y2.cat=NULL, Y2.numcat=NULL,  X=NULL, X2=NULL, Z=NULL, clus, beta.start=NULL, l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, start.imp=NULL, l2.start.imp=NULL, nburn=1000, output=1, out.iter=10) {
    if (is.null(X)) X=matrix(1,max(nrow(Y.cat),nrow(Y.con)),1)
    if (is.null(X2)) X2=matrix(1,max(nrow(Y2.cat),nrow(Y2.con)),1)
    if (is.null(Z)) Z=matrix(1,nrow(X),1)
    if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(max(0,ncol(Y.con))+max(0,(sum(Y.numcat)-length(Y.numcat)))))
    if (is.null(l2.beta.start)) l2.beta.start=matrix(0,ncol(X2),(max(0,ncol(Y2.con))+max(0,(sum(Y2.numcat)-length(Y2.numcat)))))
    if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
    if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
    if (is_tibble(Y.con)) {
      Y.con<-data.frame(Y.con)
      warning("tibbles not supported. Y.con converted to standard data.frame. ")
    }
    if (is_tibble(Y.cat)) {
      Y.cat<-data.frame(Y.cat)
      warning("tibbles not supported. Y.cat converted to standard data.frame. ")
    }
    if (is_tibble(Y2.con)) {
      Y2.con<-data.frame(Y2.con)
      warning("tibbles not supported. Y2.con converted to standard data.frame. ")
    }
    if (is_tibble(Y2.cat)) {
      Y2.cat<-data.frame(Y2.cat)
      warning("tibbles not supported. Y2.cat converted to standard data.frame. ")
    }
    if (is_tibble(X)) {
      X<-data.frame(X)
      warning("tibbles not supported. X converted to standard data.frame. ")
    }
    if (is_tibble(Z)) {
      Z<-data.frame(Z)
      warning("tibbles not supported. Z converted to standard data.frame. ")
    }
    if (is_tibble(X2)) {
      X2<-data.frame(X2)
      warning("tibbles not supported. X2 converted to standard data.frame. ")
    }
    clus<-factor(unlist(clus))
    previous_levels_clus<-levels(clus)
    levels(clus)<-0:(nlevels(clus)-1)
    ncolYcon=max(0,ncol(Y.con))
    ncolY2con=max(0,ncol(Y2.con))
    stopifnot(((!is.null(Y.con))||(!is.null(Y.cat)&!is.null(Y.numcat))),((!is.null(Y2.con))||(!is.null(Y2.cat)&!is.null(Y2.numcat))))
    if (is.null(u.start)) u.start = matrix(0, nlevels(clus), ncol(Z)*(ncolYcon+max(0,(sum(Y.numcat)-length(Y.numcat))))+(ncolY2con+max(0,(sum(Y2.numcat)-length(Y2.numcat)))))
    if (is.null(l2cov.start)) l2cov.start = diag(1, ncol(u.start))
    if (is.null(l2cov.prior)) l2cov.prior = diag(1, ncol(l2cov.start))
    if (!is.null(Y.cat)) {
      isnullcat=0
      previous_levels<-list()
      Y.cat<-data.frame(Y.cat)
      for (i in 1:ncol(Y.cat)) {
        Y.cat[,i]<-factor(Y.cat[,i])
        previous_levels[[i]]<-levels(Y.cat[,i])
        levels(Y.cat[,i])<-1:nlevels(Y.cat[,i])
      }
    } else {
      isnullcat=1
      Y.cat=-999
      Y.numcat=-999
    } 
    if (!is.null(Y2.cat)) {
      isnullcat2=0
      previous_levels2<-list()
      Y2.cat<-data.frame(Y2.cat)
      for (i in 1:ncol(Y2.cat)) {
        Y2.cat[,i]<-factor(Y2.cat[,i])
        previous_levels2[[i]]<-levels(Y2.cat[,i])
        levels(Y2.cat[,i])<-1:nlevels(Y2.cat[,i])
      }
    } else {
      isnullcat2=1
      Y2.cat=-999
      Y2.numcat=-999
    }
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    for (i in 1:ncol(X2)) {
      if (is.factor(X2[,i])) X2[,i]<-as.numeric(X2[,i])
    }
    for (i in 1:ncol(Z)) {
      if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
    }
    if (!is.null(Y.con)) {
      stopifnot(nrow(Y.con)==nrow(clus),nrow(Y.con)==nrow(X), nrow(Z)==nrow(Y.con))
    }
    if (isnullcat==0) {
      stopifnot(nrow(Y.cat)==nrow(clus),nrow(Y.cat)==nrow(X), nrow(Z)==nrow(Y.cat))
    }
    if (!is.null(Y2.con)) {
      stopifnot(nrow(Y2.con)==nrow(clus),nrow(Y2.con)==nrow(X), nrow(Z)==nrow(Y2.con))
    }
    if (isnullcat2==0) {
      stopifnot(nrow(Y2.cat)==nrow(clus),nrow(Y2.cat)==nrow(X), nrow(Z)==nrow(Y2.cat))
    }
    stopifnot(nrow(beta.start)==ncol(X), ncol(beta.start)==(ncolYcon+max(0,(sum(Y.numcat)-length(Y.numcat)))))
    stopifnot(nrow(l2.beta.start)==ncol(X2), ncol(l2.beta.start)==(ncolY2con+max(0,(sum(Y2.numcat)-length(Y2.numcat)))))
    stopifnot(ncol(u.start)==ncol(Z)*(ncolYcon+max(0,(sum(Y.numcat)-length(Y.numcat))))+(ncolY2con+max(0,(sum(Y2.numcat)-length(Y2.numcat)))))
    stopifnot(nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(beta.start))
    stopifnot(nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start))
    stopifnot(ncol(l2cov.start)==ncol(u.start), ncol(l2cov.start)==ncol(l2cov.prior), ncol(l2cov.prior)==nrow(l2cov.prior))
    betait=matrix(0,nrow(beta.start),ncol(beta.start))
    for (i in 1:nrow(beta.start)) {
      for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
    }
    beta2it=matrix(0,nrow(l2.beta.start),ncol(l2.beta.start))
    for (i in 1:nrow(l2.beta.start)) {
      for (j in 1:ncol(l2.beta.start)) beta2it[i,j]=l2.beta.start[i,j]
    }
    covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    for (i in 1:nrow(l1cov.start)) {
      for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
    }   
    uit=matrix(0,nrow(u.start),ncol(u.start))
    for (i in 1:nrow(u.start)) {
      for (j in 1:ncol(u.start)) uit[i,j]=u.start[i,j]
    }
    covuit=matrix(0,nrow(l2cov.start),ncol(l2cov.start))
    for (i in 1:nrow(l2cov.start)) {
      for (j in 1:ncol(l2cov.start)) covuit[i,j]=l2cov.start[i,j]
    }   
    if (!is.null(Y.con)) {
      colnamycon<-colnames(Y.con)
      Y.con<-data.matrix(Y.con)
      storage.mode(Y.con) <- "numeric"  
    }
    if (isnullcat==0) {
      colnamycat<-colnames(Y.cat)
      Y.cat<-data.matrix(Y.cat)
      storage.mode(Y.cat) <- "numeric"  
    }
    if (!is.null(Y2.con)) {
      colnamy2con<-colnames(Y2.con)
      Y2.con<-data.matrix(Y2.con)
      storage.mode(Y2.con) <- "numeric"  
    }
    if (isnullcat2==0) {
      colnamy2cat<-colnames(Y2.cat)
      Y2.cat<-data.matrix(Y2.cat)
      storage.mode(Y2.cat) <- "numeric"  
    }
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    colnamx2<-colnames(X2)
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"  
    stopifnot(!any(is.na(X)))
    Z<-data.matrix(Z)
    storage.mode(Z) <- "numeric"
    stopifnot(!any(is.na(Z)))
    X2<-data.matrix(X2)
    storage.mode(X2) <- "numeric"
    stopifnot(!any(is.na(X2)))
    clus <- matrix(as.integer(levels(clus))[clus], ncol=1)
    if (!is.null(Y.con)&isnullcat==0) {
      Y=cbind(Y.con,Y.cat)
      Yi=cbind(Y.con, matrix(0,nrow(Y.con),(sum(Y.numcat)-length(Y.numcat))))
    } else if (!is.null(Y.con)) {
      Y=Y.con
      Yi=Y.con
    } else {
      Y=Y.cat
      Yi=matrix(0,nrow(Y.cat),(sum(Y.numcat)-length(Y.numcat)))
    }
    n.patterns<-c(0,0)
    if (any(is.na(Y))) {
      if (ncol(Y)==1) {
        miss.pat<-matrix(c(0,1),2,1)
        n.patterns[1]<-2
      } else  {
        miss.pat<-md.pattern.mice(Y, plot=F)
        miss.pat<-miss.pat[,colnames(Y)]
        n.patterns[1]<-nrow(miss.pat)-1
      }
    } else {
      miss.pat<-matrix(0,2,ncol(Y))
      n.patterns[1]<-nrow(miss.pat)-1
    }
    
    miss.pat.id<-rep(0,nrow(Y))
    for (i in 1:nrow(Y)) {
      k <- 1
      flag <- 0
      while ((k <= n.patterns[1]) & (flag == 0)) {
        if (all(!is.na(Y[i,])==miss.pat[k,1:(ncol(miss.pat))])) {
          miss.pat.id[i] <- k
          flag <- 1
        } else {
          k <- k + 1
        }
      }
    }
    
    if (!is.null(Y2.con)&isnullcat2==0) {
      Y2=cbind(Y2.con,Y2.cat)
      Y2i=cbind(Y2.con, matrix(0,nrow(Y2.con),(sum(Y2.numcat)-length(Y2.numcat))))
    } else if (!is.null(Y2.con)) {
      Y2=Y2.con
      Y2i=Y2.con
    } else {
      Y2=Y2.cat
      Y2i=matrix(0,nrow(Y2.cat),(sum(Y2.numcat)-length(Y2.numcat)))
    } 
    if (any(is.na(Y2))) {
      if (ncol(Y2)==1) {
        miss.pat2<-matrix(c(0,1),2,1)
        n.patterns[2]<-2
      } else  {
        miss.pat2<-md.pattern.mice(Y2, plot=F)
        miss.pat2<-miss.pat2[,colnames(Y2)]
        n.patterns[2]<-nrow(miss.pat2)-1
      }
    } else {
      miss.pat2<-matrix(0,2,ncol(Y2))
      n.patterns[2]<-nrow(miss.pat2)-1
    }
    
    miss.pat.id2<-rep(0,nrow(Y2))
    for (i in 1:nrow(Y2)) {
      k <- 1
      flag <- 0
      while ((k <= n.patterns[2]) & (flag == 0)) {
        if (all(!is.na(Y2[i,])==miss.pat2[k,1:(ncol(miss.pat2))])) {
          miss.pat.id2[i] <- k
          flag <- 1
        } else {
          k <- k + 1
        }
      }
    }
    h=1
    if (isnullcat==0) {
      for (i in 1:length(Y.numcat)) {
        for (j in 1:nrow(Y)) {
          if (is.na(Y.cat[j,i])) {
            Yi[j,(ncolYcon+h):(ncolYcon+h+Y.numcat[i]-2)]=NA
          }
        } 
        h=h+Y.numcat[i]-1
      }
    }
    h=1
    if (isnullcat2==0) {
      for (i in 1:length(Y2.numcat)) {
        for (j in 1:nrow(Y2)) {
          if (is.na(Y2.cat[j,i])) {
            Y2i[j,(ncolY2con+h):(ncolY2con+h+Y2.numcat[i]-2)]=NA
          }
        } 
        h=h+Y2.numcat[i]-1
      }
    }
    
    if (output!=1) out.iter=nburn+2
    nimp=1
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+3)
    imp[1:nrow(Y),1:ncol(Y)]=Y
    imp[1:nrow(Y2),(ncol(Y)+1):(ncol(Y)+ncol(Y2))]=Y2
    imp[1:nrow(X), (ncol(Y)+ncol(Y2)+1):(ncol(Y)+ncol(Y2)+ncol(X))]=X
    imp[1:nrow(X2), (ncol(Y)+ncol(Y2)+ncol(X)+1):(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2))]=X2
    imp[1:nrow(Z), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+1):(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z))]=Z
    imp[1:nrow(clus), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+1)]=clus
    imp[1:nrow(X), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+2)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    Y2imp=Y2i
    Y2imp2=matrix(Y2imp, nrow(Y2imp),ncol(Y2imp))
    imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+ncol(Y2)+1):(ncol(Y)+ncol(Y2)+ncol(X))]=X
    imp[(nrow(X2)+1):(2*nrow(X2)), (ncol(Y)+ncol(Y2)+ncol(X)+1):(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2))]=X2
    imp[(nrow(Z)+1):(2*nrow(Z)), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+1):(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z))]=Z
    imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+1)]=clus
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+2)]=c(1:nrow(Y))
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+3)]=1  
    betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),nburn))
    beta2post<- array(0, dim=c(nrow(l2.beta.start),ncol(l2.beta.start),nburn))
    omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),nburn))
    upostall<-array(0, dim=c(nrow(u.start),ncol(u.start),nburn))
    covupost<- array(0, dim=c(nrow(l2cov.start),ncol(l2cov.start),nburn))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    if (!is.null(start.imp)) {
      start.imp<-as.matrix(start.imp)
      if ((nrow(start.imp)!=nrow(Yimp2))||(ncol(Yimp2)>ncol(start.imp))) {
        cat("start.imp dimensions incorrect. Not using start.imp as starting value for the imputed dataset.\n")
        start.imp=NULL
      } else {
        if ((nrow(start.imp)==nrow(Yimp2))&(ncol(Yimp2)<ncol(start.imp))) {
          Yimp2<-start.imp[,1:ncol(Yimp2)]
          cat("NOTE: start.imp has more columns than needed. Dropping unnecessary columns.\n")
        } else {
          Yimp2<-start.imp
        }
      }
    }
    if (is.null(start.imp)) {
      for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    }
    l2.meanobs<-colMeans(Y2i,na.rm=TRUE)
    if (!is.null(l2.start.imp)) {
      l2.start.imp<-as.matrix(l2.start.imp)
      if ((nrow(l2.start.imp)!=nrow(Y2imp2))||(ncol(Y2imp2)>ncol(l2.start.imp))) {
        cat("l2.start.imp dimensions incorrect. Not using l2.start.imp as starting value for the level 2 imputed dataset.\n")
        l2.start.imp=NULL
      } else {
        if ((nrow(l2.start.imp)==nrow(Y2imp2))&(ncol(Y2imp2)<ncol(l2.start.imp))) {
          Y2imp2<-l2.start.imp[,1:ncol(Y2imp2)]
          cat("NOTE: l2.start.imp has more columns than needed. Dropping unnecessary columns.\n")
        } else {
          Y2imp2<-l2.start.imp
        }
      }
    }
    if (is.null(l2.start.imp)) {
      for (i in 1:nrow(Y2i)) for (j in 1:ncol(Y2i)) if (is.na(Y2imp[i,j])) Y2imp2[i,j]=l2.meanobs[j]
    }
    .Call("jomo2comC", Y, Yimp, Yimp2, Y.cat, Y2, Y2imp,Y2imp2, Y2.cat, X, X2, Z, clus,betait,beta2it,uit,betapost,beta2post,upostall,covit,omegapost, covuit, covupost, nburn, l1cov.prior,l2cov.prior,Y.numcat,Y2.numcat, ncolYcon,ncolY2con, out.iter,1, miss.pat.id, n.patterns, miss.pat.id2, PACKAGE = "jomo")
    if (!is.null(Y.con)) {
      imp[(nrow(X)+1):(2*nrow(X)),1:ncol(Y.con)]=Yimp2[,1:ncol(Y.con)]
    }
    if (isnullcat==0) {
      imp[(nrow(X)+1):(2*nrow(X)),(ncolYcon+1):ncol(Y)]=Y.cat
    }
    if (!is.null(Y2.con)) {
      imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(Y2.con))]=Y2imp2[,1:ncol(Y2.con)]
    }
    if (isnullcat2==0) {
      imp[(nrow(X)+1):(2*nrow(X)),(ncolY2con+ncol(Y)+1):(ncol(Y)+ncol(Y2))]=Y2.cat
    }
    betapostmean<-apply(betapost, c(1,2), mean)
    beta2postmean<-apply(beta2post, c(1,2), mean)
    upostmean<-apply(upostall, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    covupostmean<-apply(covupost, c(1,2), mean)
    imp<-data.frame(imp)
    if (isnullcat==0) {
      for (i in 1:ncol(Y.cat)) {
        imp[,(ncolYcon+i)]<-as.factor(imp[,(ncolYcon+i)]) 
        levels(imp[,(ncolYcon+i)])<-previous_levels[[i]]
      }
    }
    if (isnullcat2==0) {
      for (i in 1:ncol(Y2.cat)) {
        imp[,(ncol(Y)+ncolY2con+i)]<-as.factor(imp[,(ncol(Y)+ncolY2con+i)]) 
        levels(imp[,(ncol(Y)+ncolY2con+i)])<-previous_levels2[[i]]
      }
    }
    imp[,(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+1)]<-factor(imp[,(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+1)])
    levels(imp[,(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z)+1)])<-previous_levels_clus
    clus<-factor(clus)
    levels(clus)<-previous_levels_clus
    if (ncolYcon>0) {
      for (j in 1:(ncolYcon)) {
        imp[,j]=as.numeric(imp[,j])
      }
    }
    if (ncolY2con>0) {
      for (j in 1:(ncolY2con)) {
        imp[,ncol(Y)+j]=as.numeric(imp[,ncol(Y)+j])
      }
    }
    for (j in (ncol(Y)+ncol(Y2)+1):(ncol(Y)+ncol(Y2)+ncol(X)+ncol(X2)+ncol(Z))) {
      imp[,j]=as.numeric(imp[,j])
    }
    if (isnullcat==0) {
      if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y.cat), sep = "")
    } else {
      colnamycat=NULL
      Y.cat=NULL
      Y.numcat=NULL
    }
    if (!is.null(Y.con)) {
      if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y.con), sep = "")
    } else {
      colnamycon=NULL
    }
    if (isnullcat2==0) {
      if (is.null(colnamy2cat)) colnamy2cat=paste("Y2cat", 1:ncol(Y2.cat), sep = "")
    } else {
      colnamy2cat=NULL
      Y2.cat=NULL
      Y2.numcat=NULL
    }
    if (!is.null(Y2.con)) {
      if (is.null(colnamy2con)) colnamy2con=paste("Y2con", 1:ncol(Y2.con), sep = "")
    } else {
      colnamy2con=NULL
    }
    if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    if (is.null(colnamx2)) colnamx2=paste("X2", 1:ncol(X2), sep = ".")
    colnames(imp)<-c(colnamycon,colnamycat,colnamy2con,colnamy2cat,colnamx,colnamx2,colnamz,"clus","id","Imputation")
    if (isnullcat==0) {
      cnycatcomp<-rep(NA,(sum(Y.numcat)-length(Y.numcat)))
      count=0
      for ( j in 1:ncol(Y.cat)) {
        for (k in 1:(Y.numcat[j]-1)) {
          cnycatcomp[count+k]<-paste(colnamycat[j],k,sep=".")
        }
        count=count+Y.numcat[j]-1
      }
      if (!is.null(Y.con)) {
        cnamycomp<-c(colnamycon,cnycatcomp)
      } else {
        cnamycomp<-c(cnycatcomp)
      }
      
    } else {
      cnamycomp<-c(colnamycon)
    }
    if (isnullcat2==0) {
      cny2catcomp<-rep(NA,(sum(Y2.numcat)-length(Y2.numcat)))
      count=0
      for ( j in 1:ncol(Y2.cat)) {
        for (k in 1:(Y2.numcat[j]-1)) {
          cny2catcomp[count+k]<-paste(colnamy2cat[j],k,sep=".")
        }
        count=count+Y2.numcat[j]-1
      }
      if (!is.null(Y2.con)) {
        cnamy2comp<-c(colnamy2con,cny2catcomp)
      } else {
        cnamy2comp<-c(cny2catcomp)
      }
      
    } else {
      cnamy2comp<-c(colnamy2con)
    }
    dimnames(betapost)[1] <- list(colnamx)
    dimnames(betapost)[2] <- list(cnamycomp)
    dimnames(beta2post)[1] <- list(colnamx2)
    dimnames(beta2post)[2] <- list(cnamy2comp)
    dimnames(omegapost)[1] <- list(cnamycomp)
    dimnames(omegapost)[2] <- list(cnamycomp)
    colnamcovu<-paste(cnamycomp,rep(colnamz,each=ncol(omegapost)),sep="*")
    colnamcovu<-c(colnamcovu,cnamy2comp)
    dimnames(covupost)[1] <- list(colnamcovu)
    dimnames(covupost)[2] <- list(colnamcovu)
    dimnames(upostall)[1]<-list(levels(clus))
    dimnames(upostall)[2]<-list(colnamcovu)
    dimnames(Yimp2)[2] <- list(cnamycomp)
    dimnames(Y2imp2)[2] <- list(cnamy2comp)
    betapostmean<-data.frame(apply(betapost, c(1,2), mean))
    beta2postmean<-data.frame(apply(beta2post, c(1,2), mean))
    upostmean<-data.frame(apply(upostall, c(1,2), mean))
    omegapostmean<-data.frame(apply(omegapost, c(1,2), mean))
    covupostmean<-data.frame(apply(covupost, c(1,2), mean))
    if (output==1) {
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(t(betapostmean))
      cat("\nThe posterior mean of the level 2 fixed effects estimates is:\n")
      print(t(beta2postmean))
      cat("\nThe posterior mean of the random effects estimates is:\n")
      print(upostmean)
      cat("\nThe posterior mean of the level 1 covariance matrix is:\n")
      print(omegapostmean)
      cat("\nThe posterior mean of the level 2 covariance matrix is:\n")
      print(covupostmean)
    }
    return(list("finimp"=imp,"collectbeta"=betapost,"collect.l2.beta"=beta2post,"collectomega"=omegapost,"collectu"=upostall, "collectcovu"=covupost, "finimp.latnorm" = Yimp2, "l2.finimp.latnorm" = Y2imp2))
  }
