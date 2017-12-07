jomo.coxph.MCMCchain <-
  function(formula, data,  beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, nburn=1000, start.imp=NULL, betaY.start=NULL, output=1, out.iter=10) {
    cat("This function is beta software. Use carefully and please report any bug to the package mantainer\n")
    stopifnot(is.data.frame(data))
    stopifnot(any(grepl("~",deparse(formula))))
    fit.cr<-coxph(formula,data=data, na.action = na.omit)
    if (is.null(betaY.start)) betaY.start<-as.numeric(coef(fit.cr))
    colnamysub<-all.vars(formula[[2]])
    Ysub<-as.matrix(data[,colnamysub])
    Ycov<-data.frame(mget(all.vars(formula[[3]]), envir =as.environment(data)))
    Y.con<-NULL
    Y.cat<-NULL
    Y.numcat<-NULL
    for (j in 1:ncol(Ycov)) {
      if (is.numeric(Ycov[,j])) {
        if (is.null(Y.con)) Y.con<-data.frame(Ycov[,j,drop=FALSE])
        else Y.con<-data.frame(Y.con,Ycov[,j,drop=FALSE])
      }
      if (is.factor(Ycov[,j])) {
        if (is.null(Y.cat)) Y.cat<-data.frame(Ycov[,j,drop=FALSE])
        else Y.cat<-data.frame(Y.cat,Ycov[,j,drop=FALSE])
        Y.numcat<-cbind(Y.numcat,nlevels(Ycov[,j]))
      }
      
    }
    terms.sub<-attr(terms(formula), "term.labels")
    split.terms<-strsplit(terms.sub,":")
    length.sub<-length(terms.sub)
    order.sub<-attr(terms(formula), "order")
    submod<-matrix(1,4,sum(order.sub))
    h<-1
    for ( j in 1:length.sub) {
      for ( k in 1:order.sub[j]) {
        current.term<-split.terms[[j]][k]
        current.term<-sub(".*I\\(","",current.term)
        current.term<-sub("\\)","",current.term)
        if (grepl("\\^",current.term)) {
          submod[3,h]<-as.integer(sub(".*\\^","",current.term))
          current.term<-sub("\\^.*","",current.term)
        } else {
          submod[3,h]<-1
        }
        if (length(which(colnames(Y.con)==current.term))==0) {
          submod[1,h]<-which(colnames(Y.cat)==current.term)
          submod[2,h]<-2
          submod[4,h]<-Y.numcat[submod[1,h]]-1
        } else {
          submod[1,h]<-which(colnames(Y.con)==current.term)
          submod[2,h]<-1
        }
        h<-h+1  
      }
    }
    X<-NULL
    if (is.null(X)) X=matrix(1,max(nrow(Y.cat),nrow(Y.con)),1)
    if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(max(0,ncol(Y.con))+max(0,(sum(Y.numcat)-length(Y.numcat)))))
    if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
    if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
    ncolYcon=max(0,ncol(Y.con))
    stopifnot(((!is.null(Y.con))||(!is.null(Y.cat)&!is.null(Y.numcat))))
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
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    if (!is.null(Y.con)) {
      stopifnot(nrow(Y.con)==nrow(X))
    }
    if (isnullcat==0) {
      stopifnot(nrow(Y.cat)==nrow(X))
    }
    stopifnot(nrow(beta.start)==ncol(X), ncol(beta.start)==(ncolYcon+max(0,(sum(Y.numcat)-length(Y.numcat)))))
    stopifnot(nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(beta.start))
    stopifnot(nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start))
    betait=matrix(0,nrow(beta.start),ncol(beta.start))
    for (i in 1:nrow(beta.start)) {
      for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
    }
    covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    for (i in 1:nrow(l1cov.start)) {
      for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
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
    colnamx<-colnames(X)
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"    
    if (!is.null(Y.con)&(isnullcat==0)) {
      Y=cbind(Y.con,Y.cat)
      Yi=cbind(Y.con, matrix(0,nrow(Y.con),(sum(Y.numcat)-length(Y.numcat))))
    } else if (!is.null(Y.con)) {
      Y=Y.con
      Yi=Y.con
    } else {
      Y=Y.cat
      Yi=matrix(0,nrow(Y.cat),(sum(Y.numcat)-length(Y.numcat)))
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
    
    if (output!=1) out.iter=nburn+2
    nimp=1
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+4)
    imp[1:nrow(Y),1:2]=Ysub
    imp[1:nrow(Y),3:(ncol(Y)+2)]=Y
    imp[1:nrow(X), (ncol(Y)+3):(ncol(Y)+ncol(X)+2)]=X
    imp[1:nrow(X), (ncol(Y)+ncol(X)+3)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+3):(ncol(Y)+ncol(X)+2)]=X
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+3)]=c(1:nrow(Y))
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+4)]=1  
    betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),(nburn)))
    betaYpost<- array(0, dim=c(1,length(betaY.start),(nburn)))
    omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),(nburn)))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    if (!is.null(start.imp)) {
      start.imp<-as.matrix(start.imp)
      if ((nrow(start.imp)!=nrow(Yimp2))||(ncol(Yimp2)!=ncol(start.imp))) {
        cat("start.imp dimensions incorrect. Not using start.imp as starting value for the imputed dataset.\n")
        start.imp=NULL
      } else {
        Yimp2<-start.imp
      }
    }
    if (is.null(start.imp)) {
        for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    }
    .Call("MCMCjomocox", Ysub, submod, order.sub, Y, Yimp, Yimp2, Y.cat, X, betaY.start,betaYpost, betait,betapost, covit,omegapost, nburn, l1cov.prior,Y.numcat,ncolYcon,out.iter, PACKAGE="jomo")
    #betapost[,,1]=bpost
    #upostall[,,1]=upost
    #omegapost[,,(1)]=opost
    #covupost[,,(1)]=cpost
    bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
    bYpost<-matrix(0,1,length(betaY.start))
    opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    imp[(nrow(Y)+1):(2*nrow(Y)),1:2]=Ysub
    if (!is.null(Y.con)) {
      imp[(nrow(Y)+1):(2*nrow(Y)),3:(ncol(Y.con)+2)]=Yimp2[,1:ncol(Y.con)]
    }
    if (isnullcat==0) {
      imp[(nrow(Y)+1):(2*nrow(Y)),(ncolYcon+3):(ncol(Y)+2)]=Y.cat
    }
    if (output==1) cat("First imputation registered.", "\n")
    
    betaYpostmean<-apply(betaYpost, c(1,2), mean)
    betapostmean<-apply(betapost, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    if (output==1) {
      cat("The posterior mean of the substantive model fixed effects estimates is:\n")
      print(betaYpostmean)
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(betapostmean)
      cat("The posterior mean of the level 1 covariance matrix is:\n")
      print(omegapostmean)
    }
    imp<-data.frame(imp)
    if (isnullcat==0) {
      for (i in 1:ncol(Y.cat)) {
        imp[,(ncolYcon+2+i)]<-as.factor(imp[,(ncolYcon+2+i)]) 
        levels(imp[,(ncolYcon+2+i)])<-previous_levels[[i]]
      }
    }
    
    if (ncolYcon>0) {
      for (j in 1:(ncolYcon)) {
        imp[,(2+j)]=as.numeric(imp[,(2+j)])
      }
    }
    
    for (j in (ncol(Y)+1):(ncol(Y)+ncol(X))) {
      imp[,(2+j)]=as.numeric(imp[,(2+j)])
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
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    if (is.null(colnamysub)) colnamysub="Ysub"
    colnames(imp)<-c(colnamysub,colnamycon,colnamycat,colnamx,"id","Imputation")
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
    dimnames(betapost)[1] <- list(colnamx)
    dimnames(betapost)[2] <- list(cnamycomp)
    dimnames(omegapost)[1] <- list(cnamycomp)
    dimnames(omegapost)[2] <- list(cnamycomp)
    dimnames(Yimp2)[2] <- list(cnamycomp)
    return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost, "finimp.latnorm" = Yimp2))  
    }
