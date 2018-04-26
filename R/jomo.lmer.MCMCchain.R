jomo.lmer.MCMCchain <-
  function(formula, data, level=rep(1,ncol(data)), beta.start=NULL, l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, a.start=NULL, a.prior=NULL, betaY.start=NULL, varY.start=NULL, covuY.start=NULL, uY.start=NULL, nburn=1000, meth="common", start.imp=NULL, start.imp.sub=NULL, l2.start.imp=NULL, output=1, out.iter=10) {
    cat("This function is beta software. Use carefully and please report any bug to the package mantainer\n")
    stopifnot(is.data.frame(data))
    stopifnot(any(grepl("~",deparse(formula))))
    fit.cr<-lmer(formula,data=data, na.action = na.omit)
    if (is.null(betaY.start)) betaY.start<-fixef(fit.cr)
    if (is.null(varY.start)) varY.start<-(summary(fit.cr)$sigma)^2
    varY.prior<-(summary(fit.cr)$sigma)^2
    colnamysub<-all.vars(formula[[2]])
    Ysub<-get(colnamysub,pos=data)
    Ycov<-data.frame(mget(all.vars(formula[[3]]), envir =as.environment(data)))
    level<-data.frame(matrix(level,1,ncol(data)))
    colnames(level)<-colnames(data)
    terms.sub<-attr(terms(formula), "term.labels")
    split.terms<-strsplit(terms.sub,":")
    length.sub<-length(terms.sub)
    order.sub<-attr(terms(formula), "order")
    submod<-matrix(1,4,sum(order.sub)-1)
    Y.con<-NULL
    Y.cat<-NULL
    Y.numcat<-NULL
    Y2.con<-NULL
    Y2.cat<-NULL
    Y2.numcat<-NULL
    Y2i<-NULL
    Y2imp2<-NULL
    for (j in 1:ncol(Ycov)) {
      if (as.data.frame(level[1, which(colnames(level)==colnames(Ycov)[j])])==1) {
        if (is.numeric(Ycov[,j])) {
          if (is.null(Y.con)) Y.con<-data.frame(Ycov[,j,drop=FALSE])
          else Y.con<-data.frame(Y.con,Ycov[,j,drop=FALSE])
        }
        if (is.factor(Ycov[,j])) {
          if (is.null(Y.cat)) Y.cat<-data.frame(Ycov[,j,drop=FALSE])
          else Y.cat<-data.frame(Y.cat,Ycov[,j,drop=FALSE])
          Y.numcat<-cbind(Y.numcat,nlevels(Ycov[,j]))
        }
      } else {
        if (is.numeric(Ycov[,j])) {
          if (is.null(Y2.con)) Y2.con<-data.frame(Ycov[,j,drop=FALSE])
          else Y2.con<-data.frame(Y2.con,Ycov[,j,drop=FALSE])
        }
        if (is.factor(Ycov[,j])) {
          if (is.null(Y2.cat)) Y2.cat<-data.frame(Ycov[,j,drop=FALSE])
          else Y2.cat<-data.frame(Y2.cat,Ycov[,j,drop=FALSE])
          Y2.numcat<-cbind(Y2.numcat,nlevels(Ycov[,j]))
        }
      }
      
      
    }
    h<-1
    for ( j in 1:length.sub) {
      for ( k in 1:order.sub[j]) {
        current.term<-split.terms[[j]][k]
        if (grepl("\\|", deparse(current.term))) {
          j.tbd<-j
          clus.name<-sub(".*\\|","",current.term)
          clus.name<-gsub(" ","",clus.name)
          current.term<-sub("\\|.*$","",current.term)
          random.terms<-strsplit(current.term,"+", fixed=T)
          length.ran<-length(random.terms[[1]])
          submod.ran<-matrix(1,3,length.ran)
          for (t in 1:length.ran) {
            ct<-gsub(" ","",random.terms[[1]][t])
            if (ct==1) {
              submod.ran[1:2,t]<-0
            } else if (length(which(colnames(Y.cat)==ct))!=0) {
              submod.ran[1,t]<-which(colnames(Y.cat)==ct)
              submod.ran[2,t]<-2
              submod.ran[3,t]<-Y.numcat[submod.ran[1,t]]-1
            } else if (length(which(colnames(Y.con)==ct))!=0) {
              submod.ran[1,t]<-which(colnames(Y.con)==ct)
              submod.ran[2,t]<-1
            }
          }
        } else {
          current.term<-sub(".*I\\(","",current.term)
          current.term<-sub("\\)","",current.term)
          if (grepl("\\^",current.term)) {
            submod[3,h]<-as.integer(sub(".*\\^","",current.term))
            current.term<-sub("\\^.*","",current.term)
          } else {
            submod[3,h]<-1
          }
          if (length(which(colnames(Y.cat)==current.term))!=0) {
            submod[1,h]<-which(colnames(Y.cat)==current.term)
            submod[2,h]<-2
            submod[4,h]<-Y.numcat[submod[1,h]]-1
          } else if (length(which(colnames(Y.con)==current.term))!=0) {
            submod[1,h]<-which(colnames(Y.con)==current.term)
            submod[2,h]<-1
          } else if (length(which(colnames(Y2.cat)==current.term))!=0) {
            submod[1,h]<-which(colnames(Y2.cat)==current.term)
            submod[2,h]<-4
            submod[4,h]<-Y2.numcat[submod[1,h]]-1
          } else if (length(which(colnames(Y2.con)==current.term))!=0) {
            submod[1,h]<-which(colnames(Y2.con)==current.term)
            submod[2,h]<-3
          }
          h<-h+1  
        }
        
      }
    }
    order.sub<-order.sub[-j.tbd]
    if (!is.null(Y.con)&sum((colnames(Y.con)==clus.name)==1)) Y.con<-data.frame(Y.con[,-which(colnames(Y.con)==clus.name), drop=FALSE])
    if (!is.null(Y.cat)&sum((colnames(Y.cat)==clus.name)==1)) Y.cat<-data.frame(Y.cat[,-which(colnames(Y.cat)==clus.name), drop=FALSE])
    Y.auxiliary<-data.frame(data[,-c(which(colnames(data)%in%colnames(Y.con)),which(colnames(data)%in%colnames(Y.cat)),which(colnames(data)%in%colnames(Y2.con)),which(colnames(data)%in%colnames(Y2.cat)),which(colnames(data)==clus.name),which(colnames(data)==colnamysub)), drop=FALSE])
    Y.aux.con<-NULL
    Y.aux.cat<-NULL
    Y.aux.numcat<-NULL
    Y2.aux.con<-NULL
    Y2.aux.cat<-NULL
    Y2.aux.numcat<-NULL
    if (ncol(Y.auxiliary)>0) {
      for (j in 1:ncol(Y.auxiliary)) {
        if (level[1, which(colnames(level)==colnames(Y.auxiliary)[j])]==1) {
          if (is.numeric(Y.auxiliary[,j])) {
            if (is.null(Y.aux.con)) Y.aux.con<-data.frame(Y.auxiliary[,j,drop=FALSE])
            else Y.aux.con<-data.frame(Y.aux.con,Y.auxiliary[,j,drop=FALSE])
          }
          if (is.factor(Y.auxiliary[,j])) {
            if (is.null(Y.aux.cat)) Y.aux.cat<-data.frame(Y.auxiliary[,j,drop=FALSE])
            else Y.aux.cat<-data.frame(Y.aux.cat,Y.auxiliary[,j,drop=FALSE])
            Y.aux.numcat<-cbind(Y.aux.numcat,nlevels(Y.auxiliary[,j]))
          }
        } else {
          if (is.numeric(Y.auxiliary[,j])) {
            if (is.null(Y2.aux.con)) Y2.aux.con<-data.frame(Y.auxiliary[,j,drop=FALSE])
            else Y2.aux.con<-data.frame(Y2.aux.con,Y.auxiliary[,j,drop=FALSE])
          }
          if (is.factor(Y.auxiliary[,j])) {
            if (is.null(Y2.aux.cat)) Y2.aux.cat<-data.frame(Y.auxiliary[,j,drop=FALSE])
            else Y2.aux.cat<-data.frame(Y2.aux.cat,Y.auxiliary[,j,drop=FALSE])
            Y2.aux.numcat<-cbind(Y2.aux.numcat,nlevels(Y.auxiliary[,j]))
          }
        }
        
        
      }
    }
    X=matrix(1,max(nrow(Y.cat),nrow(Y.con)),1)
    if (!is.null(Y2.con)|!is.null(Y2.cat)|!is.null(Y2.aux.con)|!is.null(Y2.aux.cat)) {
      #cat("Level 2 variables must be fully observed for valid inference. \n")
      X2=matrix(1,max(nrow(Y2.cat),nrow(Y2.con),nrow(Y2.aux.cat),nrow(Y2.aux.con) ),1)
      if (is.null(l2.beta.start)) l2.beta.start=matrix(0,ncol(X2),(max(as.numeric(!is.null(Y2.con)),ncol(Y2.con))+max(0,(sum(Y2.numcat)-length(Y2.numcat)))+max(as.numeric(!is.null(Y2.aux.con)),ncol(Y2.aux.con))+max(0,(sum(Y2.aux.numcat)-length(Y2.aux.numcat)))))
    }
    Z=matrix(1,nrow(X),1)
    if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(max(as.numeric(!is.null(Y.con)),ncol(Y.con))+max(0,(sum(Y.numcat)-length(Y.numcat)))+max(as.numeric(!is.null(Y.aux.con)),ncol(Y.aux.con))+max(0,(sum(Y.aux.numcat)-length(Y.aux.numcat)))))
    
    clus<-factor(data[,clus.name])
    previous_levels_clus<-levels(clus)
    levels(clus)<-0:(nlevels(clus)-1)
    if (is.null(uY.start)) uY.start<-matrix(0,nlevels(clus),ncol(VarCorr(fit.cr)[[1]]))
    if (is.null(l1cov.start)) {
      if (meth=="common") {
        l1cov.start=diag(1,ncol(beta.start))
      } else {
        l1cov.start=matrix(diag(1,ncol(beta.start)),nlevels(clus)*ncol(beta.start),ncol(beta.start),2)
      }
    }
    if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
    diagVar<-as.data.frame(VarCorr(fit.cr))[1:ncol(VarCorr(fit.cr)[[1]]),4]
    if (is.null(covuY.start)) covuY.start<-VarCorr(fit.cr)[[1]][1:ncol(VarCorr(fit.cr)[[1]]),1:ncol(VarCorr(fit.cr)[[1]])]
    covuY.prior<-VarCorr(fit.cr)[[1]][1:ncol(VarCorr(fit.cr)[[1]]),1:ncol(VarCorr(fit.cr)[[1]])]
    if (kappa(covuY.start)>10^8) {
      covuY.prior<-diag(1, ncol(VarCorr(fit.cr)[[1]]))
      covuY.start<-diag(1, ncol(VarCorr(fit.cr)[[1]]))
    }
    ncolYcon<-rep(NA,4)
    ncolY2con<-rep(NA,4)
    ncolYcon[1]=max(as.numeric(!is.null(Y.con)),ncol(Y.con))+max(as.numeric(!is.null(Y.aux.con)),ncol(Y.aux.con))
    ncolY2con[1]=max(as.numeric(!is.null(Y2.con)),ncol(Y2.con))+max(as.numeric(!is.null(Y2.aux.con)),ncol(Y2.aux.con))
    ncolYcon[2]=max(as.numeric(!is.null(Y.con)),ncol(Y.con))
    ncolY2con[2]=max(as.numeric(!is.null(Y2.con)),ncol(Y2.con))
    ncolYcon[3]=ncolYcon[1]+max(0,(sum(Y.numcat)-length(Y.numcat)))
    ncolY2con[3]=ncolY2con[1]+max(0,(sum(Y2.numcat)-length(Y2.numcat)))
    ncolYcon[4]=max(0,ncol(Y.cat))
    ncolY2con[4]=max(0,ncol(Y2.cat))
    stopifnot(((!is.null(Y.con))||(!is.null(Y.cat)&!is.null(Y.numcat)))||((!is.null(Y2.con))||(!is.null(Y2.cat)&!is.null(Y2.numcat))))
    if (is.null(u.start)) u.start = matrix(0, nlevels(clus), ncol(Z)*(ncolYcon[1]+max(0,(sum(Y.numcat)-length(Y.numcat)))+max(0,(sum(Y.aux.numcat)-length(Y.aux.numcat))))+(ncolY2con[1]+max(0,(sum(Y2.numcat)-length(Y2.numcat)))+max(0,(sum(Y2.aux.numcat)-length(Y2.aux.numcat)))))
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
    }
    if (!is.null(Y.aux.cat)) {
      isnullcataux=0
      previous_levelsaux<-list()
      Y.aux.cat<-data.frame(Y.aux.cat)
      for (i in 1:ncol(Y.aux.cat)) {
        Y.aux.cat[,i]<-factor(Y.aux.cat[,i])
        previous_levelsaux[[i]]<-levels(Y.aux.cat[,i])
        levels(Y.aux.cat[,i])<-1:nlevels(Y.aux.cat[,i])
      }
    } else {
      isnullcataux=1
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
    }
    if (!is.null(Y2.aux.cat)) {
      isnullcat2aux=0
      previous_levels2aux<-list()
      Y2.aux.cat<-data.frame(Y2.aux.cat)
      for (i in 1:ncol(Y2.aux.cat)) {
        Y2.aux.cat[,i]<-factor(Y2.aux.cat[,i])
        previous_levels2aux[[i]]<-levels(Y2.aux.cat[,i])
        levels(Y2.aux.cat[,i])<-1:nlevels(Y2.aux.cat[,i])
      }
    } else {
      isnullcat2aux=1
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
    stopifnot(nrow(beta.start)==ncol(X), ncol(beta.start)==(ncolYcon[1]+max(0,(sum(Y.numcat)-length(Y.numcat)))+max(0,(sum(Y.aux.numcat)-length(Y.aux.numcat)))))
    if (!is.null(Y2.con)||isnullcat2==0||!is.null(Y2.aux.con)||isnullcat2aux==0) stopifnot(nrow(l2.beta.start)==ncol(X2), ncol(l2.beta.start)==(ncolY2con[3]+max(0,(sum(Y2.aux.numcat)-length(Y2.aux.numcat)))))
    stopifnot(ncol(u.start)==ncol(Z)*(ncolYcon[1]+max(0,(sum(Y.numcat)-length(Y.numcat)))+max(0,(sum(Y.aux.numcat)-length(Y.aux.numcat))))+(ncolY2con[1]+max(0,(sum(Y2.numcat)-length(Y2.numcat)))+max(0,(sum(Y2.aux.numcat)-length(Y2.aux.numcat)))))
    if (meth=="common") stopifnot(nrow(l1cov.start)==ncol(l1cov.start),nrow(l1cov.prior)==nrow(l1cov.start),nrow(l1cov.start)==ncol(beta.start))
    stopifnot(nrow(l1cov.prior)==ncol(l1cov.prior))
    stopifnot(ncol(l2cov.start)==ncol(u.start), nrow(l2cov.prior)==nrow(l2cov.start), nrow(l2cov.prior)==ncol(l2cov.prior))
    betait=matrix(0,nrow(beta.start),ncol(beta.start))
    for (i in 1:nrow(beta.start)) {
      for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
    }
    if (!is.null(Y2.con)||isnullcat2==0||!is.null(Y2.aux.con)||isnullcat2aux==0) {
      beta2it=matrix(0,nrow(l2.beta.start),ncol(l2.beta.start))
      for (i in 1:nrow(l2.beta.start)) {
        for (j in 1:ncol(l2.beta.start)) beta2it[i,j]=l2.beta.start[i,j]
      }
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
    if (!is.null(Y.aux.con)) {
      colnamyauxcon<-colnames(Y.aux.con)
      Y.aux.con<-data.matrix(Y.aux.con)
      storage.mode(Y.aux.con) <- "numeric"  
    }
    if (isnullcataux==0) {
      colnamyauxcat<-colnames(Y.aux.cat)
      Y.aux.cat<-data.matrix(Y.aux.cat)
      storage.mode(Y.aux.cat) <- "numeric"  
    }
    if (!is.null(Y2.aux.con)) {
      colnamy2auxcon<-colnames(Y2.aux.con)
      Y2.aux.con<-data.matrix(Y2.aux.con)
      storage.mode(Y2.aux.con) <- "numeric"  
    }
    if (isnullcat2aux==0) {
      colnamy2auxcat<-colnames(Y2.aux.cat)
      Y2.aux.cat<-data.matrix(Y2.aux.cat)
      storage.mode(Y2.aux.cat) <- "numeric"  
    }
    Y.cat.tot<-cbind(Y.cat,Y.aux.cat)
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"    
    Z<-data.matrix(Z)
    storage.mode(Z) <- "numeric" 
    if (!is.null(Y2.con)||isnullcat2==0||!is.null(Y2.aux.con)||isnullcat2aux==0) {
      colnamx2<-colnames(X2)
      X2<-data.matrix(X2)
      storage.mode(X2) <- "numeric"
    }
    clus <- matrix(as.integer(levels(clus))[clus], ncol=1)
    Y=cbind(Y.con,Y.aux.con,Y.cat, Y.aux.cat)
    Yi=cbind(Y.con, Y.aux.con, switch(is.null(Y.cat)+1, matrix(0,nrow(Y),(sum(Y.numcat)-length(Y.numcat))), NULL), switch(is.null(Y.aux.cat)+1, matrix(0,nrow(Y.aux.cat),(sum(Y.aux.numcat)-length(Y.aux.numcat))), NULL))
    Y2=cbind(Y2.con,Y2.aux.con,Y2.cat, Y2.aux.cat)
    Y2i=cbind(Y2.con, Y2.aux.con, switch(is.null(Y2.cat)+1, matrix(0,nrow(Y2),(sum(Y2.numcat)-length(Y2.numcat))), NULL), switch(is.null(Y2.aux.cat)+1, matrix(0,nrow(Y2.aux.cat),(sum(Y2.aux.numcat)-length(Y2.aux.numcat))), NULL))
    h=1
    if (isnullcat==0) {
      for (i in 1:length(Y.numcat)) {
        for (j in 1:nrow(Y)) {
          if (is.na(Y.cat[j,i])) {
            Yi[j,(ncolYcon[1]+h):(ncolYcon[1]+h+Y.numcat[i]-2)]=NA
          }
        } 
        h=h+Y.numcat[i]-1
      }
    }
    if (isnullcataux==0) {
      for (i in 1:length(Y.aux.numcat)) {
        for (j in 1:nrow(Y)) {
          if (is.na(Y.aux.cat[j,i])) {
            Yi[j,(ncolYcon[1]+h):(ncolYcon[1]+h+Y.aux.numcat[i]-2)]=NA
          }
        } 
        h=h+Y.aux.numcat[i]-1
      }
    }
    h=1
    if (isnullcat2==0) {
      for (i in 1:length(Y2.numcat)) {
        for (j in 1:nrow(Y2)) {
          if (is.na(Y2.cat[j,i])) {
            Y2i[j,(ncolY2con[1]+h):(ncolY2con[1]+h+Y2.numcat[i]-2)]=NA
          }
        } 
        h=h+Y2.numcat[i]-1
      }
    }
    if (isnullcat2aux==0) {
      for (i in 1:length(Y2.aux.numcat)) {
        for (j in 1:nrow(Y2)) {
          if (is.na(Y2.aux.cat[j,i])) {
            Y2i[j,(ncolY2con[1]+h):(ncolY2con[1]+h+Y2.aux.numcat[i]-2)]=NA
          }
        } 
        h=h+Y2.aux.numcat[i]-1
      }
    }
    if (isnullcat==0||isnullcataux==0) {
      Y.cat.tot<-cbind(Y.cat,Y.aux.cat)
      Y.numcat.tot<-c(Y.numcat, Y.aux.numcat)
    } else {
      Y.cat.tot=-999
      Y.numcat.tot=-999
    }
    if (isnullcat2==0||isnullcat2aux==0) {
      Y2.cat.tot<-cbind(Y2.cat,Y2.aux.cat)
      Y2.numcat.tot<-c(Y2.numcat, Y2.aux.numcat)
    } else {
      Y2.cat.tot=-999
      Y2.numcat.tot=-999
    }
    ncY2<-max(0,ncol(Y2))
    Ysubimp<-Ysub
    if (is.null(a.start)) a.start=50+ncol(Y)
    if (is.null(a.prior)) a.prior=a.start
    
    if (output!=1) out.iter=nburn+2
    nimp=1
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncY2+4)
    imp[1:nrow(Y),1]=Ysub
    imp[1:nrow(Y),2:(1+ncol(Y))]=Y
    if (!is.null(Y2)) imp[1:nrow(Y2),(ncol(Y)+2):(ncol(Y)+ncY2+1)]=Y2
    imp[1:nrow(clus), (ncol(Y)+ncY2+2)]=clus
    imp[1:nrow(X), (ncol(Y)+ncY2+3)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    if (!is.null(Y2)) {
      Y2imp=Y2i
      Y2imp2=matrix(Y2imp, nrow(Y2imp),ncol(Y2imp))
    }
    imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncY2+2)]=clus
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncY2+3)]=c(1:nrow(Y))
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncY2+4)]=1  
    betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),(nburn)))
    betaYpost<- array(0, dim=c(1,length(betaY.start),(nburn)))
    if (!is.null(Y2)) {
      beta2post<- array(0, dim=c(nrow(l2.beta.start),ncol(l2.beta.start),(nburn)))
    }
    upostall<-array(0, dim=c(nrow(u.start),ncol(u.start),(nburn)))
    uYpostall<-array(0, dim=c(nrow(uY.start),ncol(uY.start),(nburn)))
    omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),(nburn)))
    varYpost<-rep(0,(nburn))
    covupost<- array(0, dim=c(nrow(l2cov.start),ncol(l2cov.start),(nburn)))
    covuYpost<-array(0, dim=c(nrow(as.matrix(covuY.start)),ncol(as.matrix(covuY.start)),(nburn)))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    if (!is.null(Y2)) {
      meanobs2<-colMeans(Y2i,na.rm=TRUE)
    }
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
    if (!is.null(l2.start.imp)) {
      l2.start.imp<-as.matrix(l2.start.imp)
      if ((nrow(l2.start.imp)!=nrow(Y2imp2))||(ncol(Y2imp2)!=ncol(l2.start.imp))) {
        cat("l2.start.imp dimensions incorrect. Not using l2.start.imp as starting value for the level 2 imputed dataset.\n")
        l2.start.imp=NULL
      } else {
        Y2imp2<-l2.start.imp
      }
    }
    if (!is.null(Y2i)&is.null(l2.start.imp)) {
      for (i in 1:nrow(Y2i)) for (j in 1:ncol(Y2i)) if (is.na(Y2imp[i,j])) Y2imp2[i,j]=meanobs2[j]
    }
    if (!is.null(start.imp.sub)) {
      start.imp<-as.matrix(start.imp)
      if (!is.vector(start.imp.sub)) {
        cat("start.imp.sub must be a vector. Not using start.imp as starting value for the imputed dataset.\n")
        start.imp.sub=NULL
      } else {
        Ysubimp=start.imp.sub
      }
    }
    if (is.null(start.imp.sub)) {
      for (i in 1:length(Ysubimp)) if (is.na(Ysubimp[i])) Ysubimp[i]=mean(Ysubimp, na.rm = TRUE)
    }
    if (!is.null(Y2)) {
      if (meth=="common") {
        .Call("MCMCjomo2lmer", Ysub, Ysubimp, submod, order.sub, submod.ran, Y, Yimp, Yimp2, Y.cat.tot, Y2, Y2imp, Y2imp2, Y2.cat.tot, X, X2, Z, clus,betaY.start,betaYpost, betait,beta2it,uit,uY.start,betapost, upostall, uYpostall, beta2post, varY.start, varYpost, covit,omegapost, covuY.start, covuYpost, covuit, covupost, nburn, varY.prior, covuY.prior, l1cov.prior,l2cov.prior,Y.numcat.tot, Y2.numcat.tot, ncolYcon,ncolY2con, out.iter, PACKAGE = "jomo")
        } else {
        .Call("MCMCjomo2lmerhr", Ysub, Ysubimp, submod, order.sub, submod.ran, Y, Yimp, Yimp2, Y.cat.tot, Y2, Y2imp, Y2imp2, Y2.cat.tot, X, X2, Z, clus,betaY.start,betaYpost, betait,beta2it,uit,uY.start,betapost, upostall, uYpostall, beta2post, varY.start, varYpost, covit,omegapost, covuY.start, covuYpost, covuit, covupost, nburn, varY.prior, covuY.prior, l1cov.prior,l2cov.prior,Y.numcat.tot, Y2.numcat.tot, ncolYcon,ncolY2con,a.start, a.prior, out.iter, PACKAGE = "jomo")
      }
    } else {
      if (meth=="common") {
        .Call("MCMCjomo1lmer", Ysub, Ysubimp, submod, order.sub, submod.ran, Y, Yimp, Yimp2, Y.cat.tot, X, Z, clus,betaY.start,betaYpost, betait,uit,uY.start,betapost, upostall, uYpostall, varY.start, varYpost, covit,omegapost, covuY.start, covuYpost, covuit, covupost, nburn, varY.prior, covuY.prior, l1cov.prior,l2cov.prior,Y.numcat.tot, ncolYcon,out.iter, PACKAGE = "jomo")

        } else {
        .Call("MCMCjomo1lmerhr", Ysub, Ysubimp, submod, order.sub, submod.ran, Y, Yimp, Yimp2, Y.cat.tot, X, Z, clus,betaY.start,betaYpost, betait,uit,uY.start,betapost, upostall, uYpostall, varY.start, varYpost, covit,omegapost, covuY.start, covuYpost, covuit, covupost, nburn, varY.prior, covuY.prior, l1cov.prior,l2cov.prior,Y.numcat.tot, ncolYcon,a.start, a.prior, out.iter, PACKAGE = "jomo")
        }
    }
    imp[(nrow(Y)+1):(2*nrow(Y)),1]=Ysubimp
    if (!is.null(Y.con)|!is.null(Y.aux.con)) {
      imp[(nrow(Y)+1):(2*nrow(Y)),2:(1+max(0,ncol(Y.con))+max(0,ncol(Y.aux.con)))]=Yimp2[,1:(max(0,ncol(Y.con))+max(0,ncol(Y.aux.con)))]
    }
    if (isnullcat==0|isnullcataux==0) {
      imp[(nrow(Y)+1):(2*nrow(Y)),(ncolYcon[1]+2):(1+ncol(Y))]=Y.cat.tot
    }
    if (!is.null(Y2.con)|!is.null(Y2.aux.con)) {
      imp[(nrow(Y2)+1):(2*nrow(Y2)),(ncol(Y)+2):(1+ncol(Y)+max(0,ncol(Y2.con))+max(0,ncol(Y2.aux.con)))]=Y2imp2[,1:(max(0,ncol(Y2.con))+max(0,ncol(Y2.aux.con)))]
    }
    if (isnullcat2==0|isnullcat2aux==0) {
      imp[(nrow(Y2)+1):(2*nrow(Y2)),(ncolY2con[1]+ncol(Y)+2):(ncol(Y)+ncY2+1)]=Y2.cat
    }
    if (output==1) cat("First imputation registered.", "\n")
    betaYpostmean<-apply(betaYpost, c(1,2), mean)
    varYpostmean<-mean(varYpost)
    covuYpostmean<-apply(covuYpost, c(1,2), mean)
    uYpostmean<-apply(uYpostall, c(1,2), mean)
    betapostmean<-apply(betapost, c(1,2), mean)
    if (!is.null(Y2)) beta2postmean<-apply(beta2post, c(1,2), mean)
    upostmean<-apply(upostall, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    covupostmean<-apply(covupost, c(1,2), mean)
    if (output==1) {
      cat("The posterior mean of the substantive model fixed effects estimates is:\n")
      print(betaYpostmean)
      cat("The posterior mean of the substantive model residual variance is:\n")
      print(varYpostmean)
      cat("The posterior mean of the substantive model random effects covariance matrix is:\n")
      print(covuYpostmean)
      cat("The posterior mean of the substantive model random effects estimates is:\n")
      print(uYpostmean)
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(betapostmean)
      if (!is.null(Y2)) {
        cat("The posterior mean of the level 2 fixed effects estimates is:\n")
        print(beta2postmean)
      }
      cat("The posterior mean of the random effects estimates is:\n")
      print(upostmean)
      cat("The posterior mean of the level 1 covariance matrix is:\n")
      print(omegapostmean)
      cat("The posterior mean of the level 2 covariance matrix is:\n")
      print(covupostmean)
    }
    imp<-data.frame(imp)
    if (isnullcat==0) {
      for (i in 1:ncol(Y.cat)) {
        imp[,(1+ncolYcon[1]+i)]<-as.factor(imp[,(1+ncolYcon[1]+i)]) 
        levels(imp[,(1+ncolYcon[1]+i)])<-previous_levels[[i]]
      }
    }
    if (isnullcataux==0) {
      for (i in 1:ncol(Y.aux.cat)) {
        imp[,(1+ncolYcon[1]+ncolYcon[4]+i)]<-as.factor(imp[,(1+ncolYcon[1]+ncolYcon[4]+i)]) 
        levels(imp[,(1+ncolYcon[1]+ncolYcon[4]+i)])<-previous_levelsaux[[i]]
      }
    }
    if (isnullcat2==0) {
      for (i in 1:ncol(Y2.cat)) {
        imp[,(1+ncol(Y)+ncolY2con[1]+i)]<-as.factor(imp[,(1+ncol(Y)+ncolY2con[1]+i)]) 
        levels(imp[,(1+ncol(Y)+ncolY2con[1]+i)])<-previous_levels2[[i]]
      }
    }
    if (isnullcat2aux==0) {
      for (i in 1:ncol(Y2.aux.cat)) {
        imp[,(1+ncol(Y)+ncolY2con[1]+ncolY2con[4]+i)]<-as.factor(imp[,(1+ncol(Y)+ncolY2con[1]+ncolY2con[4]+i)]) 
        levels(imp[,(1+ncol(Y)+ncolY2con[1]+ncolY2con[4]+i)])<-previous_levels2aux[[i]]
      }
    }
    imp[,(ncol(Y)+ncY2+2)]<-factor(imp[,(ncol(Y)+ncY2+2)])
    levels(imp[,(ncol(Y)+ncY2+2)])<-previous_levels_clus
    clus<-factor(clus)
    levels(clus)<-previous_levels_clus
    if (ncolYcon[1]>0) {
      for (j in 1:(ncolYcon[1])) {
        imp[,j+1]=as.numeric(imp[,j+1])
      }
    }
    if (ncolY2con[1]>0) {
      for (j in 1:(ncolY2con[1])) {
        imp[,ncol(Y)+j+1]=as.numeric(imp[,ncol(Y)+j+1])
      }
    }
    if (isnullcat==0) {
      if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y.cat), sep = "")
    } else {
      colnamycat=NULL
      Y.cat=NULL
      Y.numcat=NULL
    }
    if (isnullcataux==0) {
      if (is.null(colnamyauxcat)) colnamyauxcat=paste("Ycat.aux", 1:ncol(Y.aux.cat), sep = "")
    } else {
      colnamyauxcat=NULL
      Y.aux.cat=NULL
      Y.aux.numcat=NULL
    }
    if (!is.null(Y.con)) {
      if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y.con), sep = "")
    } else {
      colnamycon=NULL
    }
    if (!is.null(Y.aux.con)) {
      if (is.null(colnamyauxcon)) colnamyauxcon=paste("Ycon.aux", 1:ncol(Y.aux.con), sep = "")
    } else {
      colnamyauxcon=NULL
    }
    if (isnullcataux==0) {
      if (is.null(colnamyauxcat)) colnamyauxcat=paste("Ycat.aux", 1:ncol(Y.aux.cat), sep = "")
    } else {
      colnamyauxcat=NULL
      Y.aux.cat=NULL
      Y.aux.numcat=NULL
    }
    if (isnullcat2==0) {
      if (is.null(colnamy2cat)) colnamy2cat=paste("Y2cat", 1:ncol(Y2.cat), sep = "")
    } else {
      colnamy2cat=NULL
      Y2.cat=NULL
      Y2.numcat=NULL
    }
    if (isnullcat2aux==0) {
      if (is.null(colnamy2auxcat)) colnamy2auxcat=paste("Y2cat.aux", 1:ncol(Y2.aux.cat), sep = "")
    } else {
      colnamy2auxcat=NULL
      Y2.aux.cat=NULL
      Y2.aux.numcat=NULL
    }
    if (!is.null(Y2.con)) {
      if (is.null(colnamy2con)) colnamy2con=paste("Y2con", 1:ncol(Y2.con), sep = "")
    } else {
      colnamy2con=NULL
    }
    if (!is.null(Y2.aux.con)) {
      if (is.null(colnamy2auxcon)) colnamy2auxcon=paste("Y2con.aux", 1:ncol(Y2.aux.con), sep = "")
    } else {
      colnamy2auxcon=NULL
    }
    if (is.null(colnamysub)) colnamysub="Ysub"
    colnames(imp)<-c(colnamysub,colnamycon,colnamyauxcon,colnamycat,colnamyauxcat,colnamy2con,colnamy2auxcon,colnamy2cat,colnamy2auxcat,"clus","id","Imputation")
    if (isnullcat==0) {
      cnycatcomp<-rep(NA,(sum(Y.numcat)-length(Y.numcat)))
      count=0
      for ( j in 1:ncol(Y.cat)) {
        for (k in 1:(Y.numcat[j]-1)) {
          cnycatcomp[count+k]<-paste(colnamycat[j],k,sep=".")
        }
        count=count+Y.numcat[j]-1
      }
      
    } else {
      cnycatcomp<-NULL
    }
    if (isnullcataux==0) {
      cnyauxcatcomp<-rep(NA,(sum(Y.aux.numcat)-length(Y.aux.numcat)))
      count=0
      for ( j in 1:ncol(Y.aux.cat)) {
        for (k in 1:(Y.aux.numcat[j]-1)) {
          cnyauxcatcomp[count+k]<-paste(colnamyauxcat[j],k,sep=".")
        }
        count=count+Y.aux.numcat[j]-1
      }
      
    } else {
      cnyauxcatcomp<-NULL
    }
    cnamycomp<-c(colnamycon,colnamyauxcon,cnycatcomp,cnyauxcatcomp)
    if (isnullcat2==0) {
      cny2catcomp<-rep(NA,(sum(Y2.numcat)-length(Y2.numcat)))
      count=0
      for ( j in 1:ncol(Y2.cat)) {
        for (k in 1:(Y2.numcat[j]-1)) {
          cny2catcomp[count+k]<-paste(colnamy2cat[j],k,sep=".")
        }
        count=count+Y2.numcat[j]-1
      }
      
    } else {
      cny2catcomp<-NULL
    }
    if (isnullcat2aux==0) {
      cny2auxcatcomp<-rep(NA,(sum(Y2.aux.numcat)-length(Y2.aux.numcat)))
      count=0
      for ( j in 1:ncol(Y2.aux.cat)) {
        for (k in 1:(Y2.aux.numcat[j]-1)) {
          cny2auxcatcomp[count+k]<-paste(colnamy2auxcat[j],k,sep=".")
        }
        count=count+Y2.aux.numcat[j]-1
      }
      
    } else {
      cny2auxcatcomp<-NULL
    }
    cnamy2comp<-c(colnamy2con,colnamy2auxcon,cny2catcomp,cny2auxcatcomp)
    dimnames(betapost)[1] <- list(colnamx)
    dimnames(betapost)[2] <- list(cnamycomp)
    if (!is.null(Y2)) {
      dimnames(beta2post)[1] <- list(colnamx2)
      dimnames(beta2post)[2] <- list(cnamy2comp)
    }
    if (meth=="random") {
      dimnames(omegapost)[1] <- list(paste(rep(cnamycomp, nrow(u.start)),rep(previous_levels_clus,each=ncol(omegapost)))) 
    } else  {
      dimnames(omegapost)[1] <- list(cnamycomp)
    } 
    dimnames(omegapost)[2] <- list(cnamycomp)
    dimnames(Yimp2)[2] <- list(cnamycomp)

    return(list("finimp"=imp,"collectbeta"=betapost,"collectu"=upostall,"collectomega"=omegapost, "collectcovu"=covupost, "finimp.latnorm" = Yimp2, "l2.finimp.latnorm" = Y2imp2, "collectbetaY"=betaYpost, "collectvarY"=varYpost, "collectcovuY"=covuYpost, "collectuY"=uYpostall))
  }
