jomo1ranmixhr <-
function(Y_con, Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), Z=matrix(1,nrow(Y_cat),1), clus, betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), up=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), covp=matrix(diag(1,ncol(betap)),ncol(betap)*nrow(unique(clus)),ncol(betap),2), covu=diag(1,ncol(up)), Sp=diag(1,ncol(betap)), Sup=diag(1,ncol(up)), nburn=100, nbetween=100, nimp=5,a=ncol(betap)) {
  stopifnot(nrow(Y_con)==nrow(clus),nrow(Y_con)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))),nrow(covp)==nrow(up)*ncol(covp), nrow(covp)==nrow(up)*ncol(betap), nrow(Sp)==ncol(Sp),nrow(covp)==nrow(up)*nrow(Sp),nrow(Z)==nrow(Y_con), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
  rngflag=0
  colnamycon<-colnames(Y_con)
  colnamycat<-colnames(Y_cat)
  colnamx<-colnames(X)
  colnamz<-colnames(Z)
  Y_con<-as.matrix(Y_con,nrow(Y_con),ncol(Y_con))
  Y_cat<-as.matrix(Y_cat,nrow(Y_cat),ncol(Y_cat))
  X<-as.matrix(X,nrow(X),ncol(X))
  Z<-as.matrix(Z,nrow(Z),ncol(Z))
  clus<-as.matrix(clus,nrow(clus),ncol(clus))
  Y=cbind(Y_con,Y_cat)
  Yi=cbind(Y_con, matrix(0,nrow(Y_con),(sum(Y_numcat)-length(Y_numcat))))
  h=1
  for (i in 1:length(Y_numcat)) {
    for (j in 1:nrow(Y)) {
      if (is.na(Y_cat[j,i])) {
        Yi[j,(ncol(Y_con)+h):(ncol(Y_con)+h+Y_numcat[i]-2)]=NA
      }
    } 
    h=h+Y_numcat[i]-1
  }
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+ncol(Z)+3)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(Z), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[1:nrow(clus), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[1:nrow(X), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  Yimp=Yi
  Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(Z)+1):(2*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=1  
  betapost<- array(0, dim=c(nrow(betap),ncol(betap),(nimp-1)))
  bpost<-matrix(0,nrow(betap),ncol(betap))
  upost<-matrix(0,nrow(up),ncol(up))
  upostall<-array(0, dim=c(nrow(up),ncol(up),(nimp-1)))
  omegapost<- array(0, dim=c(nrow(covp),ncol(covp),(nimp-1)))
  opost<-matrix(0,nrow(covp),ncol(covp))
  covupost<- array(0, dim=c(nrow(covu),ncol(covu),(nimp-1)))
  cpost<-matrix(0,nrow(covu),ncol(covu))
  meanobs<-colMeans(Yi,na.rm=TRUE)
  for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
  .Call("jomo1ranmixhr", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betap,up,bpost,upost,covp,opost, covu,cpost,nburn, Sp,Sup,Y_numcat, ncol(Y_con),a,rngflag, PACKAGE = "jomo")
  #betapost[,,1]=bpost
  #upostall[,,1]=upost
  #omegapost[,,(1)]=opost
  #covupost[,,(1)]=cpost
  bpost<-matrix(0,nrow(betap),ncol(betap))
  opost<-matrix(0,nrow(covp),ncol(covp))
  upost<-matrix(0,nrow(up),ncol(up))
  cpost<-matrix(0,nrow(covu),ncol(covu))
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y_con)]=Yimp2[,1:ncol(Y_con)]
  imp[(nrow(Y)+1):(2*nrow(Y)),(ncol(Y_con)+1):ncol(Y)]=Y_cat
  cat("First imputation registered.", "\n")
  for (i in 2:nimp) {
    #Yimp2=matrix(0, nrow(Yimp),ncol(Yimp))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
    imp[(i*nrow(clus)+1):((i+1)*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=i      
    .Call("jomo1ranmixhr", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betap,up,bpost,upost,covp,opost, covu,cpost,nbetween, Sp,Sup,Y_numcat, ncol(Y_con),a,rngflag, PACKAGE = "jomo")
    betapost[,,(i-1)]=bpost
    upostall[,,(i-1)]=upost
    omegapost[,,(i-1)]=opost
    covupost[,,(i-1)]=cpost
    bpost<-matrix(0,nrow(betap),ncol(betap))
    opost<-matrix(0,nrow(covp),ncol(covp))
    upost<-matrix(0,nrow(up),ncol(up))
    cpost<-matrix(0,nrow(covu),ncol(covu))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),1:ncol(Y_con)]=Yimp2[,1:ncol(Y_con)]
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y_con)+1):ncol(Y)]=Y_cat
    cat("Imputation number ", i, "registered", "\n")
  }
  betapostmean<-apply(betapost, c(1,2), mean)
  upostmean<-apply(upostall, c(1,2), mean)
  omegapostmean<-apply(omegapost, c(1,2), mean)
  covupostmean<-apply(covupost, c(1,2), mean)
  cat("The posterior mean of the fixed effects estimates is:\n")
  print(betapostmean)
  cat("The posterior mean of the random effects estimates is:\n")
  print(upostmean)
  cat("The posterior mean of the level 1 covariance matrices is:\n")
  print(omegapostmean)
  cat("The posterior mean of the level 2 covariance matrix is:\n")
  print(covupostmean)
  imp<-data.frame(imp)
  if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y_cat), sep = "")
  if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y_con), sep = "")
  if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")  
  colnames(imp)<-c(colnamycon,colnamycat,colnamx,colnamz,"clus","id","Imputation")
  return(imp)
}
