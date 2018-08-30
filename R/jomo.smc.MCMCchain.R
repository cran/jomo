jomo.smc.MCMCchain <-
  function(formula, data, level=rep(1,ncol(data)), beta.start=NULL, l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, a.start=NULL, a.prior=NULL, betaY.start=NULL, varY.start=NULL, covuY.start=NULL, uY.start=NULL,  nburn=1000, meth="common", family="binomial", start.imp=NULL, start.imp.sub=NULL, l2.start.imp=NULL, output=1, out.iter=10, model) {
    if (model=="lm") {
      imp<-jomo.lm.MCMCchain(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, betaY.start=betaY.start, varY.start=varY.start, nburn=nburn,  start.imp=start.imp, start.imp.sub=start.imp.sub, output=output, out.iter=out.iter)
    } else if (model=="glm") {
      imp<-jomo.glm.MCMCchain(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, betaY.start=betaY.start, nburn=nburn, start.imp=start.imp, start.imp.sub=start.imp.sub, output=output, out.iter=out.iter, family=family)
    } else if (model=="polr") {
      imp<-jomo.polr.MCMCchain(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, betaY.start=betaY.start, nburn=nburn, start.imp=start.imp, start.imp.sub=start.imp.sub, output=output, out.iter=out.iter)
    }else if (model=="coxph") {
      imp<-jomo.coxph.MCMCchain(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, betaY.start=betaY.start,nburn=nburn,  start.imp=start.imp, output=output, out.iter=out.iter)
    } else if (model=="lmer") {
      imp<-jomo.lmer.MCMCchain(formula=formula, data=data, level=level, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, a.start=a.start, a.prior=a.prior, betaY.start=betaY.start, varY.start=varY.start, covuY.start=covuY.start, uY.start=uY.start, nburn=nburn, meth=meth, start.imp=start.imp, start.imp.sub=start.imp.sub, l2.start.imp=l2.start.imp, output=output, out.iter=out.iter)
    } else if (model=="glmer") {
      imp<-jomo.glmer.MCMCchain(formula=formula, data=data, level=level, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, a.start=a.start, a.prior=a.prior, betaY.start=betaY.start, covuY.start=covuY.start, uY.start=uY.start, nburn=nburn, meth=meth, start.imp=start.imp, start.imp.sub=start.imp.sub, l2.start.imp=l2.start.imp, output=output, out.iter=out.iter, family=family)
    } else if (model=="clmm") {
      imp<-jomo.clmm.MCMCchain(formula=formula, data=data, level=level, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, a.start=a.start, a.prior=a.prior, betaY.start=betaY.start, covuY.start=covuY.start, uY.start=uY.start, nburn=nburn, meth=meth, start.imp=start.imp, start.imp.sub=start.imp.sub, l2.start.imp=l2.start.imp, output=output, out.iter=out.iter)
    }else {
      cat("Invalid model specification. Models currently available: lm, glm (binomial), polr, coxph, lmer, clmm, glmer (binomial).\n")
    }
    return(imp)
  }
