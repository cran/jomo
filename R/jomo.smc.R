jomo.smc <-
  function(formula, data, level=rep(1,ncol(data)), beta.start=NULL, l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, a.start=NULL, a.prior=NULL, nburn=1000, nbetween=1000, nimp=5, meth="common", family="binomial",output=1, out.iter=10, model) {
    if (model=="lm") {
        imp<-jomo.lm(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output, out.iter=out.iter)
    } else if (model=="glm") {
      imp<-jomo.glm(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output, out.iter=out.iter, family=family)
    } else if (model=="polr") {
      imp<-jomo.polr(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output, out.iter=out.iter)
    }else if (model=="coxph") {
      imp<-jomo.coxph(formula=formula, data=data, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output, out.iter=out.iter)
    } else if (model=="lmer") {
      imp<-jomo.lmer(formula=formula, data=data, level=level, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, a.start=a.start, a.prior=a.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, meth=meth, output=output, out.iter=out.iter)
    } else if (model=="glmer") {
      imp<-jomo.glmer(formula=formula, data=data, level=level, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, a.start=a.start, a.prior=a.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, meth=meth, output=output, out.iter=out.iter, family=family)
    } else if (model=="clmm") {
      imp<-jomo.clmm(formula=formula, data=data, level=level, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, a.start=a.start, a.prior=a.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, meth=meth, output=output, out.iter=out.iter)
    }else {
      cat("Invalid model specification. Models currently available: lm, glm (binomial), polr, coxph, lmer,clmm, glmer (binomial).\n")
    }
    return(imp)
  }
