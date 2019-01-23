# Inverse function
myginv <- function (X, tol = sqrt(.Machine$double.eps)){
    ## the ginv function for symmetric matrices. it avoid some errors.
        if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
                  stop("'X' must be a numeric or complex matrix")
            if (!is.matrix(X))
                      X <- as.matrix(X)
	    if(!isSymmetric(X))
		      X <- (X + t(X))/2
            Xeigen <- eigen(X,symmetric=T)
            Positive <- Xeigen$values > max(tol * Xeigen$values[1L], 0)
            if (all(Positive))
                      Xeigen$vectors %*% (1/Xeigen$values * t(Xeigen$vectors))
            else if (!any(Positive))
                      array(0, dim(X)[2L:1L])
            else Xeigen$vectors[, Positive, drop = FALSE] %*% ((1/Xeigen$values[Positive]) *
                                              t(Xeigen$vectors[, Positive, drop = FALSE]))
}
  
mtxInvPrx <- function(d,w,LD.matrix,thrw,thrctd=1,thrp=length(d)/2,thrp2=thrp/2,thrp3=5){ ### This function approximates the inverse of D+diag(w)%*%LD.matrix%*%diag(w)
 p <- length(w)
 wsq <- w^2
 rkwsq <- rank(wsq)
 wsqcumsum <- sapply(1:p, function(ii) sum(wsq[rkwsq<=ii]))
 wstat <- wsqcumsum+sqrt((p-(1:p))*wsqcumsum)
 rkthr <- min(max(which(wstat<thrw),p-thrp),p-thrp3) ## at least throw out half of the variables
 ix1 <- (rkwsq>rkthr)
 if(sum(ix1) != 1){
   ctdmax <- apply(abs(LD.matrix-diag(diag(LD.matrix)))[ix1,],2,max)
   ctdmax[ix1] <- 0
   ix1 <- ix1|(ctdmax>=thrctd)
   ix2 <- !ix1 
 
   M0 <- sapply(1:p, function(ii) w[ii]*LD.matrix[ii,])
   M0 <- sapply(1:p, function(ii) w[ii]*M0[ii,])
   M <- diag(d)+M0
   Minv <- diag(1/d)
   if(isSymmetric(M[ix1,ix1])){
     Minv[ix1,ix1] <- myginv(M[ix1,ix1]) 
   }else Minv[ix1,ix1] <- myginv((M[ix1,ix1]+t(M[ix1,ix1]))/2)#myginv(M[ix1,ix1])
   out <- Minv
 }else out <- 'ginv'
 out 
}


vestepDirectSS <- function(posterior, para, LD.matrix, wzy, ycy, anno = matrix(1, length(wzy), 1), thrw.mtxinv = 0.1, thrctd.mtxinv = 0.8, thrp.mtxinv = length(wzy)){
	p <- length(wzy)
	gammabeta <- para$gammabeta
	lodanno <- anno%*%gammabeta  
	qlodbeta <- matrix(posterior$lodbeta,p,1)
	qbeta <- 1/(1+exp(-qlodbeta))

	### update the posterior distribution of the prior distribution of pibeta
	updatePibeta <- function(expanno, qbeta){
		shape1pibeta <- expanno + qbeta
		shape2pibeta <- 1/expanno+1-qbeta
		epibeta <- shape1pibeta/(shape1pibeta+shape2pibeta)
		elogpibeta <- digamma(shape1pibeta) - digamma(shape1pibeta+shape2pibeta)
		elogpibetac <- digamma(shape2pibeta) - digamma(shape1pibeta+shape2pibeta)
		elodpibeta <- elogpibeta - elogpibetac
		out <- list(par = list(shape1pibeta = shape1pibeta, shape2pibeta = shape2pibeta), stat = list(epibeta = epibeta,
			elodpibeta = elodpibeta, elogpibeta = elogpibeta, elogpibetac = elogpibetac))
		out
	}
  
  
	### update the posterior distribution of tau
	updateTau <- function(para, qbeta){
		vareps <- para$vareps
		nutau <- para$nutau
		At <- try(myginv(diag(as.vector(qbeta))%*%LD.matrix%*%diag(as.vector(qbeta))+diag(as.vector(qbeta-qbeta^2))+diag(1,p)/nutau))
		if (class(At) == 'try-error') {
			cat('Caught an error during full inversion, trying approximation.\n')
			At <- mtxInvPrx(d=1/nutau+as.vector(qbeta-qbeta^2), w=qbeta, LD.matrix=LD.matrix, thrw=thrw.mtxinv, thrctd=thrctd.mtxinv, thrp=thrp.mtxinv, thrp2=thrp.mtxinv/2) 
		}
		vartau <- vareps*At
		mutau <- At%*%(qbeta*wzy)
		out <- list(A = At, vartau = vartau, mutau = mutau)
		out
	}
  
	post.tau <- updateTau(para,qbeta)

	### update the posterior distribution of sbeta (signal of each SNP)
	updateSbeta <- function(para, lodprior, posttau, qbeta){
		vareps <- para$vareps
		vartau <- posttau$vartau
		mutau <- posttau$mutau
		muwzy <- LD.matrix%*%(mutau*qbeta)
		res <- (wzy - muwzy)
		post.lodbeta <- lodprior-(1/(2*vareps))*(diag(LD.matrix)*(mutau^2+diag(vartau))-2*mutau*res-2*mutau^2*qbeta+2*(LD.matrix*vartau)%*%qbeta-2*diag(LD.matrix)*diag(vartau)*qbeta)
		out <- list(lodbeta = post.lodbeta)
		out
	}
	
	post.sbeta <- updateSbeta(para, lodanno, posttau = post.tau, qbeta)
  
	out <- list(posterior = c(post.sbeta, post.tau), stat = NULL)
	out
} ### End of "vestepDirectSS"


vmstepDirectSS <- function(post, LD.matrix, wzy, ycy, init = list(gammabeta = rep(0, dim(anno)[2]), vareps = 1, nutau = 5), anno = matrix(1, length(wzy), 1), iter.max.mstep = 100){
	p <- length(wzy)
	mutau <- post[['posterior']]$mutau
	vartau <- post[['posterior']]$vartau
	lodbeta <- post[['posterior']]$lodbeta
	qbeta <- 1/(1+exp(-lodbeta))
	mubeta <- mutau*qbeta
	qdly <- (ycy-2*sum(wzy*mubeta)+sum(mubeta*(LD.matrix%*%mubeta))+sum(qbeta*((LD.matrix*vartau)%*%qbeta))+sum(diag(LD.matrix)*(qbeta-qbeta^2)*(mutau^2+diag(vartau))))
	qdltau <- sum((mutau^2+diag(vartau))*qbeta) 
 
	vareps <- qdly/(p-sum(qbeta))
	nutau <- qdltau/sum(qbeta)/vareps 
  
	loglkpibetaprior <- function(gammabeta){
		lodanno <- anno%*%matrix(gammabeta, ncol = 1)
		out <- -qbeta*lodanno+log(1+exp(lodanno))
		sum(out)
	}
	
	out.loglkpibetaprior <- optim(par = init$gammabeta, fn = loglkpibetaprior, method = "BFGS", hessian = TRUE, control = list(maxit = iter.max.mstep))
	out <- list(par = list(gammabeta = out.loglkpibetaprior$par, vareps = vareps, nutau = nutau), hess = list(pibetaprior = out.loglkpibetaprior$hessian))
	out
} ### End of "vmstepDirectSS"



vemDirectSS <- function(LD.matrix, GEMsummstats, anno = NULL, iter.max = 200, iter.max.mstep = 100, er.max = 1e-5, 
          init = list(gammabeta = c(-1, rep(0, dim(anno)[2] - 1)), vareps = sd(GEMsummstats)^2, nutau = 2), 
          posterior.init = list(lodbeta = matrix(0, length(GEMsummstats), 1))){
  thrw.mtxinv = 0.1 
  thrctd.mtxinv=0.8
  wzy <- GEMsummstats
	p <- length(wzy)
	if(is.null(anno)) anno <- matrix(1, p, 1)
	LD.matrix.inv <- myginv(LD.matrix)
	ycy <- as.numeric(t(wzy)%*%LD.matrix.inv%*%wzy)
	posterior.old <- posterior.init
	wzy <- matrix(wzy,p,1)

	vaicDirect <- function(para,post){
 		mutau <- post[['posterior']]$mutau
		vartau <- post[['posterior']]$vartau
		lodbeta <- post[['posterior']]$lodbeta
		vareps <- para[['vareps']]
		qbeta <- 1/(1+exp(-lodbeta))
		mubeta <- mutau*qbeta
		p <- length(lodbeta)
		qdly <- (ycy-2*sum(wzy*mubeta)+sum(mubeta*(LD.matrix%*%mubeta))+sum(qbeta*((LD.matrix*vartau)%*%qbeta))+sum(diag(LD.matrix)*(qbeta-qbeta^2)*(mutau^2+diag(vartau))))
		ely <- -qdly/(2*vareps) -(p/2)*log(vareps)
		mly <- -(ycy-2*sum(wzy*mubeta)+sum(mubeta*(LD.matrix%*%mubeta)))/(2*vareps) -(p/2)*log(vareps)
		pd <- 2*(mly-ely)
		vaic <- -2*mly+2*pd
		vaic
	}

	er <- Inf
	ll.all <- Inf
	ll.old <- Inf
	iter <- 1
	mstep.all <- list()

	### The algorithm stop until the hyperparameter estimation converges in relative error.
	while((abs(er)>er.max)&(iter<iter.max)){
		estep <- vestepDirectSS(posterior.old, init, LD.matrix, wzy, ycy, anno, thrw.mtxinv = thrw.mtxinv, thrctd.mtxinv = thrctd.mtxinv, thrp.mtxinv = floor(p*(1/4+(3/4)^iter)))
		mstep <- vmstepDirectSS(post = estep, LD.matrix, wzy, ycy, init = init, anno, iter.max.mstep = iter.max.mstep)
		mstep.all[[iter]] <- mstep
		posterior.old <- estep[['posterior']]
		init <- mstep$par
		ll <- vaicDirect(init, estep)
		ll.all <- c(ll.all, ll)
		er <- ll.old - ll
		er <- max(er,er/abs(ll))
		ll.old <- ll
		iter <- iter + 1
	}
	
	
	out <- list(posterior = estep[['posterior']], par = mstep$par, niter = iter, 
		converged = (er<er.max), par.all = mstep.all)
	out
} ### End of "vemDirectSS"



vestepMedSS <- function(posterior,para,LD.matrix,wzy,ycy,wzr,rcr,rcy,anno=matrix(1,length(wzy),1),thrw.mtxinv=0.1,thrctd.mtxinv=0.8,thrp.mtxinv=length(wzy)){
  p <- length(wzy)
  gammabeta <- para$gammabeta
  lodanno <- anno%*%gammabeta  
  qlodbeta <- matrix(posterior$lodbeta,p,1)
  qbeta <- 1/(1+exp(-qlodbeta))
  m1gamma <- posterior$mudelta
  
  ### update the posterior distribution of the prior distribution of pibeta
  updatePibeta <- function(expanno, qbeta){
    shape1pibeta <- expanno+qbeta
    shape2pibeta <- 1/expanno+1-qbeta
    epibeta <- shape1pibeta/(shape1pibeta+shape2pibeta)
    elogpibeta <- digamma(shape1pibeta)-digamma(shape1pibeta+shape2pibeta)
    elogpibetac <- digamma(shape2pibeta)-digamma(shape1pibeta+shape2pibeta)
    elodpibeta <- elogpibeta-elogpibetac
    out <- list(par=list(shape1pibeta=shape1pibeta,shape2pibeta=shape2pibeta),stat=list(epibeta=epibeta,elodpibeta=elodpibeta,elogpibeta=elogpibeta,elogpibetac=elogpibetac))
    out
  }
  
  ### update the posterior distribution of tau
  updateTau <- function(para,qbeta){
    vareps <- para$vareps
    nutau <- para$nutau
    At <- try(myginv(diag(as.vector(qbeta))%*%LD.matrix%*%diag(as.vector(qbeta))+diag(as.vector(qbeta-qbeta^2))+diag(1,p)/nutau))
    if (class(At) == 'try-error') {
      cat('Caught an error during full inversion, trying approximation.\n')
      At <- mtxInvPrx(d=1/nutau+as.vector(qbeta-qbeta^2), w=qbeta, LD.matrix=LD.matrix, thrw=thrw.mtxinv, thrctd=thrctd.mtxinv, thrp=thrp.mtxinv, thrp2=thrp.mtxinv/2) 
    }
    vartau <- vareps*At
    mutau <- At%*%(qbeta*(wzy-m1gamma*wzr))
    out <- list(A=At,vartau=vartau,mutau=mutau)
    out
  }
  
  post.tau <- updateTau(para,qbeta)

  updateDelta <- function(para,posttau,qbeta){
    vareps <- para$vareps
    mutau <- posttau$mutau
    mudelta <- (rcy-sum(wzr*mutau*qbeta))/rcr
    vardelta <- vareps/rcr
    out <- list(mudelta=mudelta,vardelta=vardelta)
    out
  }
  
 post.delta <- updateDelta(para,post.tau,qbeta)
  
  ### update the posterior distribution of sbeta 
  updateSbeta <- function(para,lodprior,posttau,postdelta,qbeta){
    vareps <- para$vareps
    vartau <- posttau$vartau
    mutau <- posttau$mutau
    mudelta <- postdelta$mudelta
    muwzy <- LD.matrix%*%(mutau*qbeta)+mudelta*wzr
    res <- (wzy -muwzy)
    post.lodbeta <- lodprior-(1/(2*vareps))*(diag(LD.matrix)*(mutau^2+diag(vartau))-2*mutau*res-2*mutau^2*qbeta+2*(LD.matrix*vartau)%*%qbeta-2*diag(LD.matrix)*diag(vartau)*qbeta)
    out <- list(lodbeta=post.lodbeta)
    out
  }

  post.sbeta <- updateSbeta(para,lodanno,posttau=post.tau,postdelta=post.delta,qbeta)
 
  out <- list(posterior=c(post.sbeta,post.tau,post.delta),stat=NULL)
  out
} ### End of "vestepMedSS"


vmstepMedSS <- function(post, LD.matrix, wzy, ycy, wzr, rcr, rcy, init = list(gammabeta = rep(0, dim(anno)[2]), vareps = 1, nutau = 5), 
				anno = matrix(1, length(wzy), 1), iter.max.mstep = 100){

	p <- length(wzy)
	mutau <- post[['posterior']]$mutau
	vartau <- post[['posterior']]$vartau
	mudelta <- post[['posterior']]$mudelta
	vardelta <- post[['posterior']]$vardelta
	lodbeta <- post[['posterior']]$lodbeta
	qbeta <- 1/(1+exp(-lodbeta))
	mubeta <- mutau*qbeta
	qdly <- (ycy-2*sum(wzy*mubeta)+sum(mubeta*(LD.matrix%*%mubeta))+sum(qbeta*((LD.matrix*vartau)%*%qbeta))+sum(diag(LD.matrix)*(qbeta-qbeta^2)*(mutau^2+diag(vartau)))+(mudelta^2+vardelta)*rcr-2*mudelta*(rcy-sum(wzr*mubeta)))

	qdltau <- sum((mutau^2+diag(vartau))*qbeta) 
	vareps <- qdly/(p-sum(qbeta)-1)
	nutau <- qdltau/sum(qbeta)/vareps 

	loglkpibetaprior <- function(gammabeta){
		lodanno <- anno%*%matrix(gammabeta, ncol = 1)
		out <- log(1 + exp(lodanno)) - qbeta*lodanno 
		sum(out)
	}
	out.loglkpibetaprior <- optim(par = init$gammabeta, fn = loglkpibetaprior, method = "BFGS", hessian = TRUE, control = list(maxit = iter.max.mstep))
	out <- list(par = list(gammabeta = out.loglkpibetaprior$par, vareps = vareps, nutau = nutau), hess = list(pibetaprior = out.loglkpibetaprior$hessian))
	out
} ### End of "vmstepMedSS"


vemMedSS <- function(LD.matrix, DEMsummstats, GEMsummstats, anno = NULL, iter.max = 200, iter.max.mstep = 100, er.max = 1e-5,
	init = list(gammabeta = c(-1, rep(0, dim(anno)[2] - 1)), vareps = sd(DEMsummstats)^2, nutau = 2), 
	posterior.init = list(lodbeta = matrix(0, length(DEMsummstats), 1), mudelta = 0)){
  thrw.mtxinv = 0.1
  thrctd.mtxinv = 0.8
  wzy <- DEMsummstats
  wzr <- GEMsummstats
	p <- length(wzy)
	if(is.null(anno)) anno <- matrix(1, p, 1)
	LD.matrix.inv <- myginv(LD.matrix)
  ycy <- as.numeric(t(wzy)%*%LD.matrix.inv%*%wzy)
  rcr <- as.numeric(t(wzr)%*%LD.matrix.inv%*%wzr)
  rcy <- as.numeric(t(wzr)%*%LD.matrix.inv%*%wzy)
	
	posterior.old <- posterior.init
	wzy <- matrix(wzy,p,1)

	vaicMed <- function(para,post){
		mutau <- post[['posterior']]$mutau
		vartau <- post[['posterior']]$vartau
		mudelta <- post[['posterior']]$mudelta
		vardelta <- post[['posterior']]$vardelta
		lodbeta <- post[['posterior']]$lodbeta
		vareps <- para[['vareps']]
		qbeta <- 1/(1+exp(-lodbeta))
		mubeta <- mutau*qbeta
		p <- length(lodbeta)

		qdly <- (ycy-2*sum(wzy*mubeta)+sum(mubeta*(LD.matrix%*%mubeta))+sum(qbeta*((LD.matrix*vartau)%*%qbeta))+sum(diag(LD.matrix)*(qbeta-qbeta^2)*(mutau^2+diag(vartau)))+(mudelta^2+vardelta)*rcr-2*mudelta*(rcy-sum(wzr*mubeta)))
 
		ely <- -qdly/(2*vareps) -(p/2)*log(vareps)
		mly <- -(ycy-2*sum(wzy*mubeta)+sum(mubeta*(LD.matrix%*%mubeta))+mudelta^2*rcr-2*mudelta*(rcy-sum(wzr*mubeta)))/(2*vareps) -(p/2)*log(vareps)
		pd <- 2*(mly-ely)
		vaic <- -2*mly+2*pd
		vaic
	}

	er <- Inf
	ll.all <- Inf
	ll.old <- Inf
	iter <- 1
	mstep.all <- list()
	### The algorithm stop until the hyperparameter estimation converges in relative error.
	while((abs(er)>er.max)&(iter<iter.max)){
		estep <- vestepMedSS(posterior.old, init, LD.matrix, wzy, ycy, wzr, rcr, rcy, anno, thrw.mtxinv = thrw.mtxinv, thrctd.mtxinv = thrctd.mtxinv,
			thrp.mtxinv = floor(p*(1/4+(3/4)^iter)))
		mstep <- vmstepMedSS(post = estep,LD.matrix, wzy, ycy, wzr, rcr, rcy, init = init, anno, iter.max.mstep = iter.max.mstep)
		mstep.all[[iter]] <- mstep
		posterior.old <- estep[['posterior']]
		init <- mstep$par
		ll <- vaicMed(init, estep)
		ll.all <- c(ll.all, ll)
		er <- ll.old - ll
		er <- max(er,er/abs(ll))
		ll.old <- ll
		iter <- iter + 1
  }
	out <- list(posterior = estep[['posterior']], par = mstep$par, niter = iter, converged = (er<er.max),
		par.all = mstep.all)
	out
} ### End of "vemMedSS"



getParameters <- function(omod){  ### process the output of the model without mediator
	omod[['par']]
}

getPostProb <- function(omod){  
	PP <- 1/(1+exp(- omod[['posterior']]$lodbeta))
	SNP.name <- rownames(PP)
	PP <- c(PP)
	names(PP) <- SNP.name
	idfdr <- rank(-PP, ties.method = 'max')
	fdr <- sapply(idfdr, function(ii) 1 - mean(PP[idfdr <= ii]))

  out <- PP
	out	
}


processFunMed <- function(GEMoutput, DEMoutput){
  
  mod <- list()
  mod[['model.opt']] <- list(GEM = GEMoutput, DEM = DEMoutput)
 
	mod.out <- list()
	mod.out[['Parameters']] <- lapply(mod[['model.opt']], getParameters)
	names(mod.out[['Parameters']][['GEM']]) <- c('gammaB', 'varEta', 'nuB')
	names(mod.out[['Parameters']][['DEM']]) <- c('gammaBeta', 'varEpsilon', 'nuBeta')
	mod.out[['Parameters']][['DEM']][['delta']] <- mod[['model.opt']][['DEM']][['posterior']]$mudelta
	mod.out[['PostProb']] <- lapply(mod[['model.opt']], getPostProb)
	mod.out[['Convergence']] <- lapply(mod[['model.opt']], function(xx) xx[c('niter','converged')])
  
	mod.out
}



meanAnno <- function(x,y){
  y <- (as.numeric(y)>0)
  if(sum(y)>0){ x1 <- mean(x[y])
  }else{ x1 <- 0}
  x1
}

lbAnno <- function(x, y, post, nrep=1000){  
  y <- as.numeric(y)
  n <- length(y)
  n1 <- sum(y>0)
  shm <- sapply(1:nrep, function(ii){
    idii <- sample(n,n1,replace=FALSE)
    out <- mean(x[idii])
    out })
  out <-  sum(shm >= post)/nrep
  out
}

annoEnrich <- function(iFunMedNull, annoMtx, Nperm = 5000, cores = 20){
  require(parallel)
  idanno <- colnames(annoMtx)
  ## posterior porbabilities from the null model
  post.null <- list(GEM = iFunMedNull[['PostProb']][['GEM']],
                    DEM = iFunMedNull[['PostProb']][['DEM']])
  
  ## average posterior probability of inclusion of the SNPs with the annotation, for each annotation
  post.meananno.null <- lapply(post.null, function(x) apply(annoMtx, 2, meanAnno, x = x))
  
  anno.enrichment <- mclapply(names(post.null), function(x) sapply(1:length(idanno), function(y) lbAnno(post.null[[x]], annotation.matrix[, y], post.meananno.null[[x]][y], nrep = Nperm)), mc.cores = cores, mc.preschedule = FALSE)
  names(anno.enrichment) <- names(post.null)
  names(anno.enrichment[['GEM']]) <- names(post.meananno.null[['GEM']])
  names(anno.enrichment[['DEM']]) <- names(post.meananno.null[['DEM']])
  
  return(list(enrichmentPval = anno.enrichment, avePP = post.meananno.null))
}
