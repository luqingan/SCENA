library(QUIC)

update.lambda = function(lambdas,nedges.seq,nedges){
  lambda = 0
  if(length(lambdas)==length(nedges.seq) & length(nedges.seq)>0){
    lambdas=log(lambdas)
    nedges.seq=log(nedges.seq+1)
    nedges=log(nedges+1)
    if(nedges<max(nedges.seq) & nedges>min(nedges.seq)){
      ind1 = max(which(nedges.seq>nedges))
      ind2 = min(which(nedges.seq<nedges))
      nedges1 = nedges.seq[ind1]
      nedges2 = nedges.seq[ind2]
      lambda1 = max(lambdas[ind1])
      lambda2 = min(lambdas[ind2])
      b = (nedges2-nedges1)/(lambda2-lambda1)
      a = nedges1 - b*lambda1
      lambda = (nedges-a)/b
    }
    if(nedges<min(nedges.seq)){
      lambda = max(lambdas)+1
    }
    if(nedges>max(nedges.seq)){
      lambda = min(lambdas)-1
    }
  }
  return(exp(lambda))
}


Glasso.nedge = 	function(COV,nedges,lambda.start=1,
                         max.trial=100,WEIGHT=1,OBS='full',edge.tol=0,kappa=1,Plot=FALSE,...){
  
  p = ncol(COV)
  if(OBS[1]!='full'){
    print(paste0('computing Glasso with ',nedges,' edges in O'))	
  }
  if(OBS[1]=='full'){
    OBS=diag(p)*0
    OBS=OBS==0
    print(paste0('computing Glasso with ',nedges,' edges'))
  }
  
  lambda = lambda.start
  GLASSO = QUIC(COV,rho=lambda*WEIGHT,msg=0,...)$X
  glassonedge = (sum(GLASSO[OBS]!=0)-p)/2
  trial = 0
  
  lambdas = lambda
  nedges.seq = glassonedge
  
  while(abs(glassonedge-nedges)>edge.tol & trial<max.trial){
    lambda = update.lambda(lambdas=lambdas,nedges.seq=nedges.seq,nedges=nedges)
    GLASSO = QUIC(COV,rho=lambda*WEIGHT,msg=0,...)$X
    glassonedge = (sum(GLASSO[OBS]!=0)-p)/2
    lambdas = c(lambdas,lambda)
    nedges.seq = c(nedges.seq,glassonedge)
    nedges.seq = nedges.seq[order(lambdas)]
    lambdas = sort(lambdas)
    trial=trial+1
    print(c(trial,glassonedge,lambda))
    if(Plot==TRUE){
      plot(log(lambdas),sqrt(nedges.seq),axes=FALSE,xlab=expression(log(lambda)),ylab='# edges')
      axis(1)
      axis(2,labels=sort(unique(nedges.seq)),at=sqrt(sort(unique(nedges.seq))))
      abline(v=log(lambda),h=sqrt(nedges),lty=2)
    }
  }
  
  return(list(GLASSO = GLASSO,lambda=lambda,lambdas=lambdas,glassonedge=nedges.seq))
}


graph.diff = function(THETA.ref,THETA.hat){
  E.ref = round(THETA.ref!=0)
  E.hat = round(THETA.hat!=0)
  d = ncol(E.ref)
  Q=upper.tri(diag(d))
  fdp = sum(E.hat[Q]*(1-E.ref[Q]),na.rm=TRUE)/max(c(1,sum(E.hat[Q],na.rm=TRUE)))
  fnp = sum((1-E.hat[Q])*E.ref[Q],na.rm=TRUE)/max(c(1,sum(1-E.hat[Q],na.rm=TRUE)))
  sens = sum(E.hat[Q]*E.ref[Q],na.rm=TRUE)/max(c(1,sum(E.ref[Q],na.rm=TRUE)))
  spec = sum((1-E.hat[Q])*(1-E.ref[Q]),na.rm=TRUE)/max(c(1,sum(1-E.ref[Q],na.rm=TRUE)))
  
  tot.err = 1-mean(E.hat[Q]*E.ref[Q]+(1-E.hat[Q])*(1-E.ref[Q]))
  
  return(c(fdp=fdp,fnp=fnp,sens=sens,spec=spec,tot.err=tot.err))
}

f1_score = function(THETA.ref,THETA.hat){
  E.ref = round(THETA.ref!=0)
  E.hat = round(THETA.hat!=0)
  d = ncol(E.ref)
  Q=upper.tri(diag(d))
  
  precision  = sum(E.hat[Q]*E.ref[Q],na.rm=TRUE)/max(c(1,sum(E.hat[Q],na.rm=TRUE)))
  recall = sum((1-E.hat[Q])*(1-E.ref[Q]),na.rm=TRUE)/max(c(1,sum(1-E.hat[Q],na.rm=TRUE)))
  f1 = 2*(precision*recall)/(recall+precision)
  return(f1)
}



## FUNCTIONS USING EBIC

EBIC = function(Theta.glasso,Sigma.hat,n,gamma=.5,tol=10^(-10)){
  p = ncol(Theta.glasso)
  LAMBDA = round(abs(cov2cor(Theta.glasso))<=tol)*1000000
  K = sum(LAMBDA[upper.tri(diag(p))]==0)
  Theta.mle = QUIC::QUIC(Sigma.hat,rho=LAMBDA)$X
  val = -n*(determinant(Theta.mle)$modulus-sum(diag(Sigma.hat%*%Theta.mle)))+K*(log(n)+4*gamma*log(p))
  return(val)
}


Glasso.ebic = 	function(COV,n,lambdas,MCCORES=2,gamma=.5,tol=10^(-10)){
  p = ncol(COV)
  THETAS = parallel::mclapply(lambdas,function(lambda) QUIC::QUIC(COV,rho=lambda*(1-diag(p)))$X,mc.cores=MCCORES)
  KS = sapply(1:length(THETAS),function(i) sum(abs(cov2cor(THETAS[[i]])[upper.tri(diag(p))])>tol))
  EBIC = unlist(parallel::mclapply(1:length(THETAS),function(i) EBIC(Theta.glasso=THETAS[[i]],Sigma.hat=COV,n=n,gamma=gamma,tol=tol),mc.cores=MCCORES))
  lambda.opt = lambdas[which.min(EBIC)]
  k.opt = KS[which.min(EBIC)]
  Theta.opt0 = THETAS[[which.min(EBIC)]]
  LAMBDA = round(abs(cov2cor(Theta.opt0))<=tol)*1000000
  Theta.opt = QUIC::QUIC(COV,rho=LAMBDA)$X
  par(mfrow=c(1,2))
  plot(log(lambdas),EBIC,xlab=expression(log(lambda)),type='b',pch=20,main=paste0('opt.log.lambda = ',log(lambda.opt)))
  abline(v=log(lambda.opt),col='blue')
  plot(KS,EBIC,xlab='# edges',type='b',pch=20,main=paste0('opt.n.edges = ',k.opt))
  abline(v=k.opt,col='blue')
  return(list(Theta.opt0=Theta.opt0,Theta.opt=Theta.opt,Graph.opt=round(abs(cov2cor(Theta.opt))>tol),lambda.opt=lambda.opt,k.opt=k.opt))
}


