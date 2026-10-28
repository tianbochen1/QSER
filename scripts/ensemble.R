suppressPackageStartupMessages(library(quantreg))
library(parallel)
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork=TRUE)
root <- dirname(dirname(script_path))
wd <- file.path(root, "results", "ensemble")
workers <- if (.Platform$OS.type == "windows") 1L else 8L
dir.create(wd,recursive=TRUE,showWarnings=FALSE)
taus <- c(.1,.5,.9); ns <- c(128L,512L); reps <- 400L
sig <- c(1,3)/sqrt(5)
qpop <- sapply(taus,function(tau)uniroot(function(q)mean(pnorm(q/sig))-tau,c(-10,10),tol=1e-12)$root)
gpop <- sapply(qpop,function(q)mean(dnorm(q/sig)/sig))
DES <- lapply(ns,function(n)lapply(1:(n/2-1),function(j)cbind(1,cos(2*pi*j*(1:n)/n),sin(2*pi*j*(1:n)/n))))
names(DES)<-ns
qser <- function(y){
 n<-length(y);Z<-matrix(0+0i,n,length(taus));qh<-as.numeric(quantile(y,taus,type=1));Z[1,]<-n*qh
 for(j in seq_along(DES[[as.character(n)]]))for(k in seq_along(taus)){
  cf<-rq.fit.fnb(DES[[as.character(n)]][[j]],y,taus[k])$coefficients
  Z[j+1,k]<-n/2*(cf[2]-1i*cf[3]);Z[n-j+1,k]<-Conj(Z[j+1,k])
 }
 for(k in seq_along(taus))Z[n/2+1,k]<-n/2*(quantile(y[seq(2,n,2)],taus[k],type=1)-quantile(y[seq(1,n,2)],taus[k],type=1))
 X<-apply(Z,2,function(z)Re(fft(z,inverse=TRUE))/n);X<-X[c(2:n,1),,drop=FALSE]
 PQ<-apply(X,2,function(x)Mod(fft(x-mean(x)))^2/n)
 err<-max(abs(PQ[-1,]-Mod(Z[-1,])^2/n))
 list(X=X,Z=Z,P=PQ,q=qh,error=err)
}
job <- function(model,n,rep){
 set.seed(812000+model*100000+n*500+rep)
 r<-if(model==1)0 else .8
 V<-integer(n);V[1]<-sample(0:1,1)
 for(t in 2:n)V[t]<-if(runif(1)<(1+r)/2)V[t-1] else 1-V[t-1]
 y<-sig[V+1]*rnorm(n);fit<-qser(y)
 C<-sapply(seq_along(taus),function(k)taus[k]-(y<=fit$q[k]))
 Co<-sapply(seq_along(taus),function(k)taus[k]-(y<=qpop[k]))
 PC<-apply(C,2,function(x)Mod(fft(x-mean(x)))^2/n)
 PO<-apply(Co,2,function(x)Mod(fft(x-mean(x)))^2/n)
 PY<-Mod(fft(y-mean(y)))^2/n
 check<-c(fft=fit$error,cross_scale=0,qp_scale=NA)
 if(rep<=8){
  f2<-qser(2*y);C2<-sapply(seq_along(taus),function(k)taus[k]-(2*y<=f2$q[k]))
  check['cross_scale']<-max(abs(C2-C));check['qp_scale']<-max(abs(f2$P[-1,]-4*fit$P[-1,]))/max(4*fit$P[-1,])
 }
 ii<-2:(n/2)
 list(QP=fit$P[ii,,drop=FALSE],Crossing=PC[ii,,drop=FALSE],Oracle=PO[ii,,drop=FALSE],Ordinary=PY[ii],check=check)
}
for(n in ns)for(model in 1:2){
 ans<-mclapply(1:reps,function(rep)job(model,n,rep),mc.cores=workers,mc.preschedule=TRUE)
 if(any(vapply(ans,inherits,logical(1),'try-error')))stop('A numerical job failed')
 saveRDS(list(n=n,model=model,r=if(model==1)0 else .8,taus=taus,q=qpop,g=gpop,sigma=sig,ans=ans),file.path(wd,sprintf('ensemble_%d_%d.rds',n,model)))
 f<-(1:(n/2-1))/n
 for(method in c('QP','Crossing','Oracle')){
  A<-array(unlist(lapply(ans,'[[',method)),dim=c(length(f),length(taus),reps))
  z<-expand.grid(f=f,tau=taus,rep=1:reps);z$value<-as.vector(A)
  write.csv(z,file.path(wd,sprintf('%s_%d_%d.csv',method,n,model)),row.names=FALSE)
 }
 A<-sapply(ans,'[[','Ordinary');z<-expand.grid(f=f,rep=1:reps);z$value<-as.vector(A)
 write.csv(z,file.path(wd,sprintf('Ordinary_%d_%d.csv',n,model)),row.names=FALSE)
 checks<-do.call(rbind,lapply(ans,'[[','check'));write.csv(checks,file.path(wd,sprintf('checks_%d_%d.csv',n,model)),row.names=FALSE)
}
write.csv(data.frame(tau=taus,q=qpop,g=gpop),file.path(wd,'marginal.csv'),row.names=FALSE)
