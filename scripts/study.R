suppressPackageStartupMessages(library(quantreg))
library(parallel)
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork=TRUE)
root <- dirname(dirname(script_path))
workers <- if (.Platform$OS.type == "windows") 1L else 8L
out <- file.path(root,'results','study')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
taus <- c(.1,.5,.9)
ns <- c(128L,256L,512L,1024L)
DES <- lapply(ns,function(n)lapply(seq_len(n/2-1),function(j)cbind(1,cos(2*pi*j*(1:n)/n),sin(2*pi*j*(1:n)/n))))
names(DES) <- ns
qser <- function(y){
 n <- length(y); Z <- matrix(0+0i,n,length(taus)); qh <- as.numeric(quantile(y,taus,type=1)); Z[1,] <- n*qh
 for(j in seq_along(DES[[as.character(n)]]))for(k in seq_along(taus)){
  cf <- rq.fit.fnb(DES[[as.character(n)]][[j]],y,taus[k])$coefficients
  Z[j+1,k] <- n/2*(cf[2]-1i*cf[3]); Z[n-j+1,k] <- Conj(Z[j+1,k])
 }
 for(k in seq_along(taus))Z[n/2+1,k] <- n/2*(quantile(y[seq(2,n,2)],taus[k],type=1)-quantile(y[seq(1,n,2)],taus[k],type=1))
 X <- apply(Z,2,function(z)Re(fft(z,inverse=TRUE))/n); X <- X[c(2:n,1),,drop=FALSE]
 list(x=X,q=qh,z=Z)
}
covseq <- function(x,M){x <- x-mean(x);n<-length(x);sapply(0:M,function(h)sum(x[(h+1):n]*x[1:(n-h)])/n)}
ghat <- function(y,q){bw<-bw.nrd0(y);sapply(q,function(v)mean(dnorm((y-v)/bw))/bw)}
gen_ar <- function(n){x<-numeric(n);x[1]<-rnorm(1);ee<-rnorm(n-1);for(t in 2:n)x[t]<-.6*x[t-1]+.8*ee[t-1];x}
q <- qnorm(taus);g <- dnorm(q);v <- taus*(1-taus)/g^2
truth <- sapply(seq_along(taus),function(k)sapply(1:80,function(h)integrate(function(z)exp(-q[k]^2/(1+z))/(2*pi*sqrt(1-z*z)),0,.6^h,rel.tol=1e-10,abs.tol=1e-14)$value/(taus[k]*(1-taus[k]))))
reps <- 400L
gaussian_job <- function(n,model,rep){
 set.seed(481000+model*100000+n*700+rep)
 y<-if(model==1)rnorm(n) else gen_ar(n); fit<-qser(y);gh<-ghat(y,fit$q);M<-c('128'=5,'256'=6,'512'=7)[as.character(n)]
 rows<-lapply(seq_along(taus),function(k){
  C<-taus[k]-(y<=q[k]);Ch<-taus[k]-(y<=fit$q[k]);U<-q[k]+C/g[k];Uh<-fit$q[k]+Ch/gh[k]
  a<-covseq(fit$x[,k],M);b<-covseq(Ch,M);tar<-if(model==1)rep(0,80) else truth[,k]
  aq<-c((1-(1:M)/(M+1))*a[-1]/a[1],rep(0,80-M));ac<-c((1-(1:M)/(M+1))*b[-1]/b[1],rep(0,80-M))
  data.frame(n=n,model=model,rep=rep,tau=taus[k],oracle=mean((fit$x[,k]-U)^2)/v[k],feasible=mean((fit$x[,k]-Uh)^2)/v[k],threshold=mean((fit$x[,k]-fit$q[k]-Ch/g[k])^2)/v[k],qser_spectral=2*sum((aq-tar)^2),cross_spectral=2*sum((ac-tar)^2))
 });do.call(rbind,rows)
}
allg<-list()
for(n in c(128L,256L,512L))for(model in 1:2){
 z<-mclapply(1:reps,function(rep)gaussian_job(n,model,rep),mc.cores=workers,mc.preschedule=TRUE)
 stopifnot(!any(vapply(z,inherits,logical(1),'try-error')));allg[[length(allg)+1]]<-do.call(rbind,z)
}
write.csv(do.call(rbind,allg),file.path(out,'gaussian_records.csv'),row.names=FALSE)
scale_job <- function(rep){
 set.seed(891000+rep);y<-gen_ar(256);ys<-list(y,exp(y),qt(pnorm(y),df=3));qs<-list(q,exp(q),qt(taus,df=3));gs<-list(g,g/exp(q),dt(qt(taus,df=3),df=3))
 rows<-list()
 for(m in 1:3){f<-qser(ys[[m]]);for(k in seq_along(taus)){
  C<-taus[k]-(y<=q[k]);U<-qs[[m]][k]+C/gs[[m]][k]
  rows[[length(rows)+1]]<-data.frame(rep=rep,margin=m,tau=taus[k],variance=mean((f$x[,k]-f$q[k])^2),loss=mean((f$x[,k]-U)^2)*gs[[m]][k]^2/(taus[k]*(1-taus[k])))
 }};do.call(rbind,rows)
}
z<-mclapply(1:reps,scale_job,mc.cores=workers,mc.preschedule=TRUE);stopifnot(!any(vapply(z,inherits,logical(1),'try-error')))
write.csv(do.call(rbind,z),file.path(out,'transformation_records.csv'),row.names=FALSE)
freq<-seq(0,.5,length.out=2049)
eegrows<-list();spectra<-list()
for(cond in c('open','closed')){
 y<-read.csv(file.path(root,'data',paste0(cond,'_o1.csv')))[,1];y<-y[801:(800+8192)]
 for(n in c(512L,1024L))for(b in seq_len(length(y)/n)){
  x<-y[((b-1)*n+1):(b*n)];x<-resid(lm(x~seq_along(x)));fit<-qser(x);gh<-ghat(x,fit$q);M<-n/8
  for(k in seq_along(taus)){
   Ch<-taus[k]-(x<=fit$q[k]);a<-fit$x[,k]-fit$q[k];c0<-Ch-mean(Ch);d<-gh[k]*a-c0
   ac<-covseq(a,M);cc<-covseq(c0,M);w<-1-(1:M)/(M+1)
   spq<-as.numeric(1+2*(cos(2*pi*outer(freq,1:M))%*%(w*ac[-1]/ac[1])))
   spc<-as.numeric(1+2*(cos(2*pi*outer(freq,1:M))%*%(w*cc[-1]/cc[1])))
   aa<-8/160;bb<-13/160;integ<-(sin(2*pi*bb*(1:M))-sin(2*pi*aa*(1:M)))/(pi*(1:M))
   shareq<-2*((bb-aa)+sum(w*ac[-1]/ac[1]*integ));sharec<-2*((bb-aa)+sum(w*cc[-1]/cc[1]*integ))
   eegrows[[length(eegrows)+1]]<-data.frame(condition=cond,n=n,window=b,tau=taus[k],rrmsd=sqrt(mean(d^2)/(taus[k]*(1-taus[k]))),correlation=cor(a,c0),qser_sd_uv=sqrt(ac[1]),qser_share=shareq,cross_share=sharec)
   if(n==1024L)spectra[[length(spectra)+1]]<-data.frame(condition=cond,n=n,window=b,tau=taus[k],f=freq,qser=spq,crossing=spc)
  }
 }
}
write.csv(do.call(rbind,eegrows),file.path(out,'eeg_windows.csv'),row.names=FALSE)
write.csv(do.call(rbind,spectra),file.path(out,'eeg_spectra.csv'),row.names=FALSE)
