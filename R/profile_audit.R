# Independent likelihood audit: profile prevalence, solve convex accuracy subproblem.
args<-commandArgs(TRUE);e<-new.env(parent=globalenv())
for(expr in parse('R/bias_simulation.R')){
 eval(expr,e);if(exists('round_uniform',envir=e,inherits=FALSE))break
}
profile_fit<-function(y,n,x,phat){
 fn<-function(lp){
  p<-exp(lp);z<-exp(x*log1p(-p));H<-cbind(1-z,-z)
  obj<-function(a){pr<-pmin(1-1e-14,pmax(1e-14,z+as.numeric(H%*%a)));-sum(y*log(pr)+(n-y)*log1p(-pr))}
  gr<-function(a){pr<-pmin(1-1e-14,pmax(1e-14,z+as.numeric(H%*%a)));-as.numeric(crossprod(H,(y-n*pr)/(pr*(1-pr))))}
  start<-tryCatch(qr.solve(H,y/n-z),error=function(err)c(.85,.95));start<-pmax(.501,pmin(.999,start))
  fit<-nlminb(start,obj,gr,lower=c(.5,.5),upper=c(1,1),control=list(iter.max=200,eval.max=400,rel.tol=1e-12,x.tol=1e-10))
  fit$objective
 }
 pp<-sort(unique(c(exp(seq(log(1e-7),log(1-1e-7),length.out=51)),seq(.005,.2,length.out=30),phat)))
 ll<-log(pp);values<-vapply(ll,fn,numeric(1));best<-min(values)
 mins<-which(values<=c(Inf,head(values,-1)) & values<=c(tail(values,-1),Inf))
 for(i in mins){if(i>1&&i<length(ll))best<-min(best,optimize(fn,c(ll[i-1],ll[i+1]),tol=1e-9)$objective)}
 best
}
files<-list.files(e$out,pattern='^cell_[0-9]+[.]rds$',full.names=TRUE)
ans<-parallel::mclapply(files,function(file){
 a<-readRDS(file);cc<-a$cell;d<-a$estimates
 ids<-sort(unique(c(seq_len(min(50,nrow(d))),which(d$pg>1e-3))))
 do.call(rbind,lapply(ids,function(i){
   value<-profile_fit(a$y[i,],cc$n,cc$x,d$p[i])
   data.frame(context=cc$context,C=cc$C,design=cc$design,id=i,
     joint_nll=d$nll[i],profile_nll=value,profile_gain=d$nll[i]-value,
     normalized_pg=d$pg[i]/sum(cc$n))
 }))
},mc.cores=e$cores,mc.set.seed=FALSE)
d<-do.call(rbind,ans);write.csv(d,file.path(e$out,'profile_audit.csv'),row.names=FALSE)
print(aggregate(profile_gain~context+C+design,d,max));
cat('Maximum independent-profile likelihood gain:',max(d$profile_gain),'\n')
stopifnot(max(d$profile_gain)<1e-6)
