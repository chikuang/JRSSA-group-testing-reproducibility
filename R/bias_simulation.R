# Joint prevalence/sensitivity/specificity Monte Carlo study.
# Reproduces Main Table 10; Supplementary Tables S10-S13 and Figures S5-S6.
# Usage: Rscript R/bias_simulation.R vendor/gtDesign results/bias 10000 4
# Design construction uses gtDesign. The package has no joint MLE routine.
# The vectorized likelihood and gradient below include accuracy = 0.5, which
# gtDesign's model-input validation excludes; they are checked against the
# package on the interior parameter space by scripts/check_package_use.R.
# Custom rounding retains the published fully-spent integer comparison.
# Output-to-paper mapping: PAPER_REPRODUCTION.md and paper_items.csv.
args <- commandArgs(TRUE)
if(length(args)<2) stop('Supply gtDesign source directory and output directory')
src <- normalizePath(args[1]); out <- args[2]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
B <- if(length(args)>2) as.integer(args[3]) else 5000L
cores <- if(length(args)>3) as.integer(args[4]) else 4L
if(.Platform$OS.type=="windows") cores<-1L
stopifnot(B>=2,cores>=1)
source("R/load_gtdesign.R")
load_gtdesign(src)
options(digits=12)
prob <- function(x,theta) theta[2]-(theta[2]+theta[3]-1)*exp(x*log1p(-theta[1]))
jac <- function(x,theta) {
 z <- exp(x*log1p(-theta[1]))
 cbind(x*(theta[2]+theta[3]-1)*z/(1-theta[1]),1-z,-z)
}
fisher <- function(x,n,theta) {
 pi <- prob(x,theta); F <- jac(x,theta)
 crossprod(F,F*as.numeric(n/(pi*(1-pi))))
}
# Same constrained joint MLE for every design; no generating parameter is passed.
# Log prevalence stabilises optimisation. Bounds approximate the closure of the
# manuscript parameter space: prevalence [1e-7,1-1e-7], accuracy [.5,1].
fit_mle <- function(y,n,x,extra_starts=FALSE) {
 lower <- c(log(1e-7),.5,.5); upper <- c(log(1-1e-7),1,1)
 decode <- function(v) c(exp(v[1]),v[2:3])
 fn <- function(v) {
  pi <- pmin(1-1e-14,pmax(1e-14,prob(x,decode(v))))
  -sum(y*log(pi)+(n-y)*log1p(-pi))
 }
 gr <- function(v) {
  th <- decode(v); pi <- pmin(1-1e-14,pmax(1e-14,prob(x,th)))
  J <- jac(x,th); J[,1] <- J[,1]*th[1]
  -as.numeric(crossprod(J,(y-n*pi)/(pi*(1-pi))))
 }
 # When three empirical probabilities yield a feasible saturated solution,
 # it is a global maximum of the binomial likelihood.
 saturated <- NULL
 if(length(x)==3 && all(y>0 & y<n)) {
  v <- y/n; rr <- (v[2]-v[1])/(v[3]-v[1]); rmin <- (x[2]-x[1])/(x[3]-x[1])
  if(all(diff(v)>0) && is.finite(rr) && rr>rmin && rr<1) {
   ratio <- function(p) (-expm1((x[2]-x[1])*log1p(-p)))/(-expm1((x[3]-x[1])*log1p(-p)))-rr
   rt <- tryCatch(uniroot(ratio,c(1e-7,1-1e-7),tol=1e-12)$root,error=function(e)NA_real_)
   if(is.finite(rt)) {
    z <- exp(x*log1p(-rt)); aa <- (v[2]-v[1])/(z[1]-z[2]); ss <- v[1]+aa*z[1]
    th <- c(rt,ss,1+aa-ss)
    if(all(th[2:3]>=.5 & th[2:3]<=1)) saturated <- th
   }
  }
 }
 if(!is.null(saturated) && !extra_starts) {
  vv <- c(log(saturated[1]),saturated[2:3])
  return(c(saturated,nll=fn(vv),convergence=0,pg=max(abs(gr(vv))),saturated=1))
 }
 starts <- c(.003,.02,.10,.40)
 if(extra_starts) starts <- sort(unique(c(starts,.0001,.001,.007,.04,.07,.2,.7,.95)))
 best <- NULL
 for(p in starts) {
  z <- exp(x*log1p(-p))
  # y/n = sensitivity*(1-z) - specificity*z + z.
  seed <- tryCatch(as.numeric(qr.solve(cbind(1-z,-z),y/n-z)),error=function(e)c(.85,.95))
  seed <- pmin(.999,pmax(.501,seed))
  a <- tryCatch(nlminb(c(log(p),seed),fn,gr,lower=lower,upper=upper,
      control=list(iter.max=200,eval.max=400,rel.tol=1e-11,x.tol=1e-9)),error=function(e)NULL)
  if(!is.null(a) && is.finite(a$objective) && (is.null(best)||a$objective<best$objective)) best <- a
 }
 if(is.null(best)) return(c(rep(NA_real_,4),convergence=99,pg=Inf,saturated=0))
 v <- best$par; grad <- gr(v)
 grad[v<=lower+1e-7 & grad>0] <- 0; grad[v>=upper-1e-7 & grad<0] <- 0
 pg <- max(abs(grad))
 # Polish numerical terminations with a second optimiser using the same likelihood.
 if(best$convergence!=0 || pg>1e-4) {
  b <- tryCatch(optim(v,fn,gr,method='L-BFGS-B',lower=lower,upper=upper,
      control=list(maxit=1000,factr=1e4,pgtol=1e-8)),error=function(e)NULL)
  if(!is.null(b) && b$value<=best$objective+1e-9) {
   v <- b$par; best$objective <- b$value; best$convergence <- b$convergence
   grad <- gr(v); grad[v<=lower+1e-7 & grad>0]<-0;grad[v>=upper-1e-7 & grad<0]<-0;pg<-max(abs(grad))
  }
 }
 c(decode(v),nll=best$objective,convergence=best$convergence,pg=pg,saturated=0)
}
# Fully spent integer budget. D allocation maximises determinant among nearby
# integer counts on the optimal approximate support; no exact-global claim.
round_three <- function(x,w,C,q,theta) {
 cost <- 1-q+q*x; target <- C*w/cost
 grid <- expand.grid(b=seq(max(1,floor(target[2])-6),ceiling(target[2])+6),
                     c=seq(max(1,floor(target[3])-6),ceiling(target[3])+6))
 grid$a <- (C-grid$b*cost[2]-grid$c*cost[3])/cost[1]
 grid <- grid[grid$a>=1 & abs(grid$a-round(grid$a))<1e-7,]
 ns <- cbind(round(grid$a),grid$b,grid$c)
 vals <- apply(ns,1,function(n) as.numeric(determinant(fisher(x,n,theta),TRUE)$modulus))
 as.numeric(ns[which.max(vals),])
}
round_uniform <- function(x,C,q) {
 cost <- 1-q+q*x; target <- C/sum(cost)
 ranges <- rep(list(seq(max(1,floor(target)-3),ceiling(target)+3)),3)
 grid <- do.call(expand.grid,ranges)
 aa <- (C-as.numeric(as.matrix(grid)%*%cost[-1]))/cost[1]
 ns <- cbind(aa,as.matrix(grid)); ok <- aa>=1 & abs(aa-round(aa))<1e-7
 ns <- round(ns[ok,,drop=FALSE]); val <- rowSums((ns-target)^2)
 as.numeric(ns[which.min(val),])
}
contexts <- list(cost=list(theta=c(.07,.93,.96),M=61,q=.2,budgets=c(100,500,2000)),
                 low=list(theta=c(.022,.92,.965),M=15,q=0,budgets=c(500,2000,10000)))
cells<-list();designs<-list();certs<-list()
for(ctx in names(contexts)) {
 cc<-contexts[[ctx]]; u<-seq_len(cc$M)
 FF<-t(vapply(u,gtDesign::gt_huang2020_regressor(cc$theta,cc$q),numeric(3)))
 transform<-solve(chol(crossprod(FF))); FFw<-FF%*%transform
 a<-gtDesign::calc_Dopt(u,function(x)as.numeric(FFw[x,]),ridge=0,drop_tol=0,
             tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10)
 w<-a$weights;w[w<1e-5]<-0;w<-w/sum(w); ix<-which(w>0)
 stopifnot(length(ix)==3)
 I<-crossprod(FFw,FFw*w)
 residual<-max(rowSums((FFw%*%solve(I))*FFw)/3)-1
 certs[[ctx]]<-data.frame(context=ctx,status=a$status,max_relative_violation=residual,
                        support=paste(ix,collapse=';'),budget_weights=paste(signif(w[ix],9),collapse=';'))
 stopifnot(residual<1e-4)
 for(C in cc$budgets) for(des in c('D','U4')) {
  x<-if(des=='D') ix else as.integer(round(seq(1,cc$M,length.out=4)))
  n<-if(des=='D') round_three(x,w[ix],C,cc$q,cc$theta) else round_uniform(x,C,cc$q)
  J<-fisher(x,n,cc$theta); E<-eigen(J,symmetric=TRUE)$values
  stopifnot(all(n>0),abs(sum(n*(1-cc$q+cc$q*x))-C)<1e-7,min(E)>0)
  eff<-exp((as.numeric(determinant(J/C,TRUE)$modulus)-as.numeric(determinant(crossprod(FF,FF*w),TRUE)$modulus))/3)
  cell<-list(context=ctx,C=C,design=des,x=x,n=n,theta=cc$theta,q=cc$q,M=cc$M)
  cells[[length(cells)+1]]<-cell
  designs[[length(designs)+1]]<-data.frame(context=ctx,C=C,design=des,x=x,n=n,cost=1-cc$q+cc$q*x,
      spent=sum(n*(1-cc$q+cc$q*x)),assays=sum(n),specimens=sum(n*x),D_eff=eff,
      rank=qr(J)$rank,lambda_min=min(E),condition_number=max(E)/min(E))
 }
}
write.csv(do.call(rbind,designs),file.path(out,'allocations.csv'),row.names=FALSE)
write.csv(do.call(rbind,certs),file.path(out,'design_certificates.csv'),row.names=FALSE)
# Reproducible independent cell streams; estimator starts do not use the truth.
RNGkind('L\'Ecuyer-CMRG');set.seed(20261110)
seeds<-vector('list',length(cells)); seeds[[1]]<-.Random.seed
if(length(cells)>1) for(k in 2:length(cells))seeds[[k]]<-parallel::nextRNGStream(seeds[[k-1]])
results<-list(); audits<-list()
for(k in seq_along(cells)) {
 cc<-cells[[k]]; assign('.Random.seed',seeds[[k]],envir=.GlobalEnv)
 yy<-sapply(seq_along(cc$x),function(j)rbinom(B,cc$n[j],gtDesign::gt_huang2020_pi(cc$x[j],cc$theta)))
 tm<-proc.time()[3]
 ans<-parallel::mclapply(seq_len(B),function(i)fit_mle(yy[i,],cc$n,cc$x),mc.cores=cores,mc.set.seed=FALSE)
 est<-do.call(rbind,ans)
 # Audit the first 50 (preselected) datasets and every problematic termination
 # with an expanded grid of initial values. Retain boundary fits.
 audit_ids<-sort(unique(c(seq_len(min(50,B)),which(est[,'convergence']!=0|est[,'pg']>1e-3|!is.finite(est[,1])))))
 au<-lapply(audit_ids,function(i){
   b<-fit_mle(yy[i,],cc$n,cc$x,TRUE)
   gain<-est[i,'nll']-b['nll']; delta<-max(abs(est[i,1:3]-b[1:3]))
   if(is.finite(gain) && gain>1e-7)est[i,]<<-b
   data.frame(id=i,nll_gain=gain,max_parameter_difference=delta,expanded_pg=b['pg'])
 })
 aa<-do.call(rbind,au);aa$context<-cc$context;aa$C<-cc$C;aa$design<-cc$design;audits[[k]]<-aa
 df<-data.frame(context=cc$context,C=cc$C,design=cc$design,replicate=seq_len(B),
    p=est[,1],sensitivity=est[,2],specificity=est[,3],nll=est[,'nll'],
    convergence=est[,'convergence'],pg=est[,'pg'],saturated=est[,'saturated'])
 results[[k]]<-df
 saveRDS(list(cell=cc,y=yy,estimates=df,audit=aa,seed=seeds[[k]]),file.path(out,sprintf('cell_%02d.rds',k)))
 cat(sprintf('%s C=%d %s: %d fits, %.1fs, boundary=%.3f, max audit gain=%.3g, pg>1e-3=%d\n',
   cc$context,cc$C,cc$design,B,proc.time()[3]-tm,mean(apply(est[,2:3],1,function(z)any(z<.50001|z>.99999))),
   max(aa$nll_gain,na.rm=TRUE),sum(est[,'pg']>1e-3)));flush.console()
}
res<-do.call(rbind,results);write.csv(res,file.path(out,'replicates.csv'),row.names=FALSE)
write.csv(do.call(rbind,audits),file.path(out,'optimizer_audit.csv'),row.names=FALSE)
summary<-list();diagout<-list()
for(k in seq_along(cells)) {
 cc<-cells[[k]];d<-results[[k]];stopifnot(all(is.finite(as.matrix(d[,c('p','sensitivity','specificity')]))))
 for(j in 1:3) {
  value<-d[[c('p','sensitivity','specificity')[j]]];e<-value-cc$theta[j];R<-length(e)
  bias<-mean(e);variance<-var(value);mse<-mean(e^2);rmse<-sqrt(mse)
  # Influence-function MCSE for variance and delta-method MCSE for RMSE.
  v_mcse<-sd((value-mean(value))^2)/sqrt(R)
  summary[[length(summary)+1]]<-data.frame(context=cc$context,C=cc$C,design=cc$design,
    parameter=c('prevalence','sensitivity','specificity')[j],truth=cc$theta[j],replicates=R,
    mean=mean(value),bias=bias,bias_mcse=sd(value)/sqrt(R),variance=variance,variance_mcse=v_mcse,
    rmse=rmse,rmse_mcse=sd(e^2)/sqrt(R)/(2*rmse),bias_fraction_mse=bias^2/mse,
    lower_boundary=mean(value<=c(1e-7,.5,.5)[j]+1e-6),upper_boundary=mean(value>=c(1-1e-7,1,1)[j]-1e-6))
 }
 diagout[[k]]<-data.frame(context=cc$context,C=cc$C,design=cc$design,replicates=B,
   nonzero_optimizer_code=sum(d$convergence!=0),large_projected_gradient=sum(d$pg>1e-3),
   any_accuracy_boundary=mean(d$sensitivity<.500001|d$sensitivity>.999999|d$specificity<.500001|d$specificity>.999999))
}
write.csv(do.call(rbind,summary),file.path(out,'summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagout),file.path(out,'diagnostics.csv'),row.names=FALSE)
saveRDS(list(contexts=contexts,cells=cells,seeds=seeds,replicates=B),file.path(out,'configuration.rds'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
