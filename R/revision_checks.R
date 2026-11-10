# Revision checks: source the pinned chikuang/gtDesign R files supplied as argv[1].
# No installed package or GitHub repository is modified.
args <- commandArgs(TRUE)
if (!length(args) %in% c(2L,3L)) stop('Usage: Rscript revision_checks.R /path/to/pinned/gtDesign /path/to/results [/path/to/figures]')
src <- normalizePath(args[1]); out <- normalizePath(args[2])
fig_out <- if (length(args) == 3L) normalizePath(args[3]) else out
for (file in list.files(file.path(src, 'R'), full.names=TRUE, pattern='[.]R$')) source(file)
options(digits=12)
theta <- c(.07,.93,.96); u <- 1:61; q <- .2
f <- gt_huang2020_regressor(theta,q)
F <- t(vapply(u,f,numeric(3)))
cost <- 1-q+q*u
fim <- function(n,C) crossprod(F, F*as.numeric(n*cost/C))
loss <- function(I,k) switch(k,D=exp(-as.numeric(determinant(I,TRUE)$modulus)/3),
 A=sum(diag(solve(I))),Ds=solve(I)[1,1],E=1/min(eigen(I,symmetric=TRUE,only.values=TRUE)$values))
fmt <- function(n) paste(paste0(u[n>0],':',n[n>0]),collapse=';')
res <- list(); cert <- list(); weights <- list()
for(k in c('D','A','Ds','E')) {
  cat('Solving',k,'\n')
  a <- switch(k,D=calc_Dopt(u,f,ridge=0,drop_tol=0),A=calc_Aopt(u,f,drop_tol=0),
              Ds=calc_copt(u,f,cVec=c(1,0,0),drop_tol=0),E=calc_Eopt(u,f,drop_tol=0))
  wfull <- a$weights
  w <- wfull; w[w<1e-6] <- 0; w <- w/sum(w)
  I <- crossprod(F,F*w); V <- solve(I)
  G <- switch(k,D=V/3,A=V%*%V,Ds=tcrossprod(V[,1]),E={
    eg<-eigen(I,symmetric=TRUE); tcrossprod(eg$vectors[,3])/eg$values[3]^2})
  score <- rowSums((F%*%G)*F); avg <- sum(diag(G%*%I))
  cert[[k]] <- data.frame(criterion=k,solver_status=a$status,
    discarded_mass=sum(wfull[wfull<1e-6]),max_relative_violation=max(score/avg)-1)
  weights[[k]] <- data.frame(criterion=k,pool=u,weight=w,full_weight=wfull)
  a$design <- data.frame(point=u[w>0],weight=w[w>0])
  for(C in c(100,500)) {
    n0<-floor(C*w/cost); ns<-n0
    repeat {
      eligible<-which(w>0 & cost<=C-sum(ns*cost)+1e-9)
      if(!length(eligible)) break
      j<-eligible[which.max(C*w[eligible]-ns[eligible]*cost[eligible])]
      ns[j]<-ns[j]+1
    }
    r<-round_gt_design_budget(a,u,theta,C,q,k,opts=list(cVec_Ds=c(1,0,0)),
                             fix_zero_floor=FALSE,n_index=2)
    nl<-numeric(length(u)); nl[r$design_exact[1,]]<-r$design_exact[2,]
    stopifnot(sum(ns*cost)<=C+1e-8,sum(nl*cost)<=C+1e-8,
              loss(fim(nl,C),k)<=loss(fim(ns,C),k)+1e-7)
    bench<-loss(I,k)
    res[[length(res)+1]]<-data.frame(C=C,criterion=k,
      simple=bench/loss(fim(ns,C),k),local=bench/loss(fim(nl,C),k),
      gain=bench/loss(fim(nl,C),k)-bench/loss(fim(ns,C),k),
      simple_spent=sum(ns*cost),local_spent=sum(nl*cost),
      simple_counts=fmt(ns),local_counts=fmt(nl))
    cat(k,C,'done\n')
  }
}
write.csv(do.call(rbind,res),file.path(out,'rounding_R.csv'),row.names=FALSE)
write.csv(do.call(rbind,cert),file.path(out,'certificates_R.csv'),row.names=FALSE)
write.csv(do.call(rbind,weights),file.path(out,'weights_R.csv'),row.names=FALSE)
p<-.022; eps<-.035; M<-15
bounds<-sqrt(eps*(1-eps)/(p^4*c(500*choose(M,2)^2,6632*M*(M-1)^2/4)))
print(bounds); stopifnot(bounds[1]>.16,bounds[2]>.17)
for(p in c(.001,.022,.07,.3)) for(x in 1:150) {
  r<-1-(1-p)^x-x*p*(1-p)^(x-1)
  stopifnot(abs(r-pbinom(1,x,p,lower.tail=FALSE))<1e-12)
}
# A direct comparison of the two displayed feasible specimen allocations.
u<-1:15; theta<-c(.022,.92,.965); f<-gt_huang2020_regressor(theta,0)
F<-t(vapply(u,f,numeric(3)))
nold<-numeric(15); nold[c(1,7,8,15)]<-c(288,288,1,288)
nnew<-numeric(15); nnew[c(1,4,5,15)]<-c(2209,212,274,147)
stopifnot(sum(u*nold)==6632,sum(u*nnew)==6632)
Jold<-crossprod(F,F*nold); Jnew<-crossprod(F,F*nnew)
ratio<-exp((as.numeric(determinant(Jold,TRUE)$modulus)-as.numeric(determinant(Jnew,TRUE)$modulus))/3)
cat('D-efficiency old versus displayed new allocation:',ratio,'\n')
writeLines(c(paste('SE lower bounds:',paste(bounds,collapse=', ')),
             paste('Old/new determinant efficiency:',ratio)),file.path(out,'analytic_checks.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
print(do.call(rbind,res)); print(do.call(rbind,cert))

# Independent coordinate benchmarks and a conditioned D solve for the specimen cap.
f <- gt_huang2020_regressor(theta,1)
Fc <- t(vapply(u,f,numeric(3)))
B <- solve(chol(crossprod(Fc))); FF <- Fc %*% B
fw <- function(x) as.numeric(FF[x,])
a <- calc_Dopt(u,fw,drop_tol=0,ridge=0)
I <- crossprod(FF,FF*a$weights); V <- solve(I)
Ie <- crossprod(FF,FF*(u*nnew/6632))
spec_eff <- exp((determinant(Ie,TRUE)$modulus-determinant(I,TRUE)$modulus)/3)
spec_residual <- max(rowSums((FF%*%V)*FF)/3)-1
coord <- list()
for(j in 1:3) {
  cvec <- diag(3)[,j]; best <- Inf
  for(ix in combn(15,3,simplify=FALSE)) {
    al <- tryCatch(solve(t(F[ix,]),cvec),error=function(e)NULL)
    if(is.null(al)) next
    value <- sum(abs(al))^2
    if(value<best) {best<-value; bx<-ix; bw<-abs(al)/sum(abs(al))}
  }
  Ic <- crossprod(F[bx,],F[bx,]*bw); vc<-solve(Ic,cvec)
  residual <- max((F%*%vc)^2/best)-1
  stopifnot(residual<1e-9)
  coord[[j]] <- data.frame(parameter=j-1,variance_per_assay=best,SE_500=sqrt(best/500),
    max_relative_violation=residual)
}
write.csv(do.call(rbind,coord),file.path(out,'coordinate_checks_R.csv'),row.names=FALSE)
writeLines(c(paste('Specimen-cap D efficiency of displayed allocation:',spec_eff),
             paste('D benchmark relative directional residual:',spec_residual)),
           file.path(out,'specimen_check_R.txt'))

# Regenerate deterministic scenario figures from the author's pinned R scenarios.
scenario <- new.env(parent=globalenv())
for(expr in parse(file.path('inst','paper_reproduce','stage2_regret_comparison.R'))) {
  eval(expr,scenario)
  if(exists('stage2_compare_tbl',envir=scenario,inherits=FALSE)) break
}
d <- scenario$stage2_compare_tbl
write.csv(d,file.path(out,'stage2_scenarios_R.csv'),row.names=FALSE)
plot_scenario <- function(column,label,file) {
  plot <- ggplot2::ggplot(d,ggplot2::aes(x=.data[[column]],y=forcats::fct_rev(label)))+
    ggplot2::geom_col(fill='#28618a',width=.72)+
    ggplot2::labs(x=label,y=NULL)+
    ggplot2::theme_minimal(base_size=11)+
    ggplot2::theme(panel.grid.major.y=ggplot2::element_blank())
  ggplot2::ggsave(file.path(fig_out,file),plot,width=8,height=6.5)
}
plot_scenario('k_chosen','Selected Stage II pool size','stage2_scenario_rebuilt.pdf')
plot_scenario('regret','Additional expected assays','stage2_regret_rebuilt.pdf')
plot_scenario('actual_expected_tests_under_true','Expected assays at the generating parameters','stage2_expected_rebuilt.pdf')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
