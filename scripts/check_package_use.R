# Check the package-based execution path without replacing published simulations.
# Usage: Rscript scripts/check_package_use.R vendor/gtDesign results/checks 2 1
e <- new.env(parent=globalenv())
for(expr in parse("R/bias_simulation.R")) {
  eval(expr,e)
  if(exists("round_uniform",e,inherits=FALSE))break
}
for(p in c(.00001,.022,.07,.5)) for(s in c(.51,.92,1)) for(t in c(.51,.965,1)) {
  th<-c(p,s,t);x<-1:150
  pr<-vapply(x,gtDesign::gt_huang2020_pi,numeric(1),theta=th)
  gr<-t(vapply(x,gtDesign::gt_huang2020_f,numeric(3),theta=th))
  stopifnot(max(abs(e$prob(x,th)-pr))<1e-12,max(abs(e$jac(x,th)-gr))<1e-10)
}
# A public package function rejects accuracy = 0.5; the estimator intentionally
# includes this boundary and must retain its vectorized likelihood implementation.
stopifnot(all(is.finite(e$prob(c(1,7,15),c(.022,.5,.965)))))
files<-list.files("results/bias",pattern="^cell_[0-9]+[.]rds$",full.names=TRUE)
stopifnot(length(files)==12L)
for(file in files) {
  a<-readRDS(file);cc<-a$cell
  assign(".Random.seed",a$seed,envir=.GlobalEnv)
  yy<-sapply(seq_along(cc$x),function(j)rbinom(nrow(a$y),cc$n[j],
    gtDesign::gt_huang2020_pi(cc$x[j],cc$theta)))
  stopifnot(identical(yy,a$y))
  for(i in 1:2) {
    fit<-e$fit_mle(a$y[i,],cc$n,cc$x)
    stopifnot(abs(fit["nll"]-a$estimates$nll[i])<1e-7)
  }
}
cat("Package model checks passed; all 120000 generated outcome vectors unchanged.\n")
cat("First two joint MLEs in each of 12 cells retain the saved likelihood values.\n")
