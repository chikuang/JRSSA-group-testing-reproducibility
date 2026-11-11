# Recompute all parameter blocks for Main Tables 6-7 and Supplementary Table S3.
# Usage: Rscript R/paper_designs.R [main-table-06|main-table-07|supp-table-S03]
# gtDesign supplies the model, optimizers and criterion efficiencies. Timings and
# near-tied supports can differ by solver; full weights are saved for inspection.
source("R/load_gtdesign.R")
load_gtdesign()
args <- commandArgs(TRUE)
item <- if(length(args))args[1]else"main-table-06"
stopifnot(item %in% c("main-table-06","main-table-07","supp-table-S03"))
dest <- file.path("results/paper",item)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
theta <- c(.07,.93,.96)
opts <- list(cVec_Ds=c(1,0,0),cVec_c=c(0,1,1))
rows <- list(); weights <- list()
for (M in c(61,150)) for(q in if(item=="supp-table-S03")c(0,.2,.6,.8)else c(0,.2,.8)) {
  u <- seq_len(M)
  f <- gtDesign::gt_huang2020_regressor(theta,q)
  fits <- list(D=gtDesign::calc_Dopt(u,f,drop_tol=1e-6),
    A=gtDesign::calc_Aopt(u,f,drop_tol=1e-6),
    Ds=gtDesign::calc_copt(u,f,cVec=opts$cVec_Ds,drop_tol=1e-6),
    c=gtDesign::calc_copt(u,f,cVec=opts$cVec_c,drop_tol=1e-6),
    E=gtDesign::calc_Eopt(u,f,drop_tol=1e-6))
  refs <- list(D=-fits$D$value,A=fits$A$value,Ds=fits$Ds$value,c=fits$c$value,E=fits$E$value)
  criteria <- if(item=="main-table-06") as.list(names(fits)) else if(item=="main-table-07")
    list(c("D","A"),c("D","A","Ds"),c("D","A","c","E")) else list(c("D","Ds"),c("D","Ds","c"))
  for(cr in criteria) {
    tm <- proc.time()[3]
    a <- if(length(cr)==1)fits[[cr]]else gtDesign::compute_maximin_design(u,f,
      loss_ref=refs[cr],criteria=cr,opts=opts,drop_tol=1e-6)
    elapsed <- proc.time()[3]-tm
    key <- paste(cr,collapse="-")
    # Evaluate the full unrounded weight vector using the package regressors.
    F <- t(vapply(u,f,numeric(3)))
    I <- crossprod(F,F*as.numeric(a$weights))
    eff <- gtDesign::exact_design_efficiency_maximin(I,refs,names(fits),opts)$efficiencies
    row <- data.frame(M=M,q=q,criterion=key,solver_status=a$status,minimum_efficiency=min(eff[cr]),
      elapsed_seconds=if(length(cr)==1)NA_real_ else elapsed)
    for(k in names(eff)) row[[paste0("eff_",k)]] <- eff[k]
    if(length(cr)==1)row$loss <- switch(cr,D=exp(-fits$D$value),E=fits$E$value,a$value)
    rows[[length(rows)+1]] <- row
    weights[[length(weights)+1]] <- data.frame(M=M,q=q,criterion=key,pool=u,weight=a$weights)
  }
}
write.csv(do.call(rbind,rows),file.path(dest,"criteria.csv"),row.names=FALSE)
write.csv(do.call(rbind,weights),file.path(dest,"weights.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
cat("Recomputed",item,"using gtDesign; outputs:",dest,"\n")
