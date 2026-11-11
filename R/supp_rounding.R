# Supplementary Tables S4 and S6: package-based integer design calculations.
# Usage: Rscript R/supp_rounding.R [supp-table-S04|supp-table-S06]
# Full numerical results are saved. Near ties and solver tolerances can select
# alternative integer supports; this is not a transcription of printed counts.
source("R/load_gtdesign.R");load_gtdesign()
args<-commandArgs(TRUE);item<-args[1]
stopifnot(item %in% c("supp-table-S04","supp-table-S06"))
M<-if(item=="supp-table-S04")150L else 61L
q<-if(item=="supp-table-S04")0 else .2
budgets<-if(q==0)c(5:9,50)else c(100,500,10000)
u<-seq_len(M);theta<-c(.07,.93,.96)
f<-gtDesign::gt_huang2020_regressor(theta,q)
opts<-list(cVec_Ds=c(1,0,0),cVec_c=c(0,1,1))
fits<-list(D=gtDesign::calc_Dopt(u,f,drop_tol=1e-6),A=gtDesign::calc_Aopt(u,f,drop_tol=1e-6),
 Ds=gtDesign::calc_copt(u,f,cVec=opts$cVec_Ds,drop_tol=1e-6),
 c=gtDesign::calc_copt(u,f,cVec=opts$cVec_c,drop_tol=1e-6))
if(q>0)fits$E<-gtDesign::calc_Eopt(u,f,drop_tol=1e-6)
rows<-list()
for(k in names(fits))for(C in budgets){
 a<-gtDesign::round_gt_design_budget(fits[[k]],u,theta,C,q,k,opts=opts,n_index=2,fix_zero_floor=FALSE)
 x<-a$design_exact[1,];n<-a$design_exact[2,]
 cost<-vapply(x,gtDesign::gt_huang2020_cost,numeric(1),q=q)
 stopifnot(all(n==round(n)),sum(n*cost)<=C+1e-8)
 rows[[length(rows)+1]]<-data.frame(criterion=k,budget=C,
  counts=paste(paste0(x,":",n),collapse=";"),spent=sum(n*cost),
  remaining=C-sum(n*cost),loss=a$loss_exact,efficiency=a$efficiency)
}
dest<-file.path("results/paper",item);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,rows),file.path(dest,"allocations.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
cat("Recomputed",item,"with gtDesign rounding functions.\n")
