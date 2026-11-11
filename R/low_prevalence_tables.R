# Main Table 2: compute designs and coordinate efficiencies using gtDesign.
# Main Tables 3-4: VERIFY the displayed counts, resources and information.
# The latter two tables are inputs to the verification, not regenerated counts.
source("R/load_gtdesign.R");load_gtdesign()
dir.create("results/paper",recursive=TRUE,showWarnings=FALSE)
theta<-c(.022,.92,.965);u<-1:15
f<-gtDesign::gt_huang2020_regressor(theta,0)
F<-t(vapply(u,f,numeric(3)))
# This change of basis preserves c-optimality when the contrast changes too.
# It avoids an inaccurate solver status for the unconditioned c problem.
B<-solve(chol(crossprod(F)));Fw<-F%*%B;fw<-function(x)as.numeric(Fw[x,])
fits<-list(D=gtDesign::calc_Dopt(u,f,drop_tol=1e-5),A=gtDesign::calc_Aopt(u,f,drop_tol=1e-5),
 E=gtDesign::calc_Eopt(u,f,drop_tol=1e-5,ridge=0,
   tol_gap_abs=1e-12,tol_gap_rel=1e-12,tol_feas=1e-12),
 c=gtDesign::calc_copt(u,fw,cVec=as.numeric(t(B)%*%c(0,1,1)),drop_tol=1e-5,
   tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10),
 Ds=gtDesign::calc_copt(u,f,cVec=c(1,0,0),drop_tol=1e-5))
# The change of basis improves conditioning and transforms each contrast too.
bench<-vapply(1:3,function(j){
 a<-gtDesign::calc_copt(u,fw,cVec=as.numeric(t(B)%*%diag(3)[,j]),drop_tol=0)
 a$value
},numeric(1))
rows<-lapply(names(fits),function(k){
 a<-fits[[k]];I<-crossprod(F,F*as.numeric(a$weights));ef<-bench/diag(solve(I))
 d<-a$design
 data.frame(criterion=k,solver_status=a$status,pool_sizes=paste(d$point,collapse=";"),
   weights=paste(format(d$weight,digits=10),collapse=";"),
   prevalence_efficiency=ef[1],sensitivity_efficiency=ef[2],specificity_efficiency=ef[3])
})
stopifnot(all(vapply(fits,function(a)a$status=="optimal",logical(1))))
write.csv(do.call(rbind,rows),"results/paper/main_table_02.csv",row.names=FALSE)
write.csv(data.frame(parameter=c("prevalence","sensitivity","specificity"),variance=bench),
 "results/paper/main_table_02_coordinate_benchmarks.csv",row.names=FALSE)
displayed<-list(
 `03`=list(D=rbind(c(1,4,5,15),c(2209,212,274,147)),
   A=rbind(c(1,4,15),c(870,818,166)),E=rbind(c(1,4,15),c(870,818,166)),
   c=rbind(c(1,4,15),c(859,817,167)),Ds=rbind(c(1,4,15),c(937,850,153))),
 `04`=list(D=rbind(c(1,7,15),c(167,167,166)),
   A=rbind(c(1,7,15),c(80,258,162)),E=rbind(c(1,7,15),c(80,258,162)),
   c=rbind(c(1,7,15),c(78,259,163)),Ds=rbind(c(1,7,15),c(86,263,151))))
for(id in names(displayed)){
 rows<-lapply(names(displayed[[id]]),function(k){
  a<-displayed[[id]][[k]];x<-a[1,];n<-a[2,];H<-t(vapply(x,f,numeric(3)))
  J<-crossprod(H,H*n)
  stopifnot(if(id=="03")sum(x*n)==6632 else sum(n)==500,qr(J)$rank==3)
  data.frame(criterion=k,counts=paste(paste0(x,":",n),collapse=";"),
   assays=sum(n),specimens=sum(x*n),rank=qr(J)$rank,
   SE_prevalence=sqrt(solve(J)[1,1]),SE_sensitivity=sqrt(solve(J)[2,2]),
   SE_specificity=sqrt(solve(J)[3,3]))
 })
 write.csv(do.call(rbind,rows),paste0("results/paper/main_table_",id,"_verification.csv"),row.names=FALSE)
}
capture.output(sessionInfo(),file="results/paper/low_prevalence_sessionInfo.txt")
cat("Recomputed Main Table 2; verified displayed Main Tables 3-4.\n")
