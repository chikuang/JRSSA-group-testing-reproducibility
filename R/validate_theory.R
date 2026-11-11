# Main Table 1 and Supplementary Section S13.1-S13.2.
# Outputs nonidentifiability_example.csv, jensen_bias.csv, three_size_comparison.csv.
args<-commandArgs(TRUE)
# Load the exact definitions used by the simulation, without running its loop.
e<-new.env(parent=globalenv())
for(expr in parse('R/bias_simulation.R')) {
 eval(expr,e)
 if(exists('round_uniform',envir=e,inherits=FALSE))break
}
with(e,{
 th<-c(.022,.92,.965); x<-c(1,7,15)
 stopifnot(max(abs(jac(x,th)-numDeriv::jacobian(function(t)prob(x,t),th)))<1e-7)
 # Expected score is zero, but this does not establish finite-sample unbiasedness.
 y<-c(10,55,90);n<-c(100,200,300);fit<-fit_mle(y,n,x)
 cat('Estimator check:',fit,'\n')
 # Identical distributions at two distinct parameter vectors with two sizes.
 xx<-c(1,15);v<-prob(xx,th);p2<-.025;z<-(1-p2)^xx
 A<-(v[2]-v[1])/(z[1]-z[2]);s<-v[1]+A*z[1];alternative<-c(p2,s,1+A-s)
 cat('Two-size nonidentifiability: theta=',th,'alternative=',alternative,'\n')
 cat('Maximum probability difference=',max(abs(prob(xx,th)-prob(xx,alternative))),'\n')
 stopifnot(all(alternative[2:3]>.5),all(alternative[2:3]<1),max(abs(prob(xx,th)-prob(xx,alternative)))<1e-12)
 write.csv(data.frame(parameter=c('prevalence','sensitivity','specificity'),truth=th,alternative=alternative),file.path(out,'nonidentifiability_example.csv'),row.names=FALSE)
 # Jensen's inequality: exact positive finite-sample bias in perfect-test submodel.
 m<-10;p<-.022;pr<-1-(1-p)^m
 checks<-do.call(rbind,lapply(c(100,500,2000),function(n){
  yy<-0:n;E<-sum((1-(1-yy/n)^(1/m))*dbinom(yy,n,pr))
  leading<-(m-1)/(2*m^2*n)*pr*(1-pr)^(1/m-1)
  data.frame(n=n,p=p,m=m,exact_bias=E-p,leading_bias=leading)
 }))
 stopifnot(all(checks$exact_bias>0));print(checks)
 write.csv(checks,file.path(out,'jensen_bias.csv'),row.names=FALSE)
})
# Same assay budget, contrasting structural identification and information.
with(e,{
 sets<-list(one=10,two=c(1,15),close=c(1,2,3),spread=c(1,7,15),four=c(1,6,10,15));th<-c(.022,.92,.965)
 ans<-lapply(names(sets),function(nm){
  x<-sets[[nm]];n<-switch(as.character(length(x)),`1`=500,`2`=c(250,250),
    `3`=c(167,167,166),`4`=rep(125,4))
  J<-fisher(x,n,th);rank<-qr(jac(x,th))$rank
  V<-if(rank==3)solve(J)else matrix(NA_real_,3,3)
  data.frame(design=nm,pool_sizes=paste(x,collapse=';'),counts=paste(n,collapse=';'),
   assays=sum(n),rank=qr(jac(x,th))$rank,SE_p=sqrt(V[1,1]),
   SE_sensitivity=sqrt(V[2,2]),SE_specificity=sqrt(V[3,3]),
   logdet=if(rank==3)as.numeric(determinant(J,TRUE)$modulus)else -Inf)
 })
 ans<-do.call(rbind,ans);ans$relative_D_efficiency<-exp((ans$logdet-ans$logdet[ans$design=='spread'])/3)
 stopifnot(identical(as.integer(ans$rank),c(1L,2L,3L,3L,3L)),all(ans$assays==500))
 write.csv(ans,file.path(out,'three_size_comparison.csv'),row.names=FALSE)
})
