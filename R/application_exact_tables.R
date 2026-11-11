# Main Tables 3 and 4: construct every allocation with gtDesign.
# Usage: Rscript R/application_exact_tables.R
# No manuscript allocation is an input. The specimen problem uses q=1, C=6632;
# the assay problem uses q=0, C=500. Original scripts remain in inst/.
source('R/load_gtdesign.R')
load_gtdesign()

application_setup <- function(q) {
  u <- 1:15
  theta <- c(.022,.92,.965)
  f <- gtDesign::gt_huang2020_regressor(theta,q)
  F <- t(vapply(u,f,numeric(3)))
  B <- solve(chol(crossprod(F)))
  Q <- svd(F)$v
  Fw <- F %*% B
  Fr <- F %*% Q
  fw <- function(x) as.numeric(Fw[x,])
  fr <- function(x) as.numeric(Fr[x,])
  # D is invariant to nonsingular congruence; c needs B^T c as its contrast.
  # A uses only an orthogonal rotation, preserving the original trace criterion.
  # E is solved on the original parameter scale, with tighter absolute accuracy.
  fits <- list(
    D=gtDesign::calc_Dopt(u,fw,ridge=0,drop_tol=1e-6,
      tol_gap_abs=1e-11,tol_gap_rel=1e-11,tol_feas=1e-11),
    A=gtDesign::calc_Aopt(u,fr,drop_tol=1e-6,
      tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9),
    E=gtDesign::calc_Eopt(u,f,ridge=0,drop_tol=1e-6,
      tol_gap_abs=1e-12,tol_gap_rel=1e-12,tol_feas=1e-12),
    c=gtDesign::calc_copt(u,fw,cVec=as.numeric(t(B)%*%c(0,1,1)),drop_tol=1e-6,
      tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10),
    Ds=gtDesign::calc_copt(u,fw,cVec=as.numeric(t(B)%*%c(1,0,0)),drop_tol=1e-6,
      tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10))
  list(u=u,theta=theta,q=q,f=f,F=F,B=B,Q=Q,fits=fits,
    opts=list(cVec_c=c(0,1,1),cVec_Ds=c(1,0,0)))
}

application_certificate <- function(F,w,k,opts) {
  I <- crossprod(F,F*w)
  V <- solve(I)
  G <- switch(k,D=V/3,A=V%*%V,
    E={eg<-eigen(I,symmetric=TRUE);tcrossprod(eg$vectors[,3])},
    c=tcrossprod(V%*%opts$cVec_c),Ds=tcrossprod(V[,1]))
  scores <- rowSums((F%*%G)*F)
  max(scores/sum(diag(G%*%I)))-1
}

# Enumerate the package's candidate completions, retaining its fixed floor.
# For D with 500 assays, numerical ties use the fewest specimens, then the
# lexicographically largest small-pool counts. The loss tolerance is 1e-9.
application_D_tie <- function(r,u,tolerance=1e-9) {
  ext <- as.data.frame(r$extension_table,check.names=FALSE)
  cols <- intersect(as.character(u),names(ext))
  floor <- numeric(length(u));floor[match(r$design_round1[1,],u)]<-r$design_round1[2,]
  candidates <- matrix(floor,nrow=nrow(ext),ncol=length(u),byrow=TRUE)
  candidates[,match(as.integer(cols),u)] <- candidates[,match(as.integer(cols),u),drop=FALSE]+as.matrix(ext[,cols,drop=FALSE])
  eligible <- which(ext$loss <= min(ext$loss)+tolerance)
  specimens <- as.numeric(candidates%*%u)
  order_args <- c(list(specimens[eligible]),lapply(seq_along(u),function(j)-candidates[eligible,j]))
  pick <- eligible[do.call(order,order_args)[1]]
  n <- candidates[pick,];keep<-n>0
  list(design=rbind(x=u[keep],n_prime=n[keep]),
    audit=data.frame(counts=apply(candidates[eligible,,drop=FALSE],1,function(z)
      paste(paste0(u[z>0],':',z[z>0]),collapse=';')),
      specimens=specimens[eligible],loss=ext$loss[eligible],
      selected=eligible==pick),tolerance=tolerance)
}

run_application_tables <- function(dest='results/paper') {
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  runs<-list();certificates<-list();weights<-list()
  for(q in c(1,0)) {
    setup<-application_setup(q);C<-if(q==1)6632 else 500
    rows<-list();id<-if(q==1)'03'else'04'
    for(k in names(setup$fits)) {
      a<-setup$fits[[k]]
      residual<-application_certificate(setup$F,a$weights,k,setup$opts)
      stopifnot(a$status=='optimal',is.finite(residual),residual<1e-4)
      r<-gtDesign::round_gt_design_budget(a,setup$u,setup$theta,C,q,k,
        opts=setup$opts,n_index=2,fix_zero_floor=FALSE)
      de<-r$design_exact
      if(q==0 && k=='D') {
        tie<-application_D_tie(r,setup$u);de<-tie$design
        write.csv(tie$audit,file.path(dest,'main_table_04_D_ties.csv'),row.names=FALSE)
      }
      x<-de[1,];n<-de[2,]
      H<-t(vapply(x,gtDesign::gt_huang2020_regressor(setup$theta,0),numeric(3)))
      J<-crossprod(H,H*n);M<-J/C
      I<-crossprod(setup$F,setup$F*a$weights)
      # Public package criterion comparison, evaluated on the original scale.
      refs<-list(D=-as.numeric(determinant(I,logarithm=TRUE)$modulus),
        A=sum(diag(solve(I))),E=-min(eigen(I,symmetric=TRUE)$values),
        c=drop(t(setup$opts$cVec_c)%*%solve(I,setup$opts$cVec_c)),Ds=solve(I)[1,1])
      eff<-gtDesign::exact_design_efficiency_maximin(M,refs[k],k,setup$opts)$efficiencies[[k]]
      direct_eff<-switch(k,
        D=exp((as.numeric(determinant(M,TRUE)$modulus)+refs$D)/3),
        A=refs$A/sum(diag(solve(M))),
        E=-min(eigen(M,symmetric=TRUE)$values)/refs$E,
        c=refs$c/drop(t(setup$opts$cVec_c)%*%solve(M,setup$opts$cVec_c)),
        Ds=refs$Ds/solve(M)[1,1])
      stopifnot(abs(eff-direct_eff)<1e-8)
      spent<-sum(n*(1-q+q*x));V<-solve(J)
      stopifnot(all(n==round(n)),all(n>0),spent<=C+1e-8,
        qr(H)$rank==3,eff>0.999 && eff<=1+1e-6)
      rows[[k]]<-data.frame(criterion=k,counts=paste(paste0(x,':',n),collapse=';'),
        assays=sum(n),specimens=sum(x*n),budget=C,q=q,spent=spent,remaining=C-spent,
        efficiency=eff,SE_prevalence=sqrt(V[1,1]),SE_sensitivity=sqrt(V[2,2]),
        SE_specificity=sqrt(V[3,3]),solver_status=a$status,max_relative_violation=residual)
      certificates[[paste(q,k)]]<-data.frame(table=paste0('Main Table ',as.integer(id)),
        criterion=k,q=q,solver_status=a$status,max_relative_violation=residual)
      weights[[paste(q,k)]]<-data.frame(table=paste0('Main Table ',as.integer(id)),
        criterion=k,q=q,pool=setup$u,weight=a$weights)
      runs[[paste(q,k)]]<-list(setup=setup[c('u','theta','q','B','Q')],approx=a,
        rounding=r,selected_design=de)
    }
    tab<-do.call(rbind,rows);rownames(tab)<-NULL
    write.csv(tab,file.path(dest,paste0('main_table_',id,'.csv')),row.names=FALSE)
    tex<-vapply(seq_len(nrow(tab)),function(i){
      r<-tab[i,];cr<-if(r$criterion=='Ds')'D_s'else r$criterion
      paste0('$',cr,'$ & ',gsub(':',' (',gsub(';','); ',r$counts,fixed=TRUE),fixed=TRUE),') & ',
        r$assays,' & ',r$specimens,' & ',sprintf('%.6f',r$efficiency),' \\\\')
    },character(1))
    writeLines(tex,file.path(dest,paste0('main_table_',id,'_rows.tex')))
    print(tab[c('criterion','counts','assays','specimens','efficiency')],row.names=FALSE)
  }
  write.csv(do.call(rbind,certificates),file.path(dest,'application_certificates.csv'),row.names=FALSE)
  write.csv(do.call(rbind,weights),file.path(dest,'application_weights.csv'),row.names=FALSE)
  saveRDS(runs,file.path(dest,'application_designs.rds'))
  writeLines(trimws(capture.output(sessionInfo()),which='right'),
    file.path(dest,'application_sessionInfo.txt'))
  invisible(runs)
}
if(sys.nframe()==0)run_application_tables()
