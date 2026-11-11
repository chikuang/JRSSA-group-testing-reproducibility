# Audit original inst/paper_reproduce/table6.R and table7.R against JRSSA_rev2.tex.
# Reference counts are comparison data only, never inputs to an optimizer.
# Run after Rscript run.R original and Rscript R/application_exact_tables.R.
source('R/application_exact_tables.R')
dir.create('results/paper/historical_application',recursive=TRUE,showWarnings=FALSE)
dest<-'results/paper/historical_application'
ref<-list(
 `6`=c(D='1:288;7:288;8:1;15:288',A='1:121;5:1;7:396;14:1;15:248',
 E='1:121;5:1;7:396;14:1;15:248',c='1:118;7:397;15:249',Ds='1:136;7:418;15:238'),
 `7`=c(D='1:167;7:167;15:166',A='1:80;7:258;15:162',
 E='1:80;7:258;15:162',c='1:78;7:259;15:163',Ds='1:86;7:263;15:151'))
serialize_counts<-function(de)paste(paste0(de[1,],':',de[2,]),collapse=';')
parse_counts<-function(x){a<-strsplit(x,';',fixed=TRUE)[[1]];do.call(cbind,lapply(strsplit(a,':',fixed=TRUE),as.numeric))}
f<-gtDesign::gt_huang2020_regressor(c(.022,.92,.965),0)
info<-function(s){a<-parse_counts(s);H<-t(vapply(a[1,],f,numeric(3)));crossprod(H,H*a[2,])}
opts<-list(cVec_c=c(0,1,1),cVec_Ds=c(1,0,0))
loss<-function(M,k)switch(k,D=-as.numeric(determinant(M,TRUE)$modulus),A=sum(diag(solve(M))),
 E=-min(eigen(M,symmetric=TRUE)$values),c=drop(t(opts$cVec_c)%*%solve(M,opts$cVec_c)),Ds=solve(M)[1,1])
relative<-function(s,r,k){l<-loss(info(s),k);lr<-loss(info(r),k);switch(k,D=exp((lr-l)/3),E=l/lr,lr/l)}
rows<-list()
for(table in c('6','7')){
 original<-read.csv(paste0('results/original/table',table,'/table',table,'.csv'),check.names=FALSE)
 e<-new.env(parent=globalenv())
 # The script is untouched; this diagnostic changes only E solver accuracy.
 e$calc_Eopt<-function(u,f,...){gtDesign::calc_Eopt(u,f,...,ridge=0,
  tol_gap_abs=1e-12,tol_gap_rel=1e-12,tol_feas=1e-12)}
 printed<-capture.output(sys.source(paste0('inst/paper_reproduce/table',table,'.R'),envir=e))
 writeLines(printed,file.path(dest,paste0('table',table,'_tight_E_console.txt')))
 tight<-get(paste0('table',table),e)
 write.csv(tight,file.path(dest,paste0('table',table,'_tight_E.csv')),row.names=FALSE)
 for(k in names(ref[[table]])){
  default_counts<-gsub(', ',';',original[original$Criterion==k,'Exact design (pool:count)'],fixed=TRUE)
  tight_counts<-gsub(', ',';',tight[tight$Criterion==k,'Exact design (pool:count)'],fixed=TRUE)
  adopted<-tight_counts
  if(table=='7' && k=='D'){
   a<-gtDesign::calc_Dopt(1:15,f,drop_tol=1e-5)
   r<-gtDesign::round_gt_design_budget(a,1:15,c(.022,.92,.965),500,0,'D')
   adopted<-serialize_counts(application_D_tie(r,1:15)$design)
  }
  reference<-ref[[table]][[k]]
  rows[[paste(table,k)]]<-data.frame(source_script=paste0("table",table,".R"),criterion=k,reference_counts=reference,
    default_counts=default_counts,tight_E_counts=tight_counts,with_tie_policy_counts=adopted,
    default_exact_match=identical(default_counts,reference),
    controlled_exact_match=identical(adopted,reference),
    default_efficiency_relative_to_reference=relative(default_counts,reference,k),
    controlled_efficiency_relative_to_reference=relative(adopted,reference,k))
 }
}
out<-do.call(rbind,rows);rownames(out)<-NULL
stopifnot(all(out$controlled_exact_match))
write.csv(out,file.path(dest,'original_manuscript_comparison.csv'),row.names=FALSE)
# Save the small default E weight and its zero-floor effect explicitly.
es<-list();er<-list()
for(mode in c('default','tight')){
 a<-if(mode=='default')gtDesign::calc_Eopt(1:15,f,drop_tol=1e-5)else
   gtDesign::calc_Eopt(1:15,f,drop_tol=1e-5,ridge=0,
      tol_gap_abs=1e-12,tol_gap_rel=1e-12,tol_feas=1e-12)
 es[[mode]]<-data.frame(mode=mode,pool=1:15,weight=a$weights)
 for(tab in c('6','7')){
 r<-if(tab=='6')gtDesign::round_gt_design_subject_budget(a,1:15,c(.022,.92,.965),6632,0,'E')else
  gtDesign::round_gt_design_budget(a,1:15,c(.022,.92,.965),500,0,'E')
 er[[paste(mode,tab)]]<-data.frame(mode=mode,source_script=paste0("table",tab,".R"),
  floor_after_zero_fix=serialize_counts(r$design_round1),final_counts=serialize_counts(r$design_exact))
 }
}
write.csv(do.call(rbind,es),file.path(dest,'E_solver_weights.csv'),row.names=FALSE)
write.csv(do.call(rbind,er),file.path(dest,'E_rounding_trace.csv'),row.names=FALSE)
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(dest,'sessionInfo.txt'))
cat('Original application tables (table6.R, table7.R): ',sum(out$default_exact_match),'/10 counts match at defaults; ',
 sum(out$controlled_exact_match),'/10 match with tight E precision and the declared D tie rule.\n',sep='')
