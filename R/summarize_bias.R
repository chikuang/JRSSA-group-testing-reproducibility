# Main Table 10; Supplementary Tables S10-S13 and Figures S5-S6.
# Uses saved simulation results; does not rerun or change the estimator.
# Build the publication tables and scientific figures from the saved R results.
d<-read.csv('results/bias/summary.csv');a<-read.csv('results/bias/allocations.csv')
z<-read.csv('results/bias/diagnostics.csv');rep<-read.csv('results/bias/replicates.csv')
dir.create('figures',showWarnings=FALSE);dir.create('docs',showWarnings=FALSE)
colors<-c(D='#185C89',U4='#B34C26');symbols<-c(D=16,U4=17)
for(ctx in c('cost','low')){
 pdf(paste0('figures/bias_',ctx,'.pdf'),width=9,height=8.2)
 layout(rbind(c(1,1,1),matrix(2:10,nrow=3,byrow=TRUE),c(11,11,11)),heights=c(.65,2.3,2.3,2.3,.4))
 par(mar=c(0,0,0,0));plot.new()
 text(.5,.77,if(ctx=='cost')'Cost-constrained study: prevalence 0.07, sensitivity 0.93, specificity 0.96'
       else 'Low-prevalence study: prevalence 0.022, sensitivity 0.92, specificity 0.965',cex=.95)
 legend('center',inset=0,legend=c('Rounded D-optimal','Uniform four-size'),col=colors,pch=symbols,horiz=TRUE,bty='n',cex=.9)
 par(mar=c(3.5,4.1,2.5,1),mgp=c(2.3,.7,0),cex=.85)
 for(param in c('prevalence','sensitivity','specificity'))for(metric in c('bias','variance','rmse')){
  dd<-d[d$context==ctx & d$parameter==param,];mult<-if(metric=='variance')1e4 else 100
  mc<-dd[[paste0(metric,'_mcse')]]*mult;vv<-dd[[metric]]*mult
  lim<-range(vv-1.96*mc,vv+1.96*mc,if(metric=='bias')0 else NA,na.rm=TRUE)
  plot(range(dd$C),lim,type='n',log='x',xaxt='n',xlab=if(ctx=='cost')'Cost budget C' else 'Assay budget N',
       ylab=if(metric=='variance')'Percentage points squared' else 'Percentage points',
       main=paste(tools::toTitleCase(param),switch(metric,bias='bias',variance='variance',rmse='RMSE')))
  axis(1,at=sort(unique(dd$C)),labels=format(sort(unique(dd$C)),big.mark=','))
  if(metric=='bias')abline(h=0,col='grey60',lty=3)
  for(des in c('D','U4')){
   ss<-dd$design==des;xx<-dd$C[ss]*if(des=='D').98 else 1.02
   lines(xx,vv[ss],col=colors[des],pch=symbols[des],type='b',lwd=1.5)
   lo<-vv[ss]-1.96*mc[ss];hi<-vv[ss]+1.96*mc[ss]
   segments(xx,lo,xx,hi,col=colors[des]);segments(xx/1.025,lo,xx*1.025,lo,col=colors[des]);segments(xx/1.025,hi,xx*1.025,hi,col=colors[des])
  }
 }
 par(mar=c(0,0,0,0));plot.new();text(.5,.5,sprintf('%s independent replicates per cell; bars show +/- 1.96 Monte Carlo SE.',format(unique(d$replicates)[1],big.mark=',')),cex=.8)
 dev.off()
}
# Point estimates and Monte Carlo SE, all in percentage points (variance in pp^2).
rows<-function(dd,full=FALSE){
 out<-character()
 for(i in seq_len(nrow(dd))){r<-dd[i,]
  vals<-c(if(full)format(r$C,scientific=FALSE)else if(r$context=='cost')'Cost' else 'Low',
   tools::toTitleCase(r$parameter),r$design,
   sprintf('%.3f (%.3f)',100*r$bias,100*r$bias_mcse),
   if(full)sprintf('%.3f (%.3f)',1e4*r$variance,1e4*r$variance_mcse)else sprintf('%.3f',1e4*r$variance),
   sprintf('%.3f (%.3f)',100*r$rmse,100*r$rmse_mcse))
  out<-c(out,paste0(paste(vals,collapse=' & '),' \\\\'))
 }
 out
}
writeLines(rows(d[d$C==500,]),'docs/main_bias_rows.tex')
for(ctx in c('cost','low'))writeLines(rows(d[d$context==ctx,],TRUE),paste0('docs/bias_',ctx,'_rows.tex'))
alloc<-unique(a[c('context','C','design','assays','specimens','D_eff','condition_number')])
alloc$counts<-vapply(seq_len(nrow(alloc)),function(i){r<-alloc[i,];ss<-a$context==r$context&a$C==r$C&a$design==r$design;paste(paste0(a$x[ss],':',a$n[ss]),collapse='; ')},character(1))
write.csv(alloc,'results/bias/allocation_summary.csv',row.names=FALSE)
writeLines(apply(alloc,1,function(r)paste0(paste(c(if(r['context']=='cost')'Cost' else 'Low',r['C'],r['design'],r['counts'],r['assays'],r['specimens'],sprintf('%.4f',as.numeric(r['D_eff']))),collapse=' & '),' \\\\')),'docs/allocation_rows.tex')
rep<-merge(rep,alloc[c('context','C','design','assays')]);rep$normalized_pg<-rep$pg/rep$assays
pg<-aggregate(normalized_pg~context+C+design,rep,max);z<-z[,!grepl('^normalized_pg',names(z)),drop=FALSE];z<-merge(z,pg)
write.csv(z,'results/bias/diagnostics.csv',row.names=FALSE)
writeLines(apply(z,1,function(r)paste0(paste(c(if(r['context']=='cost')'Cost' else 'Low',r['C'],r['design'],r['nonzero_optimizer_code'],sprintf('%.2f',100*as.numeric(r['any_accuracy_boundary'])),sprintf('%.1e',as.numeric(r['normalized_pg']))),collapse=' & '),' \\\\')),'docs/diagnostic_rows.tex')
# An independent-stream comparison, useful for interpreting Monte Carlo uncertainty.
D<-d[d$design=='D',];U<-d[d$design=='U4',];key<-c('context','C','parameter')
cmp<-merge(D,U,by=key,suffixes=c('_D','_U4'))
cmp$bias_difference<-cmp$bias_D-cmp$bias_U4
cmp$bias_difference_mcse<-sqrt(cmp$bias_mcse_D^2+cmp$bias_mcse_U4^2)
cmp$rmse_difference<-cmp$rmse_D-cmp$rmse_U4
cmp$rmse_difference_mcse<-sqrt(cmp$rmse_mcse_D^2+cmp$rmse_mcse_U4^2)
write.csv(cmp,'results/bias/design_comparison.csv',row.names=FALSE)
cat('Wrote tables, allocation summary, diagnostics and bias figures.\n')
