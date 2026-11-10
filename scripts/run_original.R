# Run one preserved source script in an isolated R process.
args<-commandArgs(TRUE)
if(length(args)!=1)stop("Usage: Rscript scripts/run_original.R inst/paper_reproduce/table1.R")
file<-normalizePath(args[1]); root<-normalizePath(".")
.libPaths(c(file.path(root,".library"),.libPaths()))
for(f in list.files("vendor/gtDesign/R",full.names=TRUE,pattern="[.]R$"))source(f)
name<-tools::file_path_sans_ext(basename(file));dest<-file.path(root,"results/original",name)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
setwd(dest)
pdf("plots.pdf",width=8,height=6)
sink("console.txt",split=TRUE)
env<-new.env(parent=globalenv())
if(basename(file)=="stage2_regret_comparison.R"){
  for(expr in parse(file)){
    if(is.call(expr) && identical(expr[[1]],as.name("<-")) && identical(expr[[2]],as.name("p_regret_grouped"))){
      cat("NOTE: preserved source refers to undefined regret_group_tbl in its trailing plot.\n")
      cat("The complete scenario calculations and preceding plots were reproduced.\n")
      cat("Use R/revision_checks.R for the reproducible scenario/regret figures.\n")
      break
    }
    eval(expr,envir=env)
  }
}else sys.source(file,envir=env)
for(nm in ls(env)){
 obj<-get(nm,env)
 if(is.data.frame(obj))write.csv(obj,paste0(nm,".csv"),row.names=FALSE)
 if(inherits(obj,"ggplot"))ggplot2::ggsave(paste0(nm,".pdf"),obj,width=8,height=6)
}
capture.output(sessionInfo(),file="sessionInfo.txt")
sink();dev.off()
