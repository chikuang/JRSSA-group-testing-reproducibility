# Entrypoint: Rscript run.R [bias|revision|original|all] [replicates=10000] [cores=4]
args<-commandArgs(TRUE);mode<-if(length(args))args[1]else"bias"
B<-if(length(args)>1)args[2]else"10000";cores<-if(length(args)>2)args[3]else"4"
stopifnot(mode %in% c("bias","revision","original","all"))
run<-function(script,args=character()){
 code<-system2(file.path(R.home("bin"),"Rscript"),c(script,args))
 if(code!=0)stop("Failed: ",script," (exit ",code,")")
}
if(mode %in% c("bias","all")){
 dir.create("results/bias",recursive=TRUE,showWarnings=FALSE)
 run("R/bias_simulation.R",c("vendor/gtDesign","results/bias",B,cores))
 run("R/validate_theory.R",c("vendor/gtDesign","results/bias",B,cores))
 run("R/profile_audit.R",c("vendor/gtDesign","results/bias",B,cores))
 run("R/summarize_bias.R")
}
if(mode %in% c("revision","all")){
 dir.create("results/revision",recursive=TRUE,showWarnings=FALSE)
 dir.create("figures",showWarnings=FALSE)
 # Source path passed to the original check is resolved before its source calls.
 run("R/revision_checks.R",c("vendor/gtDesign","results/revision","figures"))
}
if(mode %in% c("original","all")){
 files<-c(list.files("inst/paper_reproduce",full.names=TRUE,pattern="[.]R$"),
          list.files("inst/examples",full.names=TRUE,pattern="[.]R$"))
 for(f in files)run("scripts/run_original.R",f)
}
