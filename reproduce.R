# Execute the command associated with a current manuscript Table/Figure number.
# Usage: Rscript reproduce.R list
#        Rscript reproduce.R main-table-06
#        Rscript reproduce.R all
# "all" rebuilds displays from the saved simulation, not 120000 new MLE fits.
args<-commandArgs(TRUE)
item<-if(length(args))args[1]else"list"
index<-read.csv("paper_items.csv",stringsAsFactors=FALSE,na.strings=character())
if(item=="list"){
 print(index[c("id","paper_item","coverage")],row.names=FALSE)
}else{
 selected<-if(item=="all")index else index[index$id==item,,drop=FALSE]
 if(!nrow(selected))stop("Unknown item. Run: Rscript reproduce.R list")
 executed<-character()
 for(i in seq_len(nrow(selected))){
  r<-selected[i,]
  cat("\n",r$paper_item," [",r$coverage,"]\n",r$notes,"\n",sep="")
  if(!nzchar(r$command)){cat("Source: ",r$script,"\n",sep="");next}
  if(r$command %in% executed){cat("Shared computation already completed in this run.\n");next}
  tokens<-strsplit(r$command," ",fixed=TRUE)[[1]]
  stopifnot(tokens[1]=="Rscript")
  status<-system2(file.path(R.home("bin"),"Rscript"),vapply(tokens[-1],shQuote,character(1)))
  if(status!=0)stop("Failed to reproduce ",r$id," (exit ",status,")")
  executed<-c(executed,r$command)
 }
 cat("\nCompleted ",length(executed)," distinct commands for ",nrow(selected)," indexed items.\n",sep="")
}
