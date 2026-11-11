# Main Table 9 and Figure 4; Supplementary Table S7.
# Recompute the author's scenario inputs with gtDesign's response probability.
# This retains the manuscript's ceiling(N/x) + N*pi(x) approximation.
source("R/load_gtdesign.R")
load_gtdesign()
dir.create("results/paper", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
env <- new.env(parent = globalenv())
for (expr in parse("inst/paper_reproduce/stage2_regret_comparison.R")) {
  eval(expr, env)
  if (exists("theta_scenarios", env, inherits = FALSE)) break
}
scenarios <- env$theta_scenarios
sizes <- 1:20
theta <- c(.07, .93, .96)
N <- 10000
tests <- function(th) ceiling(N/sizes) + N*vapply(sizes, gtDesign::gt_huang2020_pi, numeric(1), theta=th)
true_tests <- tests(theta)
oracle <- which.min(true_tests)
rows <- lapply(seq_len(nrow(scenarios)), function(i) {
  r <- scenarios[i,]
  values <- tests(as.numeric(r[c("p0","p1","p2")]))
  k <- which.min(values)
  data.frame(label=r$label,p0=r$p0,p1=r$p1,p2=r$p2,k_chosen=k,
    expected_tests_under_used_theta=values[k],
    actual_expected_tests_under_true=true_tests[k],
    oracle_k=oracle,oracle_expected_tests_under_true=true_tests[oracle],
    regret=true_tests[k]-true_tests[oracle],
    pct_increase_vs_oracle=100*(true_tests[k]/true_tests[oracle]-1))
})
all <- do.call(rbind, rows)
write.csv(all,"results/paper/supp_table_S07.csv",row.names=FALSE)
order <- c(5,4,6,3,8)
grouped <- all[match(order,all$k_chosen),c("k_chosen","regret","pct_increase_vs_oracle","actual_expected_tests_under_true")]
stopifnot(!anyNA(grouped),all(sort(unique(all$k_chosen))==sort(order)))
grouped$pattern <- c("True parameter values or near truth","Mild misspecification",
 "Moderate underestimation","Strong overestimation","Severe underestimation")
write.csv(grouped,"results/paper/main_table_09.csv",row.names=FALSE)
grouped$pattern <- factor(grouped$pattern,levels=rev(grouped$pattern))
plot <- ggplot2::ggplot(grouped,ggplot2::aes(regret,pattern))+
  ggplot2::geom_point(ggplot2::aes(colour=k_chosen==5),size=2.5,show.legend=FALSE)+
  ggplot2::geom_text(ggplot2::aes(x=regret+8,label=sprintf("%.3f%%",pct_increase_vs_oracle)),hjust=0,size=3.8)+
  ggplot2::geom_vline(xintercept=0,linetype=2)+
  ggplot2::scale_colour_manual(values=c("black","red"))+
  ggplot2::scale_x_continuous(limits=c(-10,530),breaks=seq(0,500,100))+
  ggplot2::labs(title="Stage II regret under parameter misspecification",
    x="Extra expected tests relative to the oracle choice",y=NULL)+
  ggplot2::theme_bw(base_size=12)
ggplot2::ggsave("figures/main_figure_04.pdf",plot,width=8,height=6)
capture.output(sessionInfo(),file="results/paper/stage2_sessionInfo.txt")
cat("Wrote Main Table 9, Figure 4 and Supplementary Table S7.\n")
