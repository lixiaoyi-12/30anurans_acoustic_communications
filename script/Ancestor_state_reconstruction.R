library(phytools)
tree <- read.tree("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/30物种树文件/tre30re.nwk")
plotTree(tree, type="fan",ftype="i")

#continuous trait
svl <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/svl.csv", header=FALSE)
td <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/td.csv", header=FALSE)
nd <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/nd.csv")
Fmax <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/Fmax.csv")
Fmin <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/Fmin.csv")
df <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/DF.csv")
abr <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/ABR.csv")
hf <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/HF.csv")
bf <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/BF.csv")

svl <- setNames(svl[, 2], svl[, 1])
td <- setNames(td[, 2], td[, 1])


# 验证：物种名是否和树的tip匹配
all(tree$tip.label %in% names(td))  # 应返回 TRUE

# 然后正常使用
fit <- fastAnc(tree, svl)
obj <- contMap(tree, svl)
plot(obj, type="phylogram", legend=0.7*max(nodeHeights(tree)), fsize=c(0.7, 0.9))
nodelabels(text=round(fit, 2), node=(Ntip(tree)+1):(Ntip(tree)+tree$Nnode), frame="none", cex=0.6)
plot(obj, type="fan", legend=0.7*max(nodeHeights(tree)), fsize=c(0.7, 0.9))
nodelabels(text=round(anc, 2), node=(Ntip(tree)+1):(Ntip(tree)+tree$Nnode), frame="none", cex=0.6)

fit_td <- fastAnc(tree, td)
obj_td <- contMap(tree, td, plot=FALSE)
plot(obj_td, type="phylogram", legend=0.7*max(nodeHeights(tree)), fsize=c(0.7, 0.9))
nodelabels(text=round(fit_td, 2), node=(Ntip(tree)+1):(Ntip(tree)+tree$Nnode), frame="none", cex=0.6)

p1_phenogram <- phenogram(tree, svl, fsize=0.6, spread.costs=c(1,0))
phenogram(tree, svl, fsize=0.6, spread.costs=c(0.2,1))
p_td__phenogram <- phenogram(tree, td, fsize=0.6, spread.costs=c(1,0))
p1_phenogram+p_td__phenogram

node_est <- fit$ace
plot(obj, type="fan", legend=0.7*max(nodeHeights(tree)), fsize=c(0.7, 0.9))
nodelabels(text=round(node_est, 2), cex=0.4, frame="none")
node_est <- fit_td$ace
plot(obj_td, type="fan", legend=0.7*max(nodeHeights(tree)), fsize=c(0.7, 0.9))
nodelabels(text=round(node_est, 2), cex=0.4, frame="none")

# missing data
nd$species <- trimws(nd$species)
trait_vec <- setNames(nd[,2], nd$species)

nd <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/nd.csv")
bf$species <- trimws(bf$species)
trait_vec <- setNames(bf[,2], bf$species)

# 找到NA物种
na_species <- names(trait_vec)[is.na(trait_vec)]
cat("性状缺失物种：", na_species,"\n")

# 保留有数据的物种
keep_sp <- names(trait_vec)[!is.na(trait_vec)]
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, keep_sp))
trait_clean <- trait_vec[keep_sp]

# 正常运行无报错
fit <- fastAnc(tree_pruned, trait_clean, vars=TRUE, CI=TRUE)
obj <- contMap(tree_pruned, trait_clean, plot=FALSE)
plot(obj, type="fan")
phenogram(tree_pruned, trait_clean)

library(ape)
library(phytools)

# =========================
# 1. 读取树和数据
# =========================
run_discrete_asr <- function(trait, dat, tree) {
  x <- dat[[trait]]
  names(x) <- rownames(dat)
  x <- x[!is.na(x) & x != ""]
  
  sp <- intersect(tree$tip.label, names(x))
  tree_sub <- keep.tip(tree, sp)
  tree_sub$root.edge <- NULL
  
  x <- factor(x[tree_sub$tip.label])
  
  fit_ER <- fitMk(tree_sub, x, model="ER", pi="fitzjohn")
  fit_SYM <- fitMk(tree_sub, x, model="SYM", pi="fitzjohn")
  fit_ARD <- fitMk(tree_sub, x, model="ARD", pi="fitzjohn")
  
  fits <- list(ER=fit_ER, SYM=fit_SYM, ARD=fit_ARD)
  aic <- sapply(fits, AIC)
  best_model <- names(which.min(aic))
  best_fit <- fits[[best_model]]
  
  anc <- ancr(best_fit)
  
  cat("\nTrait:", trait, "\n")
  cat("N:", length(x), "\n")
  cat("States:", paste(levels(x), collapse=", "), "\n")
  print(aic)
  cat("Best model:", best_model, "\n")
  
  return(list(tree=tree_sub, x=x, fits=fits, AIC=aic, best_model=best_model, best_fit=best_fit, anc=anc))
}

call_type_asr <- run_discrete_asr("call_type", dat, tree)
calling_site_asr <- run_discrete_asr("calling_site", dat, tree)
size_sac_asr <- run_discrete_asr("size_sac", dat, tree)
