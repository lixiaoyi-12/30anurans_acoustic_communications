library(ape)
library(phytools)
library(geiger)

# Continuous traits
tree <- read.tree("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/30物种树文件/tre30re.nwk")
tree$root.edge <- NULL

dat <- read.csv("all.csv")
rownames(dat) <- dat$species
dat$Relative_tympanum <- dat$Diameter_of_tympanum / dat$Body_length

traits <- c("Body_length", "Head_length", "Head_width", "Diameter_of_tympanum",
            "Relative_tympanum", "ABR_threshold", "best_frequency",
            "high_frequency", "ND", "DF", "Fmax")

phylo_signal <- function(trait) {
  x <- dat[[trait]]; names(x) <- rownames(dat)
  x <- x[!is.na(x) & x > 0]
  tree_sub <- keep.tip(tree, intersect(tree$tip.label, names(x)))
  x <- log10(x[tree_sub$tip.label])
  K <- phylosig(tree_sub, x, method="K", test=TRUE, nsim=9999)
  L <- phylosig(tree_sub, x, method="lambda", test=TRUE)
  data.frame(Trait=trait, N=length(x), K=K$K, K_P=K$P,
             Lambda=L$lambda, Lambda_P=L$P)
}

results <- do.call(rbind, lapply(traits, phylo_signal))
results$K_FDR <- p.adjust(results$K_P, method="BH")
results$Lambda_FDR <- p.adjust(results$Lambda_P, method="BH")

write.csv(results, "phylogenetic_signal_results.csv", row.names=FALSE)


# Categorical traits

dat <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/all_categorical_trait.csv")
rownames(dat) <- dat$species
traits <- c("calling_site", "call_type", "size_sac")

discrete_signal <- function(trait) {
  x <- dat[[trait]]; names(x) <- rownames(dat)
  x <- factor(x[!is.na(x) & x != ""])
  tree_sub <- keep.tip(tree, intersect(tree$tip.label, names(x)))
  x <- droplevels(x[tree_sub$tip.label])
  ER <- fitDiscrete(tree_sub, x, model="ER")
  SYM <- fitDiscrete(tree_sub, x, model="SYM")
  ARD <- fitDiscrete(tree_sub, x, model="ARD")
  aicc <- c(ER=ER$opt$aicc, SYM=SYM$opt$aicc, ARD=ARD$opt$aicc)
  model <- names(which.min(aicc))
  L <- fitDiscrete(tree_sub, x, model=model, transform="lambda", bounds=list(lambda=c(0,1)))
  W <- fitDiscrete(tree_sub, x, model=model, transform="white")
  data.frame(Trait=trait, N=length(x), States=nlevels(x), Model=model,
             Lambda=L$opt$lambda, AICc_lambda=L$opt$aicc,
             AICc_white=W$opt$aicc, Delta_AICc=W$opt$aicc-L$opt$aicc)}

results_discrete <- do.call(rbind, lapply(traits, discrete_signal))

write.csv(results_discrete, "discrete_phylogenetic_signal.csv", row.names=FALSE)