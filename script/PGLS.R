library(nlme)

tree <- read.tree("tre30re.nwk")
data <- read.csv("all.csv", header = TRUE, row.names = 1)
clean_data <- na.omit(data[, c("best_frequency", "DF","high_frequency","Fmax")])

valid_species <- rownames(clean_data)
pruned_tree <- keep.tip(tree, valid_species)
tree <- pruned_tree
clean_data$species <- rownames(clean_data)

clean_data$log_DF <- log10(clean_data$df)
clean_data$log_abr <- log10(clean_data$abr)
clean_data$log_Fmax <- log10(clean_data$Fmax)
clean_data$log_HF <- log10(clean_data$high_frequency)

bm_model <- corBrownian(value = 1, phy = tree, form = ~species)
pgls_DF_BF <- gls(log_BF ~ log_DF, data = clean_data, correlation = bm_model, method = "ML")
pgls_HF_Fmax <- gls(log_HF ~ log_Fmax, data = clean_data, correlation = bm_model, method = "ML")

get_pgls_line <- function(fit, x, panel) {newx <- seq(min(x), max(x), length.out = 200)
  X <- cbind(1, newx)
  beta <- coef(fit)
  V <- vcov(fit)
  fitted <- as.numeric(X %*% beta)
  se <- sqrt(diag(X %*% V %*% t(X)))
  data.frame(x = newx, fit = fitted, lower = fitted - 1.96 * se, upper = fitted + 1.96 * se, panel = panel)
}

plot_data <- rbind(
  data.frame(x = clean_data$log_DF, y = clean_data$log_BF, panel = "DF ~ BF"),
  data.frame(x = clean_data$log_Fmax, y = clean_data$log_HF, panel = "Fmax ~ HF"))

line_data <- rbind(get_pgls_line(pgls_DF_BF, clean_data$log_DF, "DF ~ BF"),
  get_pgls_line(pgls_HF_Fmax, clean_data$log_Fmax, "Fmax ~ HF"))

result_DF_BF <- summary(pgls_DF_BF)$tTable
result_HF_Fmax <- summary(pgls_HF_Fmax)$tTable

beta_DF_BF <- result_DF_BF["log_DF", "Value"]
p_DF_BF <- result_DF_BF["log_DF", "p-value"]
beta_HF_Fmax <- result_HF_Fmax["log_Fmax", "Value"]
p_HF_Fmax <- result_HF_Fmax["log_Fmax", "p-value"]

format_p <- function(p) {if (p < 0.001) "P < 0.001" else paste0("P = ", sprintf("%.3f", p))}

label_data <- data.frame(panel = c("DF ~ BF", "Fmax ~ HF"),
  label = c(paste0("β = ", sprintf("%.2f", beta_DF_BF), "\n", format_p(p_DF_BF)), paste0("β = ", sprintf("%.2f", beta_HF_Fmax), "\n", format_p(p_HF_Fmax))))

panel_order <- c("DF ~ BF", "Fmax ~ HF")

plot_data$panel <- factor(plot_data$panel, levels = panel_order)
line_data$panel <- factor(line_data$panel, levels = panel_order)
label_data$panel <- factor(label_data$panel, levels = panel_order)

p <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(data = line_data, aes(x = x, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.18) +
  geom_line(data = line_data, aes(x = x, y = fit), inherit.aes = FALSE, linewidth = 0.9) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(data = label_data, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.08, vjust = 1.15, size = 4) +
  facet_wrap(~panel, nrow = 1, scales = "free")

