library(ape)
library(phylolm)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

data_file <- "all.xlsx"
tree_file <- "tre30re.nwk"
output_dir <- "phylolm_call_results_XL"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- as.data.frame(read_excel(data_file))
dat <- dat %>% select(species, Body_length, Head_length, Head_width, Diameter_of_tympanum,
  Relative_tympanum = `Diameter_of_tympanum/body_length`, ABR_threshold, best_frequency, high_frequency,
  call_type, VC, calling_site, size_sac, DF, Fmax, ND)
numeric_columns <- c("Body_length", "Head_length", "Head_width", "Diameter_of_tympanum",
                     "Relative_tympanum", "ABR_threshold", "best_frequency",
                     "high_frequency", "VC", "DF", "Fmax", "ND")
dat[numeric_columns] <- lapply(dat[numeric_columns], as.numeric)

#setting reference trait types
dat$calling_site <- relevel(factor(dat$calling_site), ref = "arboreal")
dat$call_type <- relevel(factor(dat$call_type), ref = "complex")
dat$size_sac <- relevel(factor(dat$size_sac), ref = "external")

#input tree
tree <- read.tree(tree_file)
tree$root.edge <- NULL
dat <- dat[dat$species %in% tree$tip.label, ]

z <- function(x) as.numeric(scale(log10(x)))
prepare_data <- function(response, continuous = NULL, factor = NULL) {
  vars <- c("species", response, continuous, factor)
  d <- dat %>% select(all_of(vars)) %>% drop_na()
  tree_sub <- drop.tip(tree, setdiff(tree$tip.label, d$species))
  d <- d[match(tree_sub$tip.label, d$species), ]
  rownames(d) <- d$species
  d[[paste0(response, "_z")]] <- z(d[[response]])
  for(x in continuous) d[[paste0(x, "_z")]] <- z(d[[x]])
  list(data = d, tree = tree_sub)
}

extract_result <- function(fit, response, model, predictor = NA) {
  b <- fit$coefficients
  se <- sqrt(diag(fit$vcov))
  df <- fit$n - length(b)
  t <- b / se
  p <- 2 * pt(abs(t), df = df, lower.tail = FALSE)
  ci <- qt(0.975, df)
  data.frame(response = response, model = model, predictor = predictor, term = names(b), estimate = b, std_error = se, t_value = t, p_value = p, conf_low = b - ci * se, conf_high = b + ci * se, N = fit$n, lambda = fit$optpar, AIC = fit$aic)}
response_labels <- c(ND = "Note duration", DF = "Dominant frequency", Fmax = "Maximum frequency")
term_labels <- c(
  Body_length_z = "SVL",
  Head_length_z = "HL",
  Head_width_z = "HW",
  Diameter_of_tympanum_z = "TD",
  Relative_tympanum_z = "TD/SVL",
  ABR_threshold_z = "auditory threshold",
  best_frequency_z = "BF",
  high_frequency_z = "HF",
  VC_z = "VC",
  call_typetonal = "Call type: tonal",
  call_typetrill = "Call type: trill",
  `calling_siteterrestrial` = "terrestrial",
  `calling_sitestill water` = "still water",
  `calling_siteflowing water` = "flowing water",
  size_sacinternal = "Vocal sac: internal"
)
term <- function(x) ifelse(x %in% names(term_labels), term_labels[x], x)
responses <- c("ND", "DF", "Fmax")

# univariate models
continuous_predictors <- c("SVL", "HL", "HW", "TD", "TD/SVL", "VC", "auditory_threshold", "BF", "HF")

factor_predictors <- c("call_type", "calling_site", "size_sac")
results <- list()
fits <- list()
k <- 1

for(response in responses) {
  for(predictor in continuous_predictors) {
    x <- prepare_data(response, predictor)
    fit <- phylolm(as.formula(paste0(response, "_z ~ ", predictor, "_z")), x$data, x$tree, model = "lambda")
    fits[[paste(response, predictor, sep = "_")]] <- fit
    results[[k]] <- extract_result(fit, response, paste(response, "~", predictor), predictor)
    k <- k + 1}
  for(predictor in factor_predictors) {
    x <- prepare_data(response, factor = predictor)
    fit <- phylolm(as.formula(paste0(response, "_z ~ ", predictor)), x$data, x$tree, model = "lambda")
    fits[[paste(response, predictor, sep = "_")]] <- fit
    results[[k]] <- extract_result(fit, response, paste(response, "~", predictor), predictor)
    k <- k + 1}}

univariate <- bind_rows(results)

write.csv(univariate, file.path(output_dir, "call_univariate_all_coefficients.csv"), row.names = FALSE)
saveRDS(fits, file.path(output_dir, "call_univariate_all_models.rds"))

# forest plot of univariate models
univariate_plot <- univariate %>%
  filter(term != "(Intercept)") %>%
  mutate(
    response_label = response_labels[response],
    term_label = term(term),
    group = case_when(
      predictor %in% c("SVL", "HL", "HW", "TD", "TD/SVL") ~ "Morphology",
      predictor %in% c("auditory_threshold", "BF", "HF") ~ "Hearing",
      predictor == "call_type" ~ "Call type",
      predictor %in% c("calling_site", "VC") ~ "Ecology",
      predictor == "size_sac" ~ "Vocal sac"),
    significance = ifelse(p_value < 0.05, "P < 0.05", "P ≥ 0.05")) %>%
  group_by(response) %>%
  mutate(p_adjust_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

p1 <- ggplot(univariate_plot, aes(estimate, term_label, colour = group, shape = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = term_label)) +
  geom_point(size = 3) +
  facet_wrap(~response_label, nrow = 1) +
  scale_shape_manual(values = c("P < 0.05" = 17, "P ≥ 0.05" = 16)) +
  labs(x = "Standardized effect size", y = NULL, colour = "Predictor", shape = "Significance") +
  theme_classic()

ggsave(file.path(output_dir, "call_Univariate_plot.pdf"), p1, width = 18, height = 8)


# multivariate models
model_specs <- expand.grid(base = c("SVL", "TD"), hearing = c("auditory_threshold", "BF", "HF"),stringsAsFactors = FALSE)

results <- list()
fits <- list()
k <- 1

for(response in responses) {
  for(i in 1:nrow(model_specs)) {
    base <- model_specs$base[i]
    hearing <- model_specs$hearing[i]
    x <- prepare_data(response, c(base, hearing), "calling_site")
    formula <- as.formula(paste0(response, "_z ~ ", base, "_z + calling_site + ", hearing, "_z"))
    fit <- phylolm(formula, x$data, x$tree, model = "lambda")
    model_label <- paste(base, "+ calling site +", hearing)
    fits[[paste(response, base, hearing, sep = "_")]] <- fit
    result <- extract_result(fit, response, model_label)
    result$base_predictor <- base
    result$hearing_predictor <- hearing
    result$model_label <- model_label
    results[[k]] <- result
    k <- k + 1}}

multivariable <- bind_rows(results)
write.csv(multivariable, file.path(output_dir, "call_Multivariable_models.csv"), row.names = FALSE)
saveRDS(fits, file.path(output_dir, "call_Multivariable_models.rds"))

# forest plot of multivariate models
multivariable_plot <- multivariable %>%
  filter(term != "(Intercept)") %>%
  mutate(
    response_label = response_labels[response],
    term_label = term(term),
    group = case_when(
      grepl("SVL|TD", term) ~ "Morphology",
      grepl("auditory_threshold|BF|HF", term) ~ "Hearing",
      grepl("^calling_site", term) ~ "Ecology"
    ),
    significance = ifelse(p_value < 0.05, "P < 0.05", "P ≥ 0.05"),
    panel_label = paste(response_label, model_label, sep = "\n")
  )

p2 <- ggplot(multivariable_plot, aes(estimate, term_label, colour = group, shape = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = term_label)) +
  geom_point(size = 2.8) +
  facet_wrap(~panel_label, ncol = 6, scales = "free_y") +
  scale_shape_manual(values = c("P < 0.05" = 17, "P ≥ 0.05" = 16)) +
  theme_classic()

# Body length + call type + calling site + hearing
hearing_predictors <- c("auditory_threshold", "BF", "HF")
results <- list()
fits <- list()
k <- 1

for(response in responses) {
  for(hearing in hearing_predictors) {
    x <- prepare_data(response, c("SVL", hearing), c("call_type", "calling_site"))
    formula <- as.formula(paste0(
      response, "_z ~ SVL_z + call_type + calling_site + ", hearing, "_z"))
    fit <- phylolm(formula, x$data, x$tree, model = "lambda")
    model_label <- paste("SVL + call type + calling site +", hearing)
    fits[[paste(response, hearing, sep = "_")]] <- fit
    result <- extract_result(fit, response, model_label)
    result$hearing_predictor <- hearing
    result$model_label <- model_label
    results[[k]] <- result
    k <- k + 1}}

calltype_models <- bind_rows(results)

calltype_plot <- calltype_models %>%
  filter(term != "(Intercept)") %>%
  mutate(
    response_label = response_labels[response],
    term_label = term(term),
    group = case_when(
      term == "SVL_z" ~ "Morphology",
      grepl("^calling_site", term) ~ "Ecology",
      grepl("^call_type", term) ~ "Call type",
      grepl("auditory_threshold|BF|HF", term) ~ "Hearing"
    ),
    significance = ifelse(p_value < 0.05, "P < 0.05", "P ≥ 0.05"),
    panel_label = paste(response_label, hearing_predictor, sep = "\n")
  )

p3 <- ggplot(calltype_plot, aes(estimate, term_label, colour = group, shape = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = term_label)) +
  geom_point(size = 2.8) +
  facet_wrap(~panel_label, ncol = 3, scales = "free_y") +
  scale_shape_manual(values = c("P < 0.05" = 17, "P ≥ 0.05" = 16)) +
  theme_classic()
