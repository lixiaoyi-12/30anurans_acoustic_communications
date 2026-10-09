library(ggplot2)
library(dplyr)

dat <- read.csv("~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/听力数据/听力数据汇总-6species.csv", na.strings = c("", "NA", "无"))

threshold_cols <- grep("^threshold_", names(dat), value = TRUE)

dat[threshold_cols] <- lapply(
  dat[threshold_cols],
  function(x) as.numeric(as.character(x))
)
dat$n <- rowSums(!is.na(dat[threshold_cols]))
dat$mean_threshold <- rowMeans(dat[threshold_cols], na.rm = TRUE)
dat$mean_threshold[dat$n == 0] <- NA
dat$sd <- apply(dat[threshold_cols], 1, sd, na.rm = TRUE)
dat$se <- dat$sd / sqrt(dat$n)

smooth_list <- lapply(split(dat, dat$Species), function(x) {
  x <- x[order(x$frequency), ]
  m <- x[!is.na(x$mean_threshold), ]
  if (nrow(m) < 2)
    return(NULL)
  xnew <- seq(min(m$frequency), max(m$frequency), length.out = 500)
  mean_smooth <- spline(m$frequency, m$mean_threshold, xout = xnew, method = "natural"
  )$y
  s <- x[!is.na(x$se), ]
  if (nrow(s) >= 2) {se_smooth <- spline(s$frequency, s$se, xout = xnew, method = "natural"
    )$y
    se_smooth <- pmax(se_smooth, 0)
  } else {
    se_smooth <- rep(NA, length(xnew))
  }
  data.frame(
    Species = m$Species[1],
    frequency = xnew,
    mean = mean_smooth,
    lower = mean_smooth - se_smooth,
    upper = mean_smooth + se_smooth
  )
})

smooth_dat <- bind_rows(smooth_list)

ggplot() +
  geom_ribbon(data = smooth_dat, aes(x = frequency, ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(data = smooth_dat, aes(x = frequency,y = mean),linewidth = 0.8) +
  geom_point(data = dat %>% filter(!is.na(mean_threshold)), aes(x = frequency, y = mean_threshold), size = 1.5) + 
  facet_wrap(~ Species, axes = "all", axis.labels = "all", ncol = 4) + theme_bw() +
  theme(strip.background = element_rect(fill = NA), strip.text = element_text(face = "italic"))

