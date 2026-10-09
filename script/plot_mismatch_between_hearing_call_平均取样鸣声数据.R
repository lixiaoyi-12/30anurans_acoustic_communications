library(ggplot2)

hearing_file <- "~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/听力数据/Rhacophorus_rhodopus_hearing.txt"
call_file <- "~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/鸣声数据/all/Rhacophorus_rhodopus/Rhacophorus_rhodopus_summary_average.txt"

x1 <- sub("\t+$", "", readLines(hearing_file))
x2 <- sub("\t+$", "", readLines(call_file))
hearing <- read.delim(text = x1, header = TRUE, check.names = FALSE)
call <- read.delim(text = x2, header = TRUE, check.names = FALSE)

# ABR
threshold_cols <- grep("^threshold_A", names(hearing), value = TRUE)
hearing[threshold_cols] <- lapply(hearing[threshold_cols], function(x) suppressWarnings(as.numeric(as.character(x))))
hearing$frequency <- as.numeric(hearing$frequency)
m <- as.matrix(hearing[, threshold_cols])

abr <- data.frame(frequency = hearing$frequency, mean = rowMeans(m, na.rm = TRUE), se = apply(m, 1, function(x) {x <- x[!is.na(x)]; if(length(x) > 1) sd(x)/sqrt(length(x)) else 0}))
abr <- abr[is.finite(abr$mean), ]
abr$lower <- abr$mean - abr$se
abr$upper <- abr$mean + abr$se

# Call
call_cols <- grep("^threshold_A", names(call), value = TRUE)
call[call_cols] <- lapply(call[call_cols], function(x) suppressWarnings(as.numeric(as.character(x))))
call$frequency <- as.numeric(call$frequency)
m2 <- as.matrix(call[, call_cols])

call_mean <- data.frame(
  frequency = call$frequency,
  mean = rowMeans(m2, na.rm = TRUE),
  se = apply(m2, 1, function(x) {x <- x[!is.na(x)]; if(length(x) > 1) sd(x) / sqrt(length(x)) else 0})
)
call_mean <- call_mean[is.finite(call_mean$mean), ]
call_mean$lower <- call_mean$mean - call_mean$se
call_mean$upper <- call_mean$mean + call_mean$se

# smooth
smooth_curve <- function(x, y, xout) spline(x, y, xout = xout, method = "natural")$y
xf_abr <- seq(min(abr$frequency), max(abr$frequency), by = 0.01)
xf_call <- seq(min(call_mean$frequency), max(call_mean$frequency), by = 0.01)

abr_s <- data.frame(frequency = xf_abr)
abr_s$mean <- smooth_curve(abr$frequency, abr$mean, xf_abr)
abr_s$lower <- smooth_curve(abr$frequency, abr$lower, xf_abr)
abr_s$upper <- smooth_curve(abr$frequency, abr$upper, xf_abr)

call_s <- data.frame(frequency = xf_call)
call_s$mean <- smooth_curve(call_mean$frequency, call_mean$mean, xf_call)
call_s$lower <- smooth_curve(call_mean$frequency, call_mean$lower, xf_call)
call_s$upper <- smooth_curve(call_mean$frequency, call_mean$upper, xf_call)

left_min <- floor(min(abr_s$lower, na.rm = TRUE)/10)*10
left_max <- ceiling(max(abr_s$upper, na.rm = TRUE)/10)*10

right_min <- floor(min(call_s$lower, na.rm = TRUE)/10)*10
right_max <- ceiling(max(call_s$upper, na.rm = TRUE)/10)*10

a <- (left_max-left_min)/(right_min-right_max)
b <- left_min-a*right_max

call_s$y <- a*call_s$mean+b
call_s$ymin <- a*call_s$upper+b
call_s$ymax <- a*call_s$lower+b

p <- ggplot() +
  geom_ribbon(data = abr_s, aes(frequency, ymin = lower, ymax = upper), fill = "#3A9D5D", alpha = 0.2) +
  geom_line(data = abr_s, aes(frequency, mean), color = "#3A9D5D", linewidth = 0.8) +
  geom_point(data = abr, aes(frequency, mean), color = "#3A9D5D", size = 2) +
  geom_ribbon(data = call_s, aes(frequency, ymin = ymin, ymax = ymax), fill = "grey50", alpha = 0.2) +
  geom_line(data = call_s, aes(frequency, y), color = "grey50", linewidth = 0.8) +
  geom_point(data = call_mean, aes(frequency, a*mean+b), color = "grey50", shape = 21, fill = "white", size = 2) +
  scale_x_continuous(limits = c(0.2, 8)) +
  scale_y_reverse("ABR threshold (dB)", limits = c(left_max, left_min), sec.axis = sec_axis(~ (. - b)/a, name = "Relative Amplitude (dB)")) +
  labs(x = "Frequency (kHz)") +
  theme_classic(base_size = 14)

p
