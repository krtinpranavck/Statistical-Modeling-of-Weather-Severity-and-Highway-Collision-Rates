# Member 9: Integrated Statistical Model (Section 9)
# Project : Statistical Modeling of Weather Severity and Highway Collision Rates
# Usage   : Rscript member_9_integrated_model.R path/to/UK_Accident.csv output_folder
#           (or open in RStudio and click Source: the default paths below point to
#           the Downloads folder, so no setwd() is needed)
# Needs   : base R and the recommended MASS package only (same as Member 6).
# Source files are never modified. All counts describe reported injury accidents.
#
# What this script does
#   Part A  Load the file and remove exact duplicates (same rule as Members 4 and 6).
#   Part B  Build the predictors used in Table 9.1.
#   Part C  Severity model: logistic regression, two outcomes (KSI, Fatal),
#           two models each (Adj-A without hour of day, Adj-B with hour of day).
#           Outputs: odds ratios, fit (AIC, McFadden R2, AUC), likelihood-ratio
#           tests, collinearity (GVIF).
#   Part D  Frequency model: negative binomial regression on daily counts.
#           Outputs: model comparison (log-likelihood, AIC, BIC) and rate ratios.
#   Part E  Figures 9.1 and 9.2 and a text log.
#
# Values from the Python run that this script should reproduce (Section 9):
#   Severity sample 1,408,977 records; 3,286 days; 1,469,991 records after de-duplication.
#   KSI Adj-B odds ratios: dark lit 1.167, dark unlit or none 1.062, rain 0.821,
#   frost or ice 0.554, speed 60 mph 1.609, 1 vehicle 1.942, hour 00-03 1.318.
#   Daily model: log-likelihood -18,345.1 (calendar) and -18,220.5 (full),
#   LR chi-square 249.1 on 3 df, rain/snow/fog share IRR 1.028, darkness share IRR 0.871.

args  <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "/Users/jaisatyamalla/Downloads/UK_Accident.csv"
out   <- if (length(args) >= 2) args[2] else "/Users/jaisatyamalla/Downloads/member9_results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"),  showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)
options(stringsAsFactors = FALSE, scipen = 7, width = 110)
invisible(Sys.setlocale("LC_TIME", "C"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("Install the recommended MASS package.")

LOG <- character(0)
note <- function(...) { s <- paste0(...); cat(s, sep = "\n"); LOG <<- c(LOG, s) }
show <- function(x, label = NULL) {
  s <- capture.output(print(x))
  if (!is.null(label)) s <- c(label, s)
  cat(s, sep = "\n"); LOG <<- c(LOG, s)
}
write_tab <- function(x, name) write.csv(x, file.path(out, "tables", paste0(name, ".csv")), row.names = FALSE, na = "")
z95 <- qnorm(0.975)


# ============================ PART A. LOAD DATA ==============================
note("PART A. Loading ", input)
header <- names(read.csv(input, nrows = 0, check.names = FALSE))
needed <- c("Accident_Index", "Location_Easting_OSGR", "Location_Northing_OSGR",
            "Accident_Severity", "Number_of_Vehicles", "Number_of_Casualties", "Date", "Time",
            "Speed_limit", "Road_Type", "Light_Conditions", "Weather_Conditions",
            "Road_Surface_Conditions", "Urban_or_Rural_Area", "Year")
stopifnot(all(needed %in% header))
classes <- ifelse(header %in% needed, "character", "NULL")     # read only the columns we need
raw <- read.csv(input, colClasses = classes, check.names = FALSE, na.strings = c("", "NA"))
raw_n <- nrow(raw)

# Exact duplicate rule used by Members 4 and 6 (same values in all of these fields).
key_cols  <- c("Year", "Accident_Index", "Date", "Time", "Location_Easting_OSGR",
               "Location_Northing_OSGR", "Accident_Severity", "Number_of_Vehicles", "Number_of_Casualties")
duplicate <- duplicated(raw[key_cols])
d <- raw[!duplicate, ]
removed_n <- sum(duplicate)
rm(raw, duplicate); gc(verbose = FALSE)
note("Records read: ", format(raw_n, big.mark = ","), "; exact duplicates removed: ",
     format(removed_n, big.mark = ","), "; records kept: ", format(nrow(d), big.mark = ","))
stopifnot(nrow(d) == 1469991)        # the figure used throughout the report


# ====================== PART B. BUILD THE PREDICTORS =========================
d$Year <- as.integer(d$Year)
d$date <- as.Date(d$Date, format = "%d/%m/%Y")
stopifnot(!anyNA(d$date), all(as.integer(format(d$date, "%Y")) == d$Year))

# Time -> hour (a few records have no valid time)
valid_time <- !is.na(d$Time) & grepl("^([01][0-9]|2[0-3]):[0-5][0-9]$", d$Time)
d$hour <- NA_integer_
d$hour[valid_time] <- as.integer(substr(d$Time[valid_time], 1, 2))
d$hour_band <- cut(d$hour, breaks = c(-1, 3, 7, 11, 15, 19, 23),
                   labels = c("00-03", "04-07", "08-11", "12-15", "16-19", "20-23"))
d$hour_band <- relevel(d$hour_band, ref = "12-15")

# Severity outcomes (DfT coding 1 = Fatal, 2 = Serious, 3 = Slight)
sev <- as.integer(d$Accident_Severity)
stopifnot(all(sev %in% 1:3))
d$ksi   <- sev %in% c(1L, 2L)      # fatal or serious
d$fatal <- sev == 1L

# Weather: fine, rain, snow, fog or mist. "Other" and "Unknown" become NA.
w <- tolower(d$Weather_Conditions)
weather <- rep(NA_character_, nrow(d))
weather[grepl("^fine", w)]    <- "Fine"
weather[grepl("^rain", w)]    <- "Rain"
weather[grepl("^snow", w)]    <- "Snow"
weather[grepl("fog|mist", w)] <- "Fog or mist"
d$weather <- factor(weather, levels = c("Fine", "Rain", "Snow", "Fog or mist"))

# Lighting. The file spells one label "Darkeness", so we test for the start "dark".
l    <- tolower(d$Light_Conditions)
dark <- grepl("^dark", l)
light <- rep(NA_character_, nrow(d))
light[grepl("^daylight", l)]                       <- "Daylight"
light[dark & grepl("lights present and lit", l)]   <- "Dark, lit"
light[dark & grepl("unlit|no street", l)]          <- "Dark, unlit or none"
light[dark & grepl("lighting unknown", l)]         <- "Dark, lighting unknown"
stopifnot(!anyNA(light))
d$light <- factor(light, levels = c("Daylight", "Dark, lit", "Dark, unlit or none", "Dark, lighting unknown"))
d$is_dark <- dark

# Road surface. "Normal" (1,945 records) is counted as dry. Provisional until Section 6.3 is added.
surface_map <- c("Dry" = "Dry", "Normal" = "Dry", "Wet/Damp" = "Wet or damp",
                 "Frost/Ice" = "Frost or ice", "Snow" = "Snow on surface",
                 "Flood (Over 3cm of water)" = "Flooded")
surface <- unname(surface_map[d$Road_Surface_Conditions])
stopifnot(!anyNA(surface))
d$surface <- factor(surface, levels = c("Dry", "Wet or damp", "Frost or ice", "Snow on surface", "Flooded"))

# Speed limit (10 and 15 mph records join the "20 or below" group). Reference = 30.
sp <- as.integer(d$Speed_limit)
d$speed <- cut(sp, breaks = c(-Inf, 20, 30, 40, 50, 60, Inf),
               labels = c("20 or below", "30", "40", "50", "60", "70"))
d$speed <- relevel(d$speed, ref = "30")

# Area, road type, vehicles, weekend, year
ar <- as.character(d$Urban_or_Rural_Area)
d$area <- factor(ifelse(ar == "1", "Urban", ifelse(ar == "2", "Rural", NA)), levels = c("Urban", "Rural"))
d$road <- factor(d$Road_Type, levels = c("Single carriageway", "Dual carriageway", "Roundabout",
                                         "One way street", "Slip road", "Unknown"))
stopifnot(!anyNA(d$road))
nv <- as.integer(d$Number_of_Vehicles)
d$vehicles <- factor(ifelse(nv >= 4, "4 or more", as.character(nv)), levels = c("2", "1", "3", "4 or more"))
d$dow <- as.POSIXlt(d$date)$wday                      # 0 = Sunday ... 6 = Saturday
d$weekend <- factor(ifelse(d$dow %in% c(0, 6), "Weekend", "Weekday"), levels = c("Weekday", "Weekend"))
d$year <- factor(d$Year)                              # 2005 is the reference

# Severity-model sample: known weather, valid hour, allocated area
use <- !is.na(d$weather) & !is.na(d$hour_band) & !is.na(d$area)
note("Records dropped: weather other/unknown = ", sum(is.na(d$weather)),
     "; no valid hour = ", sum(is.na(d$hour_band)), "; area not allocated = ", sum(is.na(d$area)))
note("Severity model sample: ", format(sum(use), big.mark = ","), " of ", format(nrow(d), big.mark = ","),
     " (", round(100 * mean(use), 1), "%)")
stopifnot(sum(use) == 1408977)
note("Within the sample: fatal or serious = ", round(100 * mean(d$ksi[use]), 2), "%; fatal = ",
     round(100 * mean(d$fatal[use]), 2), "%")


# ================== PART C. SEVERITY MODEL (LOGISTIC REGRESSION) =============
note("")
note("PART C. Severity model")
# Records with the same value on every predictor are combined into one row with counts.
# This gives the same coefficients as a record-level fit but is much faster.
predictors <- c("light", "weather", "surface", "speed", "area", "road", "vehicles",
                "weekend", "year", "hour_band")
m   <- d[use, c(predictors, "ksi", "fatal")]
key <- do.call(paste, c(lapply(m[predictors], as.integer), sep = "|"))
counts <- rowsum(cbind(n = 1, ksi = as.numeric(m$ksi), fatal = as.numeric(m$fatal)), key, reorder = FALSE)
first  <- !duplicated(key)                            # rowsum(reorder = FALSE) keeps first-appearance order
cells  <- cbind(m[first, predictors], as.data.frame(counts), row.names = NULL)
cells$ksi_no   <- cells$n - cells$ksi
cells$fatal_no <- cells$n - cells$fatal
rm(m, key, counts, first)
note("Combined into ", format(nrow(cells), big.mark = ","), " covariate patterns.")
stopifnot(sum(cells$n) == sum(use))

rhs_A <- "light + weather + surface + speed + area + road + vehicles + weekend + year"
rhs_B <- paste(rhs_A, "+ hour_band")
outcomes <- list(
  "Fatal or serious" = list(yes = "ksi",   no = "ksi_no"),
  "Fatal"            = list(yes = "fatal", no = "fatal_no"))

fit_logit <- function(o, rhs) {
  glm(as.formula(paste0("cbind(", o$yes, ", ", o$no, ") ~ ", rhs)), family = binomial, data = cells)
}

# Record-level log-likelihood from the fitted probabilities (the grouped glm value would
# include a term that depends on how the records were grouped, so we compute it directly).
loglik_record <- function(model, o) {
  p <- fitted(model); k <- cells[[o$yes]]; n <- cells$n
  sum(ifelse(k > 0, k * log(p), 0) + ifelse(n - k > 0, (n - k) * log(1 - p), 0))
}
loglik_null <- function(o) {
  k <- sum(cells[[o$yes]]); n <- sum(cells$n); k * log(k / n) + (n - k) * log(1 - k / n)
}
# Area under the ROC curve from grouped data. Ties get half credit.
auc_grouped <- function(p, pos, neg) {
  ord <- order(p); p <- round(p[ord], 12); pos <- pos[ord]; neg <- neg[ord]
  grp <- cumsum(c(TRUE, diff(p) != 0))
  P <- as.numeric(rowsum(pos, grp)); N <- as.numeric(rowsum(neg, grp))
  below <- cumsum(N) - N
  sum(P * (below + 0.5 * N)) / (sum(P) * sum(N))
}
# Generalised variance inflation factor (Fox and Monette), the same quantity as car::vif.
gvif_table <- function(model) {
  asg <- attr(model.matrix(model), "assign"); labs <- attr(terms(model), "term.labels")
  keep <- asg > 0; R <- cov2cor(vcov(model)[keep, keep]); a <- asg[keep]
  ld <- determinant(R, logarithm = TRUE)$modulus
  do.call(rbind, lapply(seq_along(labs), function(j) {
    inj <- a == j
    g <- exp(determinant(R[inj, inj, drop = FALSE], logarithm = TRUE)$modulus +
             determinant(R[!inj, !inj, drop = FALSE], logarithm = TRUE)$modulus - ld)
    data.frame(term = labs[j], df = sum(inj), gvif = as.numeric(g),
               gvif_adjusted = as.numeric(g)^(1 / (2 * sum(inj))))
  }))
}
# Odds-ratio table for every coefficient (Wald 95% CI)
or_table <- function(model, outcome_name, model_name) {
  asg <- attr(model.matrix(model), "assign"); labs <- c("(Intercept)", attr(terms(model), "term.labels"))
  cm <- coef(summary(model)); variable <- labs[asg + 1]
  data.frame(outcome = outcome_name, model = model_name, variable = variable,
             level = substr(rownames(cm), nchar(variable) + 1, nchar(rownames(cm))),
             odds_ratio = exp(cm[, 1]), lower = exp(cm[, 1] - z95 * cm[, 2]),
             upper = exp(cm[, 1] + z95 * cm[, 2]), p = cm[, 4], row.names = NULL)[asg > 0, ]
}

fit_rows <- list(); or_rows <- list(); lr_rows <- list(); gvif_rows <- list(); MODELS <- list()
for (on in names(outcomes)) {
  o  <- outcomes[[on]]
  mA <- fit_logit(o, rhs_A); mB <- fit_logit(o, rhs_B)
  MODELS[[on]] <- list(A = mA, B = mB)
  ll0 <- loglik_null(o)
  for (nm in c("A", "B")) {
    mod <- MODELS[[on]][[nm]]; ll <- loglik_record(mod, o)
    fit_rows[[length(fit_rows) + 1]] <- data.frame(
      outcome = on, model = ifelse(nm == "A", "Adj-A", "Adj-B (+ hour)"),
      parameters = length(coef(mod)), log_likelihood = ll, AIC = -2 * ll + 2 * length(coef(mod)),
      McFadden_R2 = 1 - ll / ll0,
      AUC = auc_grouped(fitted(mod), cells[[o$yes]], cells$n - cells[[o$yes]]))
    or_rows[[length(or_rows) + 1]] <- or_table(mod, on, ifelse(nm == "A", "Adj-A", "Adj-B (+ hour)"))
    g <- gvif_table(mod); g$outcome <- on; g$model <- ifelse(nm == "A", "Adj-A", "Adj-B (+ hour)")
    gvif_rows[[length(gvif_rows) + 1]] <- g
  }
  # Likelihood-ratio tests (the difference in deviance between nested models)
  lr <- function(small, big, label) {
    a <- anova(small, big, test = "Chisq")
    data.frame(outcome = on, block = label, chi_square = a$Deviance[2], df = a$Df[2], p = a$`Pr(>Chi)`[2])
  }
  lr_rows[[length(lr_rows) + 1]] <- lr(update(mA, . ~ . - light - weather), mA, "Lighting + weather (from Adj-A)")
  lr_rows[[length(lr_rows) + 1]] <- lr(update(mA, . ~ . - surface),         mA, "Road surface (from Adj-A)")
  lr_rows[[length(lr_rows) + 1]] <- lr(mA, mB,                                  "Hour of day (Adj-B vs Adj-A)")
}
fit_tab  <- do.call(rbind, fit_rows)
or_tab   <- do.call(rbind, or_rows)
lr_tab   <- do.call(rbind, lr_rows)
gvif_tab <- do.call(rbind, gvif_rows)
or_tab$or_ci <- sprintf("%.3f (%.3f-%.3f)", or_tab$odds_ratio, or_tab$lower, or_tab$upper)

show(fit_tab, "Table 9.2  Fit of the severity models")
show(lr_tab,  "Likelihood-ratio tests for each block")
g_ksi <- gvif_tab[gvif_tab$outcome == "Fatal or serious", ]
show(g_ksi[g_ksi$model == "Adj-B (+ hour)", c("term", "df", "gvif", "gvif_adjusted")], "GVIF (fatal or serious, Adj-B)")
show(g_ksi[g_ksi$term == "light", c("model", "gvif", "gvif_adjusted")], "GVIF for lighting without and with hour of day")
show(or_tab[or_tab$model == "Adj-B (+ hour)" & or_tab$variable %in% c("light", "weather", "surface"),
            c("outcome", "variable", "level", "or_ci")], "Table 9.3  lighting, weather and road surface (Adj-B)")
write_tab(fit_tab, "severity_model_fit"); write_tab(or_tab, "severity_odds_ratios")
write_tab(lr_tab, "severity_likelihood_ratio_tests"); write_tab(gvif_tab, "severity_gvif")


# =========== PART D. FREQUENCY MODEL (NEGATIVE BINOMIAL, DAILY COUNTS) =======
note("")
note("PART D. Frequency model")
# One row per calendar day in the years supplied (2008 is not in the file and is not zero-filled).
years    <- sort(unique(d$Year))
calendar <- do.call(c, lapply(years, function(y) seq(as.Date(paste0(y, "-01-01")), as.Date(paste0(y, "-12-31")), by = "day")))
lv       <- as.character(calendar)
day      <- factor(as.character(d$date), levels = lv)
daily <- data.frame(date = calendar)
daily$count <- as.integer(table(day))
# Daily condition measures: the SHARE of that day's accidents in each condition.
daily$dark_n    <- as.integer(table(factor(as.character(d$date[d$is_dark]), levels = lv)))
known_weather   <- !is.na(d$weather)
bad_weather     <- known_weather & d$weather != "Fine"
daily$known_w_n <- as.integer(table(factor(as.character(d$date[known_weather]), levels = lv)))
daily$bad_w_n   <- as.integer(table(factor(as.character(d$date[bad_weather]),   levels = lv)))
daily$nondry_n  <- as.integer(table(factor(as.character(d$date[d$surface != "Dry"]), levels = lv)))
stopifnot(sum(daily$count) == nrow(d), min(daily$count) > 0)
daily$dark_share   <- daily$dark_n / daily$count
daily$weather_share <- daily$bad_w_n / daily$known_w_n      # share of rain, snow or fog among accidents with known weather
daily$nondry_share <- daily$nondry_n / daily$count
daily$dark_z    <- as.numeric(scale(daily$dark_share))      # per one standard deviation
daily$weather_z <- as.numeric(scale(daily$weather_share))
daily$nondry_z  <- as.numeric(scale(daily$nondry_share))
daily$weekday <- factor(format(daily$date, "%A"),
                        levels = c("Friday", "Monday", "Tuesday", "Wednesday", "Thursday", "Saturday", "Sunday"))
daily$month <- factor(as.integer(format(daily$date, "%m")), levels = 1:12, labels = month.name)
daily$year  <- factor(format(daily$date, "%Y"))
note("Days: ", nrow(daily), "; mean daily count: ", round(mean(daily$count), 1))

pois <- glm(count ~ weekday + month + year, family = poisson(), data = daily)
disp <- sum(residuals(pois, type = "pearson")^2) / df.residual(pois)
nb0  <- MASS::glm.nb(count ~ 1, data = daily)
nb1  <- MASS::glm.nb(count ~ weekday + month + year, data = daily)
nb2  <- MASS::glm.nb(count ~ weekday + month + year + weather_z + nondry_z + dark_z, data = daily)
note("Poisson dispersion with calendar terms (Pearson chi-square / df): ", round(disp, 2))

nb_compare <- data.frame(
  model = c("Negative binomial, intercept only", "Negative binomial, calendar terms",
            "Negative binomial, calendar terms + conditions"),
  log_likelihood = c(as.numeric(logLik(nb0)), as.numeric(logLik(nb1)), as.numeric(logLik(nb2))),
  AIC = c(AIC(nb0), AIC(nb1), AIC(nb2)), BIC = c(BIC(nb0), BIC(nb1), BIC(nb2)),
  squared_correlation_observed_fitted = c(NA, cor(daily$count, fitted(nb1))^2,
                                          cor(daily$count, fitted(nb2))^2))
show(nb_compare, "Table 9.5  Comparison of daily-count models")
lr_daily <- 2 * (as.numeric(logLik(nb2)) - as.numeric(logLik(nb1)))
note("Likelihood-ratio test, conditions added to calendar terms: chi-square = ", round(lr_daily, 1),
     ", df = 3, p = ", format.pval(pchisq(lr_daily, 3, lower.tail = FALSE), digits = 3))

cm <- coef(summary(nb2))
irr_tab <- data.frame(term = rownames(cm), IRR = exp(cm[, 1]), lower = exp(cm[, 1] - z95 * cm[, 2]),
                      upper = exp(cm[, 1] + z95 * cm[, 2]), p = cm[, 4], row.names = NULL)
irr_tab$irr_ci <- sprintf("%.3f (%.3f-%.3f)", irr_tab$IRR, irr_tab$lower, irr_tab$upper)
show(irr_tab[, c("term", "irr_ci", "p")], "Table 9.6  Incidence rate ratios, full daily model")

# Collinearity among the three daily condition measures (r and VIF)
shares <- daily[, c("weather_share", "nondry_share", "dark_share")]
note("Correlation, rain/snow/fog share vs non-dry surface share: r = ", round(cor(shares$weather_share, shares$nondry_share), 2))
vif_daily <- sapply(names(shares), function(v) 1 / (1 - summary(lm(reformulate(setdiff(names(shares), v), v), data = shares))$r.squared))
show(round(vif_daily, 2), "VIF of the three daily condition measures")
write_tab(nb_compare, "daily_model_comparison"); write_tab(irr_tab, "daily_model_irr")
write_tab(daily[, c("date", "count", "dark_share", "weather_share", "nondry_share")], "daily_counts_and_shares")
write_tab(data.frame(poisson_dispersion = disp, nb_theta = nb2$theta, nb_converged = nb2$converged,
                     days = nrow(daily)), "daily_model_diagnostics")


# ======================== PART E. FIGURES AND LOG ============================
navy <- "#214B6B"; teal <- "#2C807B"
start_png <- function(name, w = 1800, h = 1000) {
  png(file.path(out, "figures", paste0(name, ".png")), width = w, height = h, res = 200)
  par(mar = c(4.5, 4.8, 2.5, 1.2), las = 1, cex = 0.9, bty = "l")
}
# Figure 9.1: forest plot of Adj-B odds ratios for lighting, weather and road surface
fp <- or_tab[or_tab$model == "Adj-B (+ hour)" & or_tab$variable %in% c("light", "weather", "surface"), ]
lev <- unique(fp$level)
start_png("fig9_1_forest", 1900, 1100)
par(mar = c(4.5, 11, 2.5, 1.2))
xr <- range(c(fp$lower, fp$upper))
plot(NA, xlim = xr, ylim = c(0.5, length(lev) + 0.5), log = "x", yaxt = "n", ylab = "",
     xlab = "Adjusted odds ratio (log scale, Adj-B)", main = "Lighting, weather and road surface")
abline(v = 1, lty = 2, col = "grey50")
axis(2, at = seq_along(lev), labels = rev(lev), las = 1, cex.axis = 0.8)
for (k in seq_along(outcomes)) {
  s <- fp[fp$outcome == names(outcomes)[k], ]; y <- length(lev) + 1 - match(s$level, lev) + c(0.15, -0.15)[k]
  col <- c(navy, teal)[k]
  segments(s$lower, y, s$upper, y, col = col, lwd = 2); points(s$odds_ratio, y, pch = 16, col = col)
}
legend("bottomright", legend = names(outcomes), col = c(navy, teal), pch = 16, lwd = 2, bty = "n")
dev.off()
# Figure 9.2: observed against fitted daily counts, and relative residuals over time
start_png("fig9_2_nb", 2000, 950)
par(mfrow = c(1, 2), mar = c(4.5, 5, 2.5, 1.2))
plot(fitted(nb2), daily$count, pch = 16, cex = 0.4, col = adjustcolor(navy, 0.4),
     xlab = "Fitted daily count", ylab = "Observed daily count", main = "(a) Observed against fitted")
abline(0, 1, col = teal, lwd = 2)
plot(daily$date, (daily$count - fitted(nb2)) / fitted(nb2), pch = 16, cex = 0.4, col = adjustcolor(navy, 0.4),
     xlab = "Date", ylab = "Relative residual (observed - fitted) / fitted", main = "(b) Relative residuals")
abline(h = 0, col = teal, lwd = 2)
dev.off()

writeLines(LOG, file.path(out, "analysis_log.txt"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Complete. Output folder:", normalizePath(out), "\n")
