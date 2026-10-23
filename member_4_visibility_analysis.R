###############################################################################
#  MEMBER 4 - VISIBILITY ANALYSIS
#  Project : Statistical Modeling of Weather Severity and Highway Collision Rates
#  Dataset : UK road accidents 2005-2014 (STATS19-style file, 33 columns)
#
#  What the project brief asks of Member 4 (Section 6 + 14):
#     - Visibility distributions
#     - Relationship plots (visibility vs crash FREQUENCY and vs SEVERITY)
#     - Statistical tests / models
#
#  IMPORTANT - READ THIS FIRST
#  This dataset has NO measured visibility column (no metres / km).
#  The only visibility-related information is:
#     * Weather_Conditions : "Fog or mist" is the only category that directly
#                            describes restricted visibility
#     * Light_Conditions   : Daylight / Darkness (lit, unlit, no lighting)
#  So visibility is analysed through TRANSPARENT PROXIES built from those two
#  columns (see section 2). If Member 1's merged weather data gives you a real
#  numeric visibility column, set VIS_NUM_COL below and section 9 will run too.
#
#  HOW TO RUN
#     1. Open this file in RStudio, click "Source" (or Ctrl+Shift+S)
#     2. A file-picker opens -> choose your big CSV (UK_Accident.csv)
#     3. Everything is saved to a folder called  Member4_output  (created in
#        your working directory): PNG plots, CSV tables, and a text log of all
#        console output that you can copy into your report.
#  Runtime on 1.5 million rows: roughly 2-6 minutes. RAM needed: ~2-3 GB.
###############################################################################


# ============================== 0. SETTINGS ==================================
DATA_FILE   <- NULL            # NULL = file picker.  Or e.g. "C:/Users/me/UK_Accident.csv"
OUT_DIR     <- "Member4_output"
DEDUP       <- TRUE            # drop exact duplicate accident rows (see section 1)
VIS_NUM_COL <- NULL            # e.g. "Visibility_km" if a real numeric column exists
FOG_DAY_Q   <- 0.90            # "high-fog day" = top 10% of days by # fog-coded accidents
SEED        <- 2026

options(scipen = 999, width = 110, digits = 6)
set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE)

for (p in c("ggplot2", "data.table", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(scales) })
# MASS (for glm.nb) ships with R itself, no install needed.

# ---------------------------- small helpers ----------------------------------
LOG <- character(0)
section <- function(title) {
  s <- c("", strrep("=", 78), title, strrep("=", 78))
  cat(s, sep = "\n"); LOG <<- c(LOG, s)
}
note <- function(...) { s <- paste0(...); cat(s, sep = "\n"); LOG <<- c(LOG, s) }
show <- function(x, label = NULL) {
  out <- capture.output(print(x))
  if (!is.null(label)) out <- c(label, out)
  cat(out, sep = "\n"); LOG <<- c(LOG, out)
}
save_tab <- function(x, name) write.csv(x, file.path(OUT_DIR, paste0(name, ".csv")), row.names = FALSE)
save_plot <- function(p, name, w = 9, h = 6) {
  print(p)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 150, bg = "white")
}
wilson <- function(x, n, conf = 0.95) {            # Wilson score interval for a proportion
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  data.frame(p = p, lo = ctr - half, hi = ctr + half)
}
chisq_report <- function(x, y, label) {             # chi-square + Cramer's V + assumption check
  tab <- table(x, y)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  ct  <- suppressWarnings(chisq.test(tab))
  V   <- sqrt(as.numeric(ct$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
  data.frame(Test = label, N = sum(tab), Chi_sq = round(as.numeric(ct$statistic), 1),
             df = as.integer(ct$parameter), p_value = format.pval(ct$p.value, digits = 3),
             Cramers_V = round(V, 4), Min_expected_count = round(min(ct$expected), 1))
}
cohens_d <- function(a, b) {
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  (mean(a) - mean(b)) / sp
}
or_table <- function(model, model_label, outcome_label) {   # odds ratios for the visibility terms only
  cf <- summary(model)$coefficients
  cf <- cf[grepl("^(Light_M|Weather_M)", rownames(cf)), , drop = FALSE]
  data.frame(Outcome = outcome_label, Model = model_label,
             Term = sub("^(Light_M|Weather_M)", "", rownames(cf)),
             Block = ifelse(grepl("^Light_M", rownames(cf)), "Lighting (ref: Daylight)", "Weather (ref: Fine)"),
             OR = exp(cf[, 1]), CI_low = exp(cf[, 1] - 1.96 * cf[, 2]),
             CI_high = exp(cf[, 1] + 1.96 * cf[, 2]), p_value = cf[, 4], row.names = NULL)
}
SEV_COLS <- c(Fatal = "#b2182b", Serious = "#ef8a62", Slight = "#67a9cf")
theme_set(theme_minimal(base_size = 12))


# ============================ 1. LOAD + AUDIT ================================
section("1. LOAD DATA + QUALITY AUDIT")
if (is.null(DATA_FILE)) DATA_FILE <- file.choose()

keep <- c("Accident_Index", "Location_Easting_OSGR", "Location_Northing_OSGR", "Police_Force",
          "Accident_Severity", "Number_of_Vehicles", "Number_of_Casualties", "Date", "Day_of_Week",
          "Time", "Speed_limit", "Road_Type", "Light_Conditions", "Weather_Conditions",
          "Road_Surface_Conditions", "Urban_or_Rural_Area", "Year")
hdr  <- names(fread(DATA_FILE, nrows = 0))
miss <- setdiff(keep, hdr)
if (length(miss) > 0) stop("These columns were not found in the file header: ",
                           paste(miss, collapse = ", "), "\nHeader is: ", paste(hdr, collapse = ", "))

vis_found <- grep("visib", hdr, ignore.case = TRUE, value = TRUE)
if (length(vis_found) > 0) {
  note("Columns with 'visib' in the name were found: ", paste(vis_found, collapse = ", "))
  if (is.null(VIS_NUM_COL)) { VIS_NUM_COL <- vis_found[1]
    note("-> using '", VIS_NUM_COL, "' as the numeric visibility column (set VIS_NUM_COL to override).") }
} else note("No column with 'visib' in its name: visibility must be analysed through proxies (Light + Weather).")

d <- fread(DATA_FILE, select = unique(c(keep, VIS_NUM_COL)), data.table = FALSE)
d$Date <- as.character(d$Date); d$Time <- as.character(d$Time)
note("Rows read: ", format(nrow(d), big.mark = ","))

# ---- duplicate check ---------------------------------------------------------
# Accident_Index alone can NOT be trusted here (index format differs between the
# source year-files), so an exact duplicate = same values in ALL these fields.
key <- paste(d$Year, d$Accident_Index, d$Date, d$Time, d$Location_Easting_OSGR,
             d$Location_Northing_OSGR, d$Accident_Severity, d$Number_of_Vehicles,
             d$Number_of_Casualties, sep = "|")
is_dup <- duplicated(key); rm(key)
note("Unique Accident_Index values      : ", format(length(unique(d$Accident_Index)), big.mark = ","))
note("Exact duplicate rows (all fields) : ", format(sum(is_dup), big.mark = ","))
yr_tab <- data.frame(Year = sort(unique(d$Year)),
                     Rows = as.integer(table(d$Year)),
                     Exact_duplicates = as.integer(tapply(is_dup, d$Year, sum)))
show(yr_tab, "Rows and exact duplicates by year:")
save_tab(yr_tab, "T0_rows_and_duplicates_by_year")
note("NOTE: DfT's published total for 2013 is 138,660 accidents; compare your 2013 row count above.")
if (DEDUP && any(is_dup)) { d <- d[!is_dup, ]; note("DEDUP = TRUE -> removed duplicates. Rows now: ", format(nrow(d), big.mark = ",")) }
rm(is_dup)

# ---- date / time -------------------------------------------------------------
best_fmt <- local({                     # pick the date format that parses AND agrees with the Year column
  idx  <- sample(nrow(d), min(50000, nrow(d)))
  fmts <- c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y", "%m/%d/%Y", "%d/%m/%y")
  sc   <- sapply(fmts, function(f) { p <- as.Date(d$Date[idx], format = f)
                                     mean(!is.na(p) & as.integer(format(p, "%Y")) == d$Year[idx], na.rm = TRUE) })
  fmts[which.max(sc)]
})
note("Date format detected: ", best_fmt)
d$Date2 <- as.Date(d$Date, format = best_fmt)
note("Rows with unparseable Date: ", sum(is.na(d$Date2)))
d$Hour  <- suppressWarnings(as.integer(sub(":.*", "", d$Time)))

# ---- severity ----------------------------------------------------------------
sev <- d$Accident_Severity
if (!is.numeric(sev)) sev <- match(tolower(trimws(sev)), c("fatal", "serious", "slight"))
d$Severity <- factor(sev, levels = 1:3, labels = c("Fatal", "Serious", "Slight"))   # DfT coding: 1 Fatal, 2 Serious, 3 Slight
d$KSI      <- as.integer(d$Severity %in% c("Fatal", "Serious"))                    # Killed or Seriously Injured
d$Fatal    <- as.integer(d$Severity == "Fatal")
note("Severity coding used: 1 = Fatal, 2 = Serious, 3 = Slight (DfT STATS19 convention).")
show(round(100 * prop.table(table(d$Severity, useNA = "ifany")), 2), "Severity % (should be roughly 1.3 / 13.6 / 85.1):")


# ===================== 2. BUILD THE VISIBILITY PROXIES =======================
section("2. VISIBILITY PROXIES (built from Light_Conditions and Weather_Conditions)")

weather_codes <- c("1" = "fine without high winds", "2" = "raining without high winds", "3" = "snowing without high winds",
                   "4" = "fine with high winds", "5" = "raining with high winds", "6" = "snowing with high winds",
                   "7" = "fog or mist", "8" = "other", "9" = "unknown")
light_codes   <- c("1" = "daylight", "4" = "darkness - lights lit", "5" = "darkness - lights unlit",
                   "6" = "darkness - no lighting", "7" = "darkness - lighting unknown")

w <- d$Weather_Conditions; l <- d$Light_Conditions
w <- if (is.numeric(w)) unname(weather_codes[as.character(w)]) else tolower(trimws(as.character(w)))
l <- if (is.numeric(l)) unname(light_codes[as.character(l)])   else tolower(trimws(as.character(l)))
w[is.na(w)] <- "unknown"; l[is.na(l)] <- "unknown"

wg <- rep("Other/Unknown", length(w))
wg[grepl("^fine", w)]    <- "Fine"
wg[grepl("^rain", w)]    <- "Rain"
wg[grepl("^snow", w)]    <- "Snow"
wg[grepl("fog|mist", w)] <- "Fog/Mist"
d$Weather_Group <- factor(wg, levels = c("Fine", "Rain", "Snow", "Fog/Mist", "Other/Unknown"))

lg <- rep("Unknown", length(l)); dark <- grepl("^dark", l)
lg[grepl("^daylight", l)]                                   <- "Daylight"
lg[dark & grepl("unknown", l)]                              <- "Dark: lighting unknown"
lg[dark & grepl("unlit|no street|no light", l)]             <- "Dark: unlit / no lighting"
lg[dark & grepl("present and lit|lights lit", l)]           <- "Dark: lit"
lg[dark & lg == "Unknown"]                                  <- "Dark: lighting unknown"
d$Light_Group <- factor(lg, levels = c("Daylight", "Dark: lit", "Dark: unlit / no lighting",
                                       "Dark: lighting unknown", "Unknown"))

show(table(Raw_Weather = d$Weather_Conditions, Group = d$Weather_Group), "Check - raw Weather text -> Weather_Group:")
show(table(Raw_Light = d$Light_Conditions, Group = d$Light_Group),       "Check - raw Light text -> Light_Group:")
note("^ Read these two tables: every raw label must land in the group you expect. If not, tell me and I'll fix the mapping.")

# versions used for INFERENCE (unknown / other categories removed)
d$Light_M   <- factor(ifelse(d$Light_Group == "Unknown", NA, as.character(d$Light_Group)),
                      levels = c("Daylight", "Dark: lit", "Dark: unlit / no lighting", "Dark: lighting unknown"))
d$Weather_M <- factor(ifelse(d$Weather_Group == "Other/Unknown", NA, as.character(d$Weather_Group)),
                      levels = c("Fine", "Rain", "Snow", "Fog/Mist"))

# combined ordinal proxy: 0 (best) ... 4 (worst)   -- an ASSUMPTION-BASED index, not a measurement
light_pts   <- c("Daylight" = 0, "Dark: lit" = 1, "Dark: unlit / no lighting" = 2)     # 'lighting unknown' -> NA
weather_pts <- c("Fine" = 0, "Rain" = 1, "Snow" = 1, "Fog/Mist" = 2)
d$Vis_Score <- unname(light_pts[as.character(d$Light_M)] + weather_pts[as.character(d$Weather_M)])
d$Vis_Class <- cut(d$Vis_Score, c(-Inf, 0, 1, 2, Inf), labels = c("0 (best)", "1", "2", "3-4 (worst)"))
note("Vis_Score = lighting points (Daylight 0, Dark-lit 1, Dark-unlit/none 2) + weather points (Fine 0, Rain 1, Snow 1, Fog/Mist 2).")
note("Because the score is a modelling choice, EVERY key result is also reported for Light and Weather separately.")
note("Rows usable for the combined score: ", format(sum(!is.na(d$Vis_Score)), big.mark = ","),
     "  (", round(100 * mean(!is.na(d$Vis_Score)), 1), "% of rows)")

# other model covariates
d$Speed_limit[!is.na(d$Speed_limit) & d$Speed_limit <= 0] <- NA      # 0 / -1 = data missing
d$Speed_f   <- cut(d$Speed_limit, c(-Inf, 20, 30, 40, 50, 60, Inf), labels = c("<=20", "30", "40", "50", "60", "70"))
ar          <- as.character(d$Urban_or_Rural_Area)
d$Area      <- factor(ifelse(ar %in% c("1", "Urban"), "Urban", ifelse(ar %in% c("2", "Rural"), "Rural", NA)), levels = c("Urban", "Rural"))
rt          <- as.character(d$Road_Type); rt[is.na(rt) | rt == ""] <- "Unknown"
d$Road_Type_f <- factor(rt); if ("Single carriageway" %in% levels(d$Road_Type_f)) d$Road_Type_f <- relevel(d$Road_Type_f, "Single carriageway")
d$Veh_cat   <- cut(d$Number_of_Vehicles, c(0, 1, 2, 3, Inf), labels = c("1", "2", "3", "4+"))
d$Hour_band <- cut(d$Hour, c(-1, 5, 9, 15, 19, 23), labels = c("Night 00-05", "Morning 06-09", "Midday 10-15", "Evening 16-19", "Late 20-23"))
d$Weekend   <- factor(ifelse(as.integer(format(d$Date2, "%u")) >= 6, "Weekend", "Weekday"), levels = c("Weekday", "Weekend"))
d$Year_f    <- factor(d$Year)


# ==================== 3. VISIBILITY DISTRIBUTIONS (Objective) =================
section("3. VISIBILITY DISTRIBUTIONS")
dist_tab <- function(g, nm) { t <- table(g, useNA = "ifany"); data.frame(Variable = nm, Category = names(t),
                              Accidents = as.integer(t), Percent = round(100 * as.integer(t) / sum(t), 2)) }
T1 <- rbind(dist_tab(d$Light_Group, "Light_Group"), dist_tab(d$Weather_Group, "Weather_Group"),
            dist_tab(d$Vis_Score, "Vis_Score (0 best - 4 worst)"))
show(T1, "Table 1. Distribution of the visibility proxies (all accident records)"); save_tab(T1, "T1_visibility_distributions")

bar_dist <- function(g, title, xlab, fill = "#4d4d4d") {
  df <- as.data.frame(table(g)); df$pct <- 100 * df$Freq / sum(df$Freq)
  ggplot(df, aes(g, Freq)) + geom_col(fill = fill) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3.5) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    labs(title = title, x = xlab, y = "Number of accident records")
}
save_plot(bar_dist(d$Light_Group,   "Accident records by lighting condition",   "Light_Group") +
            theme(axis.text.x = element_text(angle = 20, hjust = 1)), "P1_distribution_lighting")
save_plot(bar_dist(d$Weather_Group, "Accident records by weather group",        "Weather_Group"), "P2_distribution_weather")
save_plot(bar_dist(factor(d$Vis_Score), "Accident records by combined visibility score (0 = best, 4 = worst)", "Vis_Score", "#2b6a99"),
          "P3_distribution_vis_score")

# heat-tile of Light x Weather
tile <- as.data.frame(table(Light = d$Light_Group, Weather = d$Weather_Group))
save_plot(ggplot(tile, aes(Weather, Light, fill = log10(Freq + 1))) + geom_tile(colour = "white") +
            geom_text(aes(label = comma(Freq)), size = 3.3) + scale_fill_gradient(low = "#f7fbff", high = "#2171b5", guide = "none") +
            labs(title = "Lighting x weather (record counts)", x = NULL, y = NULL), "P4_light_by_weather_tiles")

# exposure warning plot: when do dark accidents happen?
hr <- as.data.frame(table(Hour = d$Hour, Light = d$Light_Group)); hr <- hr[hr$Light != "Unknown", ]
save_plot(ggplot(hr, aes(as.integer(as.character(Hour)), Freq, colour = Light)) + geom_line(linewidth = 1) +
            scale_y_continuous(labels = comma) + labs(title = "Accident records by hour of day and lighting condition",
            subtitle = "Counts follow traffic (rush hours) as much as visibility - raw counts are NOT risk", x = "Hour of day", y = "Records", colour = NULL),
          "P5_hour_by_lighting", 10, 6)

# fog seasonality
fog_m <- as.data.frame(table(Month = factor(format(d$Date2[d$Weather_Group == "Fog/Mist"], "%b"), levels = month.abb)))
save_plot(ggplot(fog_m, aes(Month, Freq)) + geom_col(fill = "#7f7f7f") + scale_y_continuous(labels = comma) +
            labs(title = "Fog/mist-coded accident records by calendar month (all years pooled)", y = "Records"), "P6_fog_by_month")


# =============== 4. SEVERITY vs VISIBILITY - TABLES AND TESTS =================
section("4. SEVERITY vs VISIBILITY: crosstabs, chi-square, trend test")
note("H0: severity is NOT associated with visibility (proxy).   H1: severity IS associated with visibility.")
note("Rows with 'Unknown/Other' weather or 'Unknown' light are excluded from inference; lighting-unknown-but-dark is kept as its own level.")

sev_table <- function(g, label) {
  ok <- !is.na(g) & !is.na(d$Severity)
  tab <- table(Group = g[ok], Severity = d$Severity[ok])
  df  <- as.data.frame.matrix(tab); df$N <- rowSums(tab)
  df$Pct_Fatal   <- round(100 * df$Fatal / df$N, 2)
  df$Pct_Serious <- round(100 * df$Serious / df$N, 2)
  ci <- wilson(df$Fatal + df$Serious, df$N)
  df$Pct_KSI <- round(100 * ci$p, 2); df$KSI_CI_low <- round(100 * ci$lo, 2); df$KSI_CI_high <- round(100 * ci$hi, 2)
  cbind(Variable = label, Group = rownames(df), df, row.names = NULL)
}
T2 <- rbind(sev_table(d$Light_M, "Lighting"), sev_table(d$Weather_M, "Weather"), sev_table(d$Vis_Class, "Vis_Class"))
show(T2, "Table 2. Severity by visibility proxy (KSI = Fatal + Serious, with 95% Wilson CI)"); save_tab(T2, "T2_severity_by_visibility")

T3 <- rbind(chisq_report(d$Light_M,   d$Severity, "Severity x Lighting"),
            chisq_report(d$Weather_M, d$Severity, "Severity x Weather"),
            chisq_report(d$Vis_Class, d$Severity, "Severity x Combined visibility class"))
show(T3, "Table 3. Chi-square tests of independence (with effect size)"); save_tab(T3, "T3_chisquare_tests")
note("Reading Table 3: with >1 million rows the p-value will be ~0 for ANY tiny difference. Judge importance by Cramer's V:")
note("   Rough benchmarks for these tables (Cohen, df* = 2): V ~0.07 small, ~0.21 medium, ~0.35 large. Min_expected_count must be >= 5.")
show(round(chisq.test(table(d$Light_M, d$Severity))$stdres, 1),   "Standardised residuals, Lighting x Severity (|value| > 3 = cell drives the association):")
show(round(chisq.test(table(d$Weather_M, d$Severity))$stdres, 1), "Standardised residuals, Weather x Severity:")

# 2x2 contrasts on the raw rates
ct2 <- function(g1, g2, gv, lab) {
  ok <- !is.na(gv) & gv %in% c(g1, g2)
  a <- sum(d$KSI[ok & gv == g1]); n1 <- sum(ok & gv == g1)
  b <- sum(d$KSI[ok & gv == g2]); n2 <- sum(ok & gv == g2)
  pt <- prop.test(c(a, b), c(n1, n2), correct = FALSE)
  or <- (a / (n1 - a)) / (b / (n2 - b)); se <- sqrt(1 / a + 1 / (n1 - a) + 1 / b + 1 / (n2 - b))
  data.frame(Contrast = lab, KSI_pct_1 = round(100 * a / n1, 2), KSI_pct_2 = round(100 * b / n2, 2),
             Crude_OR = round(or, 3), OR_CI_low = round(exp(log(or) - 1.96 * se), 3),
             OR_CI_high = round(exp(log(or) + 1.96 * se), 3), p_value = format.pval(pt$p.value, digits = 3))
}
T4 <- rbind(ct2("Fog/Mist", "Fine", d$Weather_M, "Fog/Mist vs Fine (KSI)"),
            ct2("Dark: unlit / no lighting", "Daylight", d$Light_M, "Dark unlit/none vs Daylight (KSI)"),
            ct2("Dark: lit", "Daylight", d$Light_M, "Dark lit vs Daylight (KSI)"))
show(T4, "Table 4. Key 2x2 contrasts (crude; two-proportion z-test)"); save_tab(T4, "T4_key_contrasts_crude")

# ordinal trend across the combined score
sc <- as.data.frame(table(Score = d$Vis_Score[!is.na(d$Vis_Score)]))
sc$KSI <- as.integer(tapply(d$KSI[!is.na(d$Vis_Score)], d$Vis_Score[!is.na(d$Vis_Score)], sum))
sc$Fatal <- as.integer(tapply(d$Fatal[!is.na(d$Vis_Score)], d$Vis_Score[!is.na(d$Vis_Score)], sum))
sc$Pct_KSI <- 100 * sc$KSI / sc$Freq; ci <- wilson(sc$KSI, sc$Freq); sc$lo <- 100 * ci$lo; sc$hi <- 100 * ci$hi
sc$Pct_Fatal <- 100 * sc$Fatal / sc$Freq
t5 <- data.frame(Score = as.character(sc$Score), Records = sc$Freq, KSI = sc$KSI, Pct_KSI = round(sc$Pct_KSI, 2), Pct_Fatal = round(sc$Pct_Fatal, 2))
show(t5, "Table 5. KSI and Fatal rate by combined visibility score (0 = best ... 4 = worst)")
trend_ksi   <- prop.trend.test(sc$KSI,   sc$Freq, score = as.numeric(as.character(sc$Score)))
trend_fatal <- prop.trend.test(sc$Fatal, sc$Freq, score = as.numeric(as.character(sc$Score)))
note("Cochran-Armitage trend test, KSI rate ~ visibility score  : X-squared = ", round(trend_ksi$statistic, 1),
     ", p = ", format.pval(trend_ksi$p.value, digits = 3))
note("Cochran-Armitage trend test, Fatal rate ~ visibility score: X-squared = ", round(trend_fatal$statistic, 1),
     ", p = ", format.pval(trend_fatal$p.value, digits = 3))
save_tab(sc, "T5_trend_by_vis_score")

# ---- plots: severity vs visibility
stack_plot <- function(T, var, title) {
  x <- T[T$Variable == var, ]
  long <- rbind(data.frame(Group = x$Group, Severity = "Serious", Pct = x$Pct_Serious),
                data.frame(Group = x$Group, Severity = "Fatal",   Pct = x$Pct_Fatal))
  long$Group <- factor(long$Group, levels = x$Group); long$Severity <- factor(long$Severity, levels = c("Serious", "Fatal"))
  ggplot(long, aes(Group, Pct, fill = Severity)) + geom_col(width = 0.7) +
    geom_text(data = x, aes(Group, Pct_KSI, label = sprintf("KSI %.1f%%", Pct_KSI), fill = NULL), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = SEV_COLS) + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = title, subtitle = "Share of accidents that were Fatal or Serious (KSI); the remainder were Slight",
         x = NULL, y = "% of accidents in the group") + theme(axis.text.x = element_text(angle = 15, hjust = 1))
}
save_plot(stack_plot(T2, "Lighting",  "Accident severity by lighting condition"), "P7_severity_by_lighting")
save_plot(stack_plot(T2, "Weather",   "Accident severity by weather group"),      "P8_severity_by_weather")
save_plot(stack_plot(T2, "Vis_Class", "Accident severity by combined visibility class (0 best - 3-4 worst)"), "P9_severity_by_vis_class")

save_plot(ggplot(sc, aes(as.numeric(as.character(Score)), Pct_KSI)) +
            geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.08) + geom_line(colour = "grey40") + geom_point(size = 3, colour = "#b2182b") +
            labs(title = "RELATIONSHIP PLOT: KSI rate vs combined visibility score",
                 subtitle = "Points = % Fatal/Serious with 95% CI; higher score = worse visibility proxy",
                 x = "Visibility score (0 = daylight + fine ... 4 = dark unlit + fog)", y = "% of accidents Fatal or Serious"),
          "P10_relationship_KSI_vs_vis_score")


# ================= 5. CASUALTIES (second severity measure) ===================
section("5. NUMBER OF CASUALTIES vs VISIBILITY (non-parametric + t-test)")
note("Number_of_Casualties is a highly right-skewed count (max 93), so Kruskal-Wallis / Wilcoxon are the primary tests;")
note("Welch t-tests are shown because the brief asks for them and n is huge (CLT), together with Cohen's d.")
cas_tab <- function(g, label) {
  ok <- !is.na(g); parts <- split(d$Number_of_Casualties[ok], droplevels(g[ok]))
  df <- do.call(rbind, lapply(parts, function(z) data.frame(N = length(z), Mean = mean(z), Median = median(z),
                                                            SD = sd(z), Pct_2plus = 100 * mean(z >= 2))))
  df$SE <- df$SD / sqrt(df$N); cbind(Variable = label, Group = rownames(df), round(df, 3), row.names = NULL)
}
T6 <- rbind(cas_tab(d$Light_M, "Lighting"), cas_tab(d$Weather_M, "Weather"), cas_tab(d$Vis_Class, "Vis_Class"))
show(T6, "Table 6. Casualties per accident by visibility proxy"); save_tab(T6, "T6_casualties_by_visibility")

kw <- function(g, label) {
  ok <- !is.na(g); k <- kruskal.test(d$Number_of_Casualties[ok] ~ droplevels(g[ok]))
  data.frame(Test = label, N = sum(ok), H = round(as.numeric(k$statistic), 1), df = as.integer(k$parameter),
             p_value = format.pval(k$p.value, digits = 3), Epsilon_sq = round(as.numeric(k$statistic) / (sum(ok) - 1), 5))
}
T7 <- rbind(kw(d$Light_M, "Casualties ~ Lighting"), kw(d$Weather_M, "Casualties ~ Weather"), kw(d$Vis_Class, "Casualties ~ Vis_Class"))
show(T7, "Table 7. Kruskal-Wallis tests (Epsilon_sq: <0.01 negligible, 0.01-0.06 small)"); save_tab(T7, "T7_kruskal_wallis")

two_group <- function(mask, grp, g1, g2, lab) {
  z1 <- d$Number_of_Casualties[mask & grp == g1]; z2 <- d$Number_of_Casualties[mask & grp == g2]
  tt <- t.test(z1, z2); wt <- wilcox.test(z1, z2, exact = FALSE)
  data.frame(Contrast = lab, Mean_1 = round(mean(z1), 3), Mean_2 = round(mean(z2), 3),
             Diff = round(diff(c(mean(z2), mean(z1))), 3), CI_low = round(tt$conf.int[1], 3), CI_high = round(tt$conf.int[2], 3),
             Welch_p = format.pval(tt$p.value, digits = 3), Wilcoxon_p = format.pval(wt$p.value, digits = 3),
             Cohens_d = round(cohens_d(z1, z2), 3))
}
dl <- ifelse(is.na(d$Light_M), NA, ifelse(startsWith(as.character(d$Light_M), "Dark"), "Dark", "Daylight"))
m1 <- !is.na(dl); m2 <- !is.na(d$Weather_M) & d$Weather_M %in% c("Fog/Mist", "Fine")
T8 <- rbind(two_group(m1, dl, "Dark", "Daylight", "Dark vs Daylight"),
            two_group(m2, as.character(d$Weather_M), "Fog/Mist", "Fine", "Fog/Mist vs Fine"))
show(T8, "Table 8. Welch t-test + Wilcoxon on casualties per accident (Diff = group1 - group2; |d| < 0.2 = small)"); save_tab(T8, "T8_casualty_two_group_tests")

for (v in c("Lighting", "Weather", "Vis_Class")) {
  x <- T6[T6$Variable == v, ]; x$Group <- factor(x$Group, levels = x$Group)
  save_plot(ggplot(x, aes(Group, Mean)) + geom_col(fill = "#6a6a6a", width = 0.7) +
              geom_errorbar(aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE), width = 0.2) +
              labs(title = paste("Mean casualties per accident by", v), x = NULL, y = "Mean casualties (95% CI)") +
              theme(axis.text.x = element_text(angle = 15, hjust = 1)), paste0("P11_casualties_", tolower(v)))
}


# ============== 6. LOGISTIC REGRESSION: severity ~ visibility ================
section("6. LOGISTIC MODELS: is visibility associated with severity AFTER adjusting for confounders?")
note("Outcomes: (a) KSI = Fatal or Serious vs Slight;  (b) Fatal vs Serious/Slight.")
note("Models:  Crude  |  Adj-A = + speed limit, urban/rural, road type, # vehicles, weekend, year  |  Adj-B = Adj-A + hour-of-day band")
note("Odds ratios > 1 mean HIGHER odds of the severe outcome than the reference (Daylight / Fine). Fitted on aggregated cells (exactly")
note("equivalent to individual-level logistic regression, but memory-light).")

dt <- as.data.table(d[, c("KSI", "Fatal", "Light_M", "Weather_M", "Speed_f", "Area", "Road_Type_f", "Veh_cat", "Hour_band", "Weekend", "Year_f")])
dt <- dt[complete.cases(dt)]
note("Rows used in the models (complete cases): ", format(nrow(dt), big.mark = ","), " of ", format(nrow(d), big.mark = ","))
agg <- dt[, .(ksi = sum(KSI), fat = sum(Fatal), n = .N),
          by = .(Light_M, Weather_M, Speed_f, Area, Road_Type_f, Veh_cat, Hour_band, Weekend, Year_f)]
agg <- as.data.frame(agg); rm(dt)
note("Aggregated to ", format(nrow(agg), big.mark = ","), " covariate patterns.")

f_vis  <- "Light_M + Weather_M"
f_A    <- "Speed_f + Area + Road_Type_f + Veh_cat + Weekend + Year_f"
f_B    <- paste(f_A, "+ Hour_band")
fit <- function(resp, rhs) glm(as.formula(paste0(resp, " ~ ", rhs)), family = binomial, data = agg)
agg$ksi_no <- agg$n - agg$ksi; agg$fat_no <- agg$n - agg$fat
outcomes <- list(KSI   = list(resp = "cbind(ksi, ksi_no)", lab = "Fatal or Serious vs Slight"),
                 Fatal = list(resp = "cbind(fat, fat_no)", lab = "Fatal vs Serious/Slight"))
OR_ALL <- list(); LRT_ALL <- list(); MODELS <- list()
for (o in names(outcomes)) {
  r <- outcomes[[o]]$resp; L <- outcomes[[o]]$lab
  m0 <- fit(r, f_vis);                     mA <- fit(r, paste(f_vis, "+", f_A));  mB <- fit(r, paste(f_vis, "+", f_B))
  MODELS[[o]] <- list(crude = m0, adjA = mA, adjB = mB)
  OR_ALL[[length(OR_ALL) + 1]] <- or_table(m0, "Crude", L)
  OR_ALL[[length(OR_ALL) + 1]] <- or_table(mA, "Adj-A", L)
  OR_ALL[[length(OR_ALL) + 1]] <- or_table(mB, "Adj-B (+hour)", L)
  for (nm in c("adjA", "adjB")) {
    mm <- MODELS[[o]][[nm]]; m_no <- update(mm, . ~ . - Light_M - Weather_M)
    a <- anova(m_no, mm, test = "Chisq")
    LRT_ALL[[length(LRT_ALL) + 1]] <- data.frame(Outcome = L, Model = ifelse(nm == "adjA", "Adj-A", "Adj-B (+hour)"),
              Test = "LR test: do Light+Weather add anything?", Chi_sq = round(a$Deviance[2], 1), df = a$Df[2],
              p_value = format.pval(a$`Pr(>Chi)`[2], digits = 3))
  }
}
OR_TAB <- do.call(rbind, OR_ALL); LRT_TAB <- do.call(rbind, LRT_ALL)
OR_OUT <- OR_TAB; OR_OUT[c("OR", "CI_low", "CI_high")] <- round(OR_OUT[c("OR", "CI_low", "CI_high")], 3)
OR_OUT$p_value <- format.pval(OR_OUT$p_value, digits = 3)
show(OR_OUT, "Table 9. Odds ratios (95% Wald CI) for the visibility terms"); save_tab(OR_OUT, "T9_logistic_odds_ratios")
show(LRT_TAB, "Table 10. Likelihood-ratio tests for the whole visibility block (H0: no association given covariates)"); save_tab(LRT_TAB, "T10_LR_tests")
note("Interpretation guide: report OR with its CI. An OR of 1.30 for Fog/Mist = 30% higher ODDS of a Fatal/Serious outcome vs Fine weather,")
note("holding the other variables fixed. This is ASSOCIATION, not causation. Compare Crude vs Adjusted: a big change means confounding.")
note("Caution: hour-of-day and darkness are strongly collinear (night hours are almost always dark), so the Adj-B lighting ORs are")
note("identified mostly from the hours where daylight/darkness changes with the season. Discuss both Adj-A and Adj-B.")

fp <- OR_TAB; fp$Term <- factor(fp$Term, levels = rev(unique(fp$Term))); fp$Model <- factor(fp$Model, levels = c("Crude", "Adj-A", "Adj-B (+hour)"))
save_plot(ggplot(fp, aes(OR, Term, colour = Model)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
            geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.25, orientation = "y", position = position_dodge(0.6)) +
            geom_point(position = position_dodge(0.6), size = 2) + scale_x_log10() + facet_grid(Block ~ Outcome, scales = "free_y", space = "free_y") +
            labs(title = "Odds ratios for severe outcomes by visibility proxy (log scale, 95% CI)",
                 subtitle = "Right of the dashed line = higher odds than reference (Daylight / Fine)", x = "Odds ratio (log scale)", y = NULL),
          "P12_forest_plot_odds_ratios", 11, 7)

# sanity check on the main adjusted model: Hosmer-style calibration is meaningless at this N, so report fit info instead
mB <- MODELS$KSI$adjB
note("Main model (KSI, Adj-B): null deviance = ", round(mB$null.deviance, 1), ", residual deviance = ", round(mB$deviance, 1),
     ", AIC = ", round(AIC(mB), 1), ", McFadden pseudo-R2 = ", round(1 - mB$deviance / mB$null.deviance, 4))
note("(low pseudo-R2 is normal: most severity variation is not explained by any recorded covariate)")


# ============= 7. CRASH FREQUENCY vs VISIBILITY (daily level) ================
section("7. CRASH FREQUENCY vs VISIBILITY (daily counts)")
note("LIMITATION (important for your report): this file contains ONLY accidents. Days with no accidents, exposure (traffic volume),")
note("and actual weather on non-accident days are not in it, so 'visibility of the day' can only be inferred from the accidents themselves.")
note("Results in this section are therefore EXPLORATORY associations, not causal effects of visibility on crash frequency.")
note("There are no 2008 records in the dataset, so 2008 is simply absent from the daily series (no zero-filling).")

dd <- as.data.table(d[!is.na(d$Date2), c("Date2", "Weather_Group", "Light_Group", "Vis_Score", "KSI")])
dd[, Light_chr := as.character(Light_Group)]
daily <- dd[, .(n_acc = .N, n_fog = sum(Weather_Group == "Fog/Mist"), n_ksi = sum(KSI),
                mean_vis = mean(Vis_Score, na.rm = TRUE),
                share_dark = sum(startsWith(Light_chr, "Dark")) / max(1, sum(Light_chr != "Unknown"))), by = Date2]
daily <- as.data.frame(daily[order(Date2)]); rm(dd)
daily <- daily[is.finite(daily$mean_vis), ]
daily$Year_f  <- factor(format(daily$Date2, "%Y")); daily$Month <- factor(format(daily$Date2, "%m"))
daily$Weekday <- factor(format(daily$Date2, "%u")); daily$non_fog <- daily$n_acc - daily$n_fog
daily$mean_vis_z <- as.numeric(scale(daily$mean_vis))
note("Days in daily table: ", nrow(daily))
T11 <- data.frame(Variable = c("Daily accident records", "Daily fog/mist-coded records", "Daily mean visibility score"),
                  Mean = c(mean(daily$n_acc), mean(daily$n_fog), mean(daily$mean_vis)),
                  Median = c(median(daily$n_acc), median(daily$n_fog), median(daily$mean_vis)),
                  SD = c(sd(daily$n_acc), sd(daily$n_fog), sd(daily$mean_vis)),
                  Min = c(min(daily$n_acc), min(daily$n_fog), min(daily$mean_vis)),
                  Max = c(max(daily$n_acc), max(daily$n_fog), max(daily$mean_vis)))
T11[, -1] <- round(T11[, -1], 3); show(T11, "Table 11. Daily summary"); save_tab(T11, "T11_daily_summary")

# 7a. correlation + scatter
cs <- cor.test(daily$mean_vis, daily$n_acc, method = "spearman", exact = FALSE)
cp <- cor.test(daily$mean_vis, daily$n_acc, method = "pearson")
note("Daily mean visibility score vs daily accident count:  Spearman rho = ", round(cs$estimate, 3), " (p ", format.pval(cs$p.value, digits = 3),
     "),  Pearson r = ", round(cp$estimate, 3), " [95% CI ", round(cp$conf.int[1], 3), " to ", round(cp$conf.int[2], 3), "]")
save_plot(ggplot(daily, aes(mean_vis, n_acc)) + geom_point(alpha = 0.25, colour = "#3b3b3b") +
            geom_smooth(method = "loess", formula = y ~ x, colour = "#b2182b", se = TRUE) +
            labs(title = "RELATIONSHIP PLOT: daily accident count vs daily mean visibility score",
                 subtitle = paste0("Each point = one day. Spearman rho = ", round(cs$estimate, 2), ". Exploratory (accident-only data)"),
                 x = "Mean visibility score of that day's accidents (higher = poorer visibility proxy)", y = "Accident records that day"),
          "P13_relationship_daily_count_vs_vis")

# 7b. count model
pois <- glm(n_acc ~ mean_vis_z + Weekday + Month + Year_f, family = poisson, data = daily)
disp <- sum(residuals(pois, "pearson")^2) / pois$df.residual
note("Poisson dispersion (Pearson chi2 / df) = ", round(disp, 2), "  (1 = Poisson OK; >>1 = over-dispersed -> use negative binomial)")
nb <- suppressWarnings(MASS::glm.nb(n_acc ~ mean_vis_z + Weekday + Month + Year_f, data = daily))
cf <- summary(nb)$coefficients["mean_vis_z", ]
T12 <- data.frame(Model = "Negative binomial: daily count ~ mean visibility (z) + weekday + month + year",
                  IRR_per_1SD = round(unname(exp(cf[1])), 4), CI_low = round(unname(exp(cf[1] - 1.96 * cf[2])), 4),
                  CI_high = round(unname(exp(cf[1] + 1.96 * cf[2])), 4), p_value = format.pval(unname(cf[4]), digits = 3))
show(T12, "Table 12. Association of daily accident count with the daily visibility score (IRR = incidence rate ratio)")
save_tab(T12, "T12_NB_daily_count_vs_visibility")
note("IRR > 1: days whose accidents are recorded in poorer-visibility conditions have proportionally more accidents (per +1 SD of the score).")

# 7c. high-fog days
thr <- max(1, quantile(daily$n_fog, FOG_DAY_Q))
daily$FogDay <- factor(ifelse(daily$n_fog >= thr, "High-fog day", "Other day"), levels = c("Other day", "High-fog day"))
note("'High-fog day' = at least ", thr, " fog/mist-coded accidents (top ", round(100 * (1 - FOG_DAY_Q)), "% of days). Days flagged: ", sum(daily$FogDay == "High-fog day"))
note("To avoid a circular comparison, the outcome is the count of NON-fog-coded accidents (n_acc - n_fog).")
T13 <- do.call(rbind, lapply(split(daily$non_fog, daily$FogDay), function(z)
  data.frame(Days = length(z), Mean = round(mean(z), 1), Median = median(z), SD = round(sd(z), 1))))
T13 <- cbind(Group = rownames(T13), T13, row.names = NULL); show(T13, "Table 13. Non-fog accident records per day")
tt <- t.test(non_fog ~ FogDay, data = daily); wt <- wilcox.test(non_fog ~ FogDay, data = daily, exact = FALSE)
note("Welch t-test p = ", format.pval(tt$p.value, digits = 3), ";  Wilcoxon p = ", format.pval(wt$p.value, digits = 3),
     ";  Cohen's d = ", round(cohens_d(daily$non_fog[daily$FogDay == "High-fog day"], daily$non_fog[daily$FogDay == "Other day"]), 3))
nb2 <- suppressWarnings(MASS::glm.nb(non_fog ~ FogDay + Weekday + Month + Year_f, data = daily))
cf2 <- summary(nb2)$coefficients["FogDayHigh-fog day", ]
T14 <- data.frame(Model = "NB: non-fog daily count ~ high-fog day + weekday + month + year", IRR = round(unname(exp(cf2[1])), 4),
                  CI_low = round(unname(exp(cf2[1] - 1.96 * cf2[2])), 4), CI_high = round(unname(exp(cf2[1] + 1.96 * cf2[2])), 4),
                  p_value = format.pval(unname(cf2[4]), digits = 3))
show(T14, "Table 14. Adjusted comparison, high-fog days vs other days"); save_tab(T13, "T13_high_fog_day_summary"); save_tab(T14, "T14_NB_high_fog_days")
save_plot(ggplot(daily, aes(FogDay, non_fog, fill = FogDay)) + geom_boxplot(outlier.alpha = 0.3, show.legend = FALSE) +
            scale_fill_manual(values = c("grey75", "#9ecae1")) +
            labs(title = "Non-fog accident records per day: high-fog days vs other days", x = NULL, y = "Accident records per day"),
          "P14_boxplot_fog_days")


# ============ 8. (OPTIONAL) REAL NUMERIC VISIBILITY, IF AVAILABLE ============
if (!is.null(VIS_NUM_COL)) {
  section(paste0("8. NUMERIC VISIBILITY COLUMN: ", VIS_NUM_COL))
  vis <- suppressWarnings(as.numeric(d[[VIS_NUM_COL]])); vis[!is.finite(vis) | vis < 0] <- NA
  note("Non-missing values: ", sum(!is.na(vis)), " (", round(100 * mean(!is.na(vis)), 1), "%). CHECK the units and any sentinel codes (e.g. 9999) before trusting this.")
  show(summary(vis), "Summary of numeric visibility:"); show(round(quantile(vis, c(.01, .05, .25, .5, .75, .95, .99), na.rm = TRUE), 3), "Quantiles:")
  vdf <- data.frame(vis = vis, Severity = d$Severity, KSI = d$KSI)
  save_plot(ggplot(vdf[!is.na(vdf$vis), ], aes(vis)) + geom_histogram(bins = 50, fill = "#4d4d4d") + scale_y_continuous(labels = comma) +
              labs(title = "Distribution of numeric visibility", x = VIS_NUM_COL, y = "Records"), "P15_hist_numeric_visibility")
  save_plot(ggplot(vdf[!is.na(vdf$vis), ], aes(Severity, vis, fill = Severity)) + geom_boxplot(outlier.alpha = 0.1, show.legend = FALSE) +
              scale_fill_manual(values = SEV_COLS) + labs(title = "Visibility by accident severity", y = VIS_NUM_COL), "P16_box_visibility_by_severity")
  kwv <- kruskal.test(vis ~ Severity, data = vdf); note("Kruskal-Wallis, visibility across severity: H = ", round(kwv$statistic, 1), ", p = ", format.pval(kwv$p.value, digits = 3))
  sp <- cor.test(vdf$vis, vdf$KSI, method = "spearman", exact = FALSE); note("Spearman rho (visibility vs KSI 0/1) = ", round(sp$estimate, 4), ", p = ", format.pval(sp$p.value, digits = 3))
  ix <- which(!is.na(vis)); ix <- ix[sample.int(length(ix), min(length(ix), 500000))]
  ds <- d[ix, ]; ds$vis_z <- as.numeric(scale(vis[ix]))
  g0 <- glm(KSI ~ vis_z, family = binomial, data = ds)
  g1 <- glm(KSI ~ vis_z + Speed_f + Area + Road_Type_f + Veh_cat + Weekend + Year_f + Hour_band, family = binomial, data = ds)
  for (nm in c("g0", "g1")) { s <- summary(get(nm))$coefficients["vis_z", ]
    note(ifelse(nm == "g0", "Crude", "Adjusted"), " OR for KSI per +1 SD of visibility: ", round(exp(s[1]), 4),
         " [", round(exp(s[1] - 1.96 * s[2]), 4), ", ", round(exp(s[1] + 1.96 * s[2]), 4), "], p = ", format.pval(s[4], digits = 3)) }
  dv <- aggregate(vis ~ Date2, data = data.frame(Date2 = d$Date2, vis = vis), FUN = mean)
  dn <- aggregate(n ~ Date2, data = data.frame(Date2 = d$Date2, n = 1), FUN = sum)
  dv <- merge(dv, dn, by = "Date2")
  save_plot(ggplot(dv, aes(vis, n)) + geom_point(alpha = 0.25) + geom_smooth(method = "loess", formula = y ~ x, colour = "#b2182b") +
              labs(title = "Daily accident count vs daily mean numeric visibility", x = paste("Daily mean", VIS_NUM_COL), y = "Accident records"),
            "P17_scatter_daily_count_vs_numeric_visibility")
  note("Daily Spearman rho (mean visibility vs count) = ", round(cor(dv$vis, dv$n, method = "spearman", use = "complete.obs"), 3))
} else {
  section("8. NUMERIC VISIBILITY COLUMN")
  note("Skipped - no numeric visibility column available. (Sections 2-7 use the Light/Weather proxies.)")
}


# ============================ 9. KEY RESULTS SUMMARY ==========================
section("9. KEY NUMBERS TO QUOTE IN YOUR REPORT")
g <- function(tab, var, grp, col) tab[tab$Variable == var & tab$Group == grp, col]
note("Records analysed (after de-duplication): ", format(nrow(d), big.mark = ","))
note("Fog/Mist share of all records: ", round(100 * mean(d$Weather_Group == "Fog/Mist"), 2), "%;  Dark (any) share of records with known lighting: ",
     round(100 * mean(startsWith(as.character(d$Light_Group[d$Light_Group != "Unknown"]), "Dark")), 1), "%")
note("KSI rate:  Fine ", g(T2, "Weather", "Fine", "Pct_KSI"), "%  vs  Fog/Mist ", g(T2, "Weather", "Fog/Mist", "Pct_KSI"),
     "%;   Daylight ", g(T2, "Lighting", "Daylight", "Pct_KSI"), "%  vs  Dark-unlit/none ", g(T2, "Lighting", "Dark: unlit / no lighting", "Pct_KSI"), "%")
best_or <- OR_TAB[OR_TAB$Outcome == "Fatal or Serious vs Slight" & OR_TAB$Model == "Adj-A", ]
for (i in seq_len(nrow(best_or))) note("Adj-A OR (KSI) ", best_or$Term[i], ": ", round(best_or$OR[i], 3), " [", round(best_or$CI_low[i], 3), ", ", round(best_or$CI_high[i], 3), "]")
note("Cochran-Armitage trend (KSI vs score): p = ", format.pval(trend_ksi$p.value, digits = 3))
note("Daily NB IRR per +1 SD visibility-score: ", T12$IRR_per_1SD, " [", T12$CI_low, ", ", T12$CI_high, "]")

writeLines(LOG, file.path(OUT_DIR, "Member4_full_console_log.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_package_versions.txt"))
cat("\nDONE. Everything was saved in: ", normalizePath(OUT_DIR), "\n")
cat("  - PNG plots (P1..P14), CSV tables (T0..T14), Member4_full_console_log.txt, sessionInfo_package_versions.txt\n")
