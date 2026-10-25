data <- read.csv(file.choose())

dim(data)
names(data)
str(data)
summary(data)

# ============================================================
# 1. LOAD REQUIRED PACKAGE
# ============================================================
# install.packages("ggplot2")    <- run this once to install package
library(ggplot2)

# ============================================================
# 2. BASIC DATASET INFORMATION
# ============================================================

cat("============================================================\n")
cat("DATASET OVERVIEW\n")
cat("============================================================\n")

cat("Number of rows:", nrow(data), "\n")
cat("Number of columns:", ncol(data), "\n")

cat("Unique accidents:",
    length(unique(data$Accident_Index)), "\n")

cat("Study years:",
    min(data$Year, na.rm = TRUE),
    "to",
    max(data$Year, na.rm = TRUE), "\n")

cat("\nColumn names:\n")
print(names(data))

# ============================================================
# 3. REMOVE UNNECESSARY IMPORTED INDEX COLUMN
# ============================================================

# X is only the row/index column created during CSV import.
# It is not an accident variable.

data$X <- NULL

# ============================================================
# 4. CHECK MISSING VALUES
# ============================================================

cat("\nMISSING VALUE SUMMARY\n")
missing_values <- colSums(is.na(data))
print(
  missing_values[missing_values > 0]
)


# ============================================================
# 5. CONVERT DATE AND TIME
# ============================================================

# Date is stored as DD/MM/YYYY
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
# Time is stored as HH:MM.
# Some records have missing/blank time values.
data$Hour <- suppressWarnings(
  as.numeric(substr(data$Time, 1, 2))
)

# ============================================================
# 6. BASIC ACCIDENT DESCRIPTIVE STATISTICS
# ============================================================
cat("\nACCIDENT CHARACTERISTICS\n")
cat("\nNumber of Vehicles:\n")
print(summary(data$Number_of_Vehicles))

cat("\nNumber of Casualties:\n")
print(summary(data$Number_of_Casualties))

cat("\nSpeed Limit:\n")
print(summary(data$Speed_limit))

# ============================================================
# 7. DETAILED DESCRIPTIVE STATISTICS TABLE
# ============================================================

descriptive_stats <- data.frame(
  Variable = c(
    "Number of Vehicles",
    "Number of Casualties",
    "Speed Limit"
  ),
  
  Mean = c(
    mean(data$Number_of_Vehicles, na.rm = TRUE),
    mean(data$Number_of_Casualties, na.rm = TRUE),
    mean(data$Speed_limit, na.rm = TRUE)
  ),
  
  Median = c(
    median(data$Number_of_Vehicles, na.rm = TRUE),
    median(data$Number_of_Casualties, na.rm = TRUE),
    median(data$Speed_limit, na.rm = TRUE)
  ),
  
  SD = c(
    sd(data$Number_of_Vehicles, na.rm = TRUE),
    sd(data$Number_of_Casualties, na.rm = TRUE),
    sd(data$Speed_limit, na.rm = TRUE)
  ),
  
  Minimum = c(
    min(data$Number_of_Vehicles, na.rm = TRUE),
    min(data$Number_of_Casualties, na.rm = TRUE),
    min(data$Speed_limit, na.rm = TRUE)
  ),
  
  Maximum = c(
    max(data$Number_of_Vehicles, na.rm = TRUE),
    max(data$Number_of_Casualties, na.rm = TRUE),
    max(data$Speed_limit, na.rm = TRUE)
  )
)

descriptive_stats[, 2:6] <- round(descriptive_stats[, 2:6], 3)
cat("\nDescriptive Statistics Table:\n")
print(descriptive_stats)

# ============================================================
# 8. ACCIDENT SEVERITY FREQUENCY TABLE
# ============================================================

cat("\nACCIDENT SEVERITY\n")
severity_freq <- table(data$Accident_Severity)
severity_percent <- prop.table(severity_freq) * 100
severity_table <- data.frame(
  Severity = names(severity_freq),
  Frequency = as.vector(severity_freq),
  Percentage = round(as.vector(severity_percent),2)
)
print(severity_table)

# ============================================================
# 9. SEVERITY BAR CHART
# ============================================================

ggplot(data,aes(x = factor(Accident_Severity))) + geom_bar() +
  labs(title = "Distribution of Accident Severity", x = "Accident Severity", y = "Number of Accident Records") +
  theme_minimal()

# ============================================================
# 10. ACCIDENTS BY YEAR
# ============================================================

cat("\nACCIDENTS BY YEAR\n")
year_freq <- table(data$Year)
print(year_freq)
# Year bar chart

ggplot(data,aes(x = factor(Year))) +
  geom_bar() +
  labs(title = "Number of Accident Records by Year",x = "Year",y = "Number of Accident Records") +
  theme_minimal()

# ============================================================
# 11. ACCIDENTS BY DAY OF WEEK
# ============================================================

cat("\nACCIDENTS BY DAY OF WEEK\n")
day_freq <- table(data$Day_of_Week)
print(day_freq)
# Day-of-week plot

ggplot(data,aes(x = factor(Day_of_Week))) +
  geom_bar() +
  labs(title = "Accident Records by Day of Week",x = "Day of Week",y = "Number of Accident Records") +
  theme_minimal()
# ============================================================
# 12. ACCIDENTS BY HOUR
# ============================================================

cat("\nACCIDENTS BY HOUR\n")

hour_freq <- table(data$Hour,useNA = "ifany")
print(hour_freq)
# Remove missing hours for plotting
hour_data <- data[!is.na(data$Hour),]

ggplot(hour_data,aes(x = Hour)) +
  geom_bar() +
  scale_x_continuous(
    breaks = 0:23
  ) +
  labs(
    title = "Accident Records by Hour of Day",x = "Hour of Day",y = "Number of Accident Records"
  ) +
  theme_minimal()

# ============================================================
# 13. CREATE TIME PERIODS
# ============================================================

data$Time_Period <- cut(
  data$Hour,
  breaks = c(-1,5,11,17,23),
  labels = c("Night","Morning","Afternoon","Evening")
)

cat("\nACCIDENTS BY TIME PERIOD\n")

time_period_freq <- table(data$Time_Period,useNA = "ifany")
print(time_period_freq)

# Time-period plot

ggplot(data,aes(x = Time_Period)) +
  geom_bar() +
  labs(title = "Accident Records by Time Period",x = "Time Period",y = "Number of Accident Records"
  ) +
  theme_minimal()


# ============================================================
# 14. DAILY ACCIDENT FREQUENCY
# ============================================================

cat("\nDAILY ACCIDENT FREQUENCY\n")
daily_accidents <- aggregate(
  Accident_Index ~ Date,
  data = data,
  FUN = length
)

names(daily_accidents)[2] <-"Accident_Records"
head(daily_accidents)
# Daily time-series plot
ggplot(daily_accidents,
  aes(x = Date,y = Accident_Records)
) +
  geom_line() +
  labs(title = "Daily Accident Frequency",x = "Date",y = "Number of Accident Records"
  ) +
  theme_minimal()

# ============================================================
# 15. MONTHLY ACCIDENT FREQUENCY
# ============================================================

data$Month <- format(data$Date,"%Y-%m")

monthly_accidents <- as.data.frame(table(data$Month))

names(monthly_accidents) <- c("Month","Accident_Records")
ggplot(monthly_accidents,
  aes(x = Month,y = Accident_Records,group = 1)
) +
  geom_line() +
  labs(title = "Monthly Accident Frequency",x = "Month",y = "Number of Accident Records"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)
  )

# ============================================================
# 16. LOCATION ANALYSIS
# ============================================================

cat("LOCATION ANALYSIS\n")

cat("\nNumber of unique districts:\n")
print(length(unique(data$Local_Authority_.District.)))
# Top 10 districts
district_freq <- sort(
  table(data$Local_Authority_.District.),
  decreasing = TRUE
)

cat("\nTop 10 districts by accident records:\n")
print(
  head(district_freq, 10)
)

# ============================================================
# 17. URBAN / RURAL PATTERN
# ============================================================
cat("\nURBAN / RURAL DISTRIBUTION\n")

urban_rural_freq <- table(data$Urban_or_Rural_Area)
print(urban_rural_freq)
urban_rural_percent <-prop.table(urban_rural_freq) * 100
print(round(urban_rural_percent,2))
ggplot(data,aes(x = factor(Urban_or_Rural_Area))
) +
  geom_bar() +
  labs(
    title = "Accident Records by Urban or Rural Area",x = "Urban / Rural Area Code",y = "Number of Accident Records"
  ) +
  theme_minimal()

# ============================================================
# 18. ROAD TYPE DISTRIBUTION
# ============================================================

cat("\nROAD TYPE DISTRIBUTION\n")
road_type_freq <- sort(table(data$Road_Type),decreasing = TRUE)
print(road_type_freq)
ggplot(data,aes(x = reorder(Road_Type,Road_Type,FUN = function(x) -length(x)))
) +
  geom_bar() +
  labs(
    title = "Accident Records by Road Type",x = "Road Type",y = "Number of Accident Records"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text( angle = 45, hjust = 1))

# ============================================================
# 20. NUMBER OF VEHICLES DISTRIBUTION
# ============================================================

ggplot(
  data,
  aes(
    x = Number_of_Vehicles
  )
) +
  
  geom_histogram(
    binwidth = 1
  ) +
  
  labs(
    title = "Distribution of Number of Vehicles",
    x = "Number of Vehicles",
    y = "Frequency"
  ) +
  
  theme_minimal()

# ============================================================
# 21. NUMBER OF CASUALTIES DISTRIBUTION
# ============================================================

ggplot(
  data,
  aes(
    x = Number_of_Casualties
  )
) +
  
  geom_histogram(
    binwidth = 1
  ) +
  
  labs(
    title = "Distribution of Number of Casualties",
    x = "Number of Casualties",
    y = "Frequency"
  ) +
  
  theme_minimal()
# ============================================================
# 22. FINAL SUMMARY VALUES
# ============================================================

cat("\nMEMBER 2 SUMMARY\n")

cat(
  "\nTotal records:",
  nrow(data)
)

cat(
  "\nUnique accidents:",
  length(unique(data$Accident_Index))
)

cat(
  "\nStudy period:",
  min(data$Year, na.rm = TRUE),
  "-",
  max(data$Year, na.rm = TRUE)
)

cat(
  "\nAverage vehicles per record:",
  round(
    mean(
      data$Number_of_Vehicles,
      na.rm = TRUE
    ),
    3
  )
)

cat(
  "\nAverage casualties per record:",
  round(
    mean(
      data$Number_of_Casualties,
      na.rm = TRUE
    ),
    3
  )
)

cat(
  "\nAverage speed limit:",
  round(
    mean(
      data$Speed_limit,
      na.rm = TRUE
    ),
    2
  )
)

cat("\n\nMost frequent severity:\n")
print(
  names(
    which.max(severity_freq)
  )
)

cat("\nYear with most accident records:\n")
print(
  names(
    which.max(year_freq)
  )
)

cat("\nDay with most accident records:\n")
print(
  names(
    which.max(day_freq)
  )
)

cat("\nHour with most accident records:\n")
print(
  names(
    which.max(hour_freq)
  )
)

cat("\nMost frequent road type:\n")
print(
  names(
    which.max(road_type_freq)
  )
)

