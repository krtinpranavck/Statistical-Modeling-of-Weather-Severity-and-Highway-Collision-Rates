# Member 6: Traffic and Time Factors
# Usage: Rscript member_6_traffic_time_analysis.R path/to/UK_Accident.csv output_folder
# Requires base/recommended R only (MASS ships with most R installations).
# Source files are never modified. All counts describe reported injury accidents.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "UK_Accident.csv"
out <- if (length(args) >= 2) args[2] else "member6_results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)
options(stringsAsFactors = FALSE, scipen = 7)
Sys.setlocale("LC_TIME", "C")
if (!requireNamespace("MASS", quietly = TRUE)) stop("Install the recommended MASS package.")
write_tab <- function(x, name) write.csv(x, file.path(out, "tables", paste0(name, ".csv")), row.names = FALSE, na = "")
cat("Reading source CSV...\n")
header <- names(read.csv(input, nrows = 0, check.names = FALSE))
needed <- c("Accident_Index", "Year", "Date", "Time", "Location_Easting_OSGR",
 "Location_Northing_OSGR", "Accident_Severity", "Number_of_Vehicles",
 "Number_of_Casualties", "Day_of_Week", "Light_Conditions", "1st_Road_Class")
stopifnot(all(needed %in% header))
classes <- rep("NULL", length(header)); names(classes) <- header
classes[needed] <- "character"
raw <- read.csv(input, colClasses = classes, check.names = FALSE, na.strings = c("", "NA"))
for (nm in c("Year", "Accident_Severity", "Number_of_Vehicles", "Number_of_Casualties", "Day_of_Week", "1st_Road_Class")) raw[[nm]] <- as.integer(raw[[nm]])
key <- c("Year", "Accident_Index", "Date", "Time", "Location_Easting_OSGR",
 "Location_Northing_OSGR", "Accident_Severity", "Number_of_Vehicles", "Number_of_Casualties")
duplicate <- duplicated(raw[key])
years <- sort(unique(raw$Year))
audit <- data.frame(year = years, raw_records = as.integer(table(factor(raw$Year, levels = years))),
 duplicates_removed = as.integer(table(factor(raw$Year[duplicate], levels = years))))
audit$clean_records <- audit$raw_records - audit$duplicates_removed
write_tab(audit, "audit_by_year")
d <- raw[!duplicate, ]; raw_n <- nrow(raw); removed_n <- sum(duplicate)
rm(raw, duplicate); gc(verbose = FALSE)
d$date <- as.Date(d$Date, format = "%d/%m/%Y")
stopifnot(!anyNA(d$date), all(as.integer(format(d$date, "%Y")) == d$Year))
valid_time <- !is.na(d$Time) & grepl("^([01][0-9]|2[0-3]):[0-5][0-9]$", d$Time)
d$hour <- NA_integer_; d$hour[valid_time] <- as.integer(substr(d$Time[valid_time], 1, 2))
d$day_night <- ifelse(grepl("^Daylight", d$Light_Conditions), "Daylight",
 ifelse(grepl("^Dark", d$Light_Conditions), "Darkness", NA_character_))
stopifnot(!anyNA(d$day_night), all(d$Accident_Severity %in% 1:3))
d$day_night <- factor(d$day_night, levels = c("Daylight", "Darkness"))
d$dow <- as.POSIXlt(d$date)$wday
stopifnot(all(d$dow + 1L == d$Day_of_Week))
d$weekend <- factor(ifelse(d$dow %in% c(0, 6), "Weekend", "Weekday"), levels = c("Weekday", "Weekend"))
d$ksi <- d$Accident_Severity %in% c(1, 2)
d$period <- cut(d$hour, breaks = c(-1, 5, 11, 17, 23), labels = c("00-05", "06-11", "12-17", "18-23"))
wilson <- function(k, n) {
 z <- qnorm(.975); den <- 1 + z*z/n
 mid <- (k/n + z*z/(2*n))/den
 half <- z * sqrt(k/n*(1-k/n)/n + z*z/(4*n*n))/den
 cbind(lower = mid-half, upper = mid+half)
}
group_summary <- function(g, levels) {
 good <- !is.na(g); f <- factor(g[good], levels = levels)
 tab <- table(f, factor(d$Accident_Severity[good], levels = 1:3))
 z <- data.frame(group = levels, records = as.integer(rowSums(tab)), fatal = as.integer(tab[, 1]), serious = as.integer(tab[, 2]), slight = as.integer(tab[, 3]))
 z$ksi <- z$fatal + z$serious; z$share_pct <- 100*z$records/sum(z$records)
 z$ksi_pct <- 100*z$ksi/z$records; ci <- wilson(z$ksi, z$records)
 z$ksi_lower_pct <- 100*ci[,1]; z$ksi_upper_pct <- 100*ci[,2]
 z$fatal_pct <- 100*z$fatal/z$records
 z
}
hourly <- group_summary(d$hour, 0:23)
dn <- group_summary(d$day_night, c("Daylight", "Darkness"))
period <- group_summary(d$period, levels(d$period))
week <- group_summary(d$weekend, c("Weekday", "Weekend"))
dow_levels <- c(1:6,0)
dow <- group_summary(d$dow, dow_levels)
dow$group <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
write_tab(hourly, "hourly_summary"); write_tab(dn, "day_night_summary")
write_tab(period, "six_hour_summary"); write_tab(week, "weekday_weekend_summary")
# Vehicles involved is an accident descriptor, not a traffic-flow estimate.
d$vehicle_group <- ifelse(is.na(d$Number_of_Vehicles) | d$Number_of_Vehicles < 1, NA_character_,
 ifelse(d$Number_of_Vehicles >= 4, "4+", as.character(d$Number_of_Vehicles)))
d$vehicle_group <- factor(d$vehicle_group, levels=c("1","2","3","4+"))
vehicles <- group_summary(d$vehicle_group, levels(d$vehicle_group))
write_tab(vehicles,"vehicles_involved_summary")

# Complete calendar only within years actually supplied. 2008 is not made zero.
calendar <- do.call(c, lapply(years, function(y) seq(as.Date(paste0(y,"-01-01")), as.Date(paste0(y,"-12-31")), by="day")))
daily <- data.frame(date = calendar)
daily$count <- as.integer(table(factor(as.character(d$date), levels = as.character(calendar))))
daily$dow <- factor(as.POSIXlt(calendar)$wday, levels = dow_levels, labels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))
daily$weekend <- factor(ifelse(daily$dow %in% c("Saturday", "Sunday"), "Weekend", "Weekday"), levels = c("Weekday", "Weekend"))
daily$month <- factor(as.integer(format(calendar, "%m")), levels=1:12)
daily$year <- factor(format(calendar, "%Y"))
daily$daylight <- as.integer(table(factor(as.character(d$date[d$day_night == "Daylight"]), levels=as.character(calendar))))
daily$darkness <- daily$count - daily$daylight
daily$timed <- as.integer(table(factor(as.character(d$date[valid_time]), levels=as.character(calendar))))
stopifnot(sum(daily$count) == nrow(d), sum(daily$timed) == sum(valid_time))
dow$calendar_days <- as.integer(table(daily$dow))
dow$mean_per_calendar_day <- dow$records/dow$calendar_days
week$calendar_days <- as.integer(table(daily$weekend))
week$mean_per_calendar_day <- week$records/week$calendar_days
week$sd_daily_count <- as.numeric(tapply(daily$count,daily$weekend,sd))
write_tab(dow, "day_of_week_summary"); write_tab(week, "weekday_weekend_summary")
write_tab(daily, "daily_counts")

# Severity association. Model-based intervals assume independent records.
tests <- list()
for (nm in c("day_night", "period", "weekend", "vehicle_group")) {
 tab <- table(d[[nm]], d$Accident_Severity)
 ct <- chisq.test(tab, correct = FALSE)
 tests[[nm]] <- data.frame(test=nm, n=sum(tab), chi_square=as.numeric(ct$statistic), df=as.numeric(ct$parameter), p=ct$p.value,
 cramers_v=sqrt(as.numeric(ct$statistic)/(sum(tab)*min(nrow(tab)-1,ncol(tab)-1))), minimum_expected=min(ct$expected))
}
write_tab(do.call(rbind, tests), "severity_chi_square")
contrasts <- list()
for (nm in c("day_night", "weekend")) {
 z <- if(nm == "day_night") dn else week
 a <- z$ksi[2]; b <- z$records[2]-a; c0 <- z$ksi[1]; e <- z$records[1]-c0
 lor <- log(a/b/(c0/e)); se <- sqrt(1/a+1/b+1/c0+1/e)
 contrasts[[nm]] <- data.frame(contrast=paste(z$group[2],"vs",z$group[1]), odds_ratio=exp(lor), lower=exp(lor-1.96*se), upper=exp(lor+1.96*se), difference_percentage_points=z$ksi_pct[2]-z$ksi_pct[1])
}
write_tab(do.call(rbind,contrasts),"ksi_crude_odds_ratios")

# A whole-week cluster bootstrap assesses record clustering and serial dependence.
# The point estimate remains the full-sample crude odds ratio. Fixed seed.
week_id <- as.character(d$date - ((d$dow+6) %% 7))
bootdata <- data.frame(week=week_id,daynight=as.character(d$day_night),ksi=as.integer(d$ksi),nonksi=as.integer(!d$ksi))
weekly <- aggregate(cbind(ksi,nonksi) ~ week+daynight, bootdata, sum)
wks <- sort(unique(weekly$week)); wb <- matrix(0,nrow=length(wks),ncol=4)
for (j in 1:2) {
 z <- weekly[weekly$daynight==c("Daylight","Darkness")[j],]
 wb[match(z$week,wks),2*j-1] <- z$ksi; wb[match(z$week,wks),2*j] <- z$nonksi
}
set.seed(62026); boot_or <- replicate(1999, {s<-colSums(wb[sample.int(nrow(wb),nrow(wb),replace=TRUE),]); (s[3]/s[4])/(s[1]/s[2])})
write_tab(data.frame(contrast="Darkness vs Daylight", bootstrap_units=nrow(wb), replicates=1999,
 odds_ratio=contrasts$day_night$odds_ratio, lower=unname(quantile(boot_or,.025)),upper=unname(quantile(boot_or,.975))),"ksi_week_cluster_bootstrap")
rm(bootdata,weekly,week_id); gc(verbose=FALSE)

# Daily frequency: Poisson check then negative binomial, adjusted for season/year.
cat("Fitting daily count models...\n")
pois <- glm(count ~ weekend + month + year, family=poisson(), data=daily)
disp <- sum(residuals(pois,type="pearson")^2)/df.residual(pois)
nb <- MASS::glm.nb(count ~ weekend + month + year, data=daily)
nb_dow <- MASS::glm.nb(count ~ dow + month + year, data=daily)
coef_table <- function(model) {
 cm <- coef(summary(model)); data.frame(term=rownames(cm), coefficient=cm[,1], se=cm[,2], ratio=exp(cm[,1]), lower=exp(cm[,1]-qnorm(.975)*cm[,2]), upper=exp(cm[,1]+qnorm(.975)*cm[,2]), p=cm[,4], row.names=NULL)
}
write_tab(coef_table(nb),"negative_binomial_weekend")
write_tab(coef_table(nb_dow),"negative_binomial_weekday")

# Newey-West HAC sandwich SE for NB regression coefficients, conditional on theta.
# Chronological lags never join across the missing year; lag 7 captures weekly dependence.
nb_hac <- function(model, dat, maxlag=7) {
 X<-model.matrix(model); mu<-fitted(model); theta<-model$theta; y<-model$y
 score<-X*as.numeric((y-mu)/(1+mu/theta)); weights<-as.numeric(mu/(1+mu/theta))
 bread<-solve(crossprod(X,X*weights)); meat<-crossprod(score)
 n<-nrow(score)
 for (lag in 1:maxlag) {
   ii<-(lag+1):n; jj<-1:(n-lag)
   keep<-as.integer(dat$date[ii]-dat$date[jj])==lag
   S<-crossprod(score[ii[keep],,drop=FALSE],score[jj[keep],,drop=FALSE])
   meat<-meat+(1-lag/(maxlag+1))*(S+t(S))
 }
 V<-bread%*%meat%*%bread*n/(n-ncol(X)); se<-sqrt(diag(V)); b<-coef(model)
 data.frame(term=names(b),coefficient=b,se=se,ratio=exp(b),lower=exp(b-qnorm(.975)*se),upper=exp(b+qnorm(.975)*se),p=2*pnorm(-abs(b/se)),row.names=NULL)
}
write_tab(nb_hac(nb,daily),"negative_binomial_weekend_hac7")
write_tab(nb_hac(nb,daily,14),"negative_binomial_weekend_hac14")
write_tab(nb_hac(nb_dow,daily),"negative_binomial_weekday_hac7")
write_tab(data.frame(poisson_dispersion=disp,nb_theta=nb$theta,nb_pearson_dispersion=sum(residuals(nb,type="pearson")^2)/df.residual(nb),nb_converged=nb$converged,days=nrow(daily),minimum_daily_count=min(daily$count)),"count_model_diagnostics")

# All roads form the primary analysis, consistent with Section 6.2.
# A road-class restriction must be a separate, team-agreed analysis.
scope <- group_summary(d[["1st_Road_Class"]],sort(unique(d[["1st_Road_Class"]])))
write_tab(scope,"road_class_scope")
flow <- data.frame(metric=c("raw_records","duplicates_removed","clean_records","invalid_time","invalid_date","covered_days","time_valid_records","unknown_lighting","year_2008_records"),
 value=c(raw_n,removed_n,nrow(d),sum(!valid_time),sum(is.na(d$date)),nrow(daily),sum(valid_time),sum(is.na(d$day_night)),sum(d$Year==2008)))
write_tab(flow,"analysis_flow")
write_tab(data.frame(metric="invalid_number_of_vehicles",value=sum(is.na(d$vehicle_group))),"vehicle_data_quality")

# Publication charts, base R. Uniform navy/teal palette and labelled denominators.
navy <- "#214B6B"; teal <- "#2C807B"; grey <- "#DDE4E9"
start_png <- function(name,w=1800,h=1000) {png(file.path(out,"figures",paste0(name,".png")),width=w,height=h,res=200); par(mar=c(4.5,4.8,2.5,1.2),las=1,cex=.9,col.axis="#303030",col.lab="#303030",bty="l")}
start_png("Figure_6_4_1_hourly_counts")
bp<-barplot(hourly$records/1000,names.arg=sprintf("%02d",0:23),col=navy,border=NA,xlab="Hour of recorded accident",ylab="Accident records (thousands)",main="Accidents by hour of day",cex.names=.8)
dev.off()
start_png("Figure_6_4_2_day_night",1800,950)
par(mfrow=c(1,2),mar=c(4.5,4.5,2.5,1))
b<-barplot(dn$records/1000,names.arg=dn$group,col=c(navy,teal),border=NA,ylab="Accident records (thousands)",ylim=c(0,max(dn$records/1000)*1.18),main="Recorded frequency")
text(b,dn$records/1000,labels=paste0(formatC(dn$share_pct,digits=1,format="f"),"%"),pos=3)
b<-barplot(dn$ksi_pct,names.arg=dn$group,col=c(navy,teal),border=NA,ylab="Fatal or serious accidents (%)",ylim=c(0,max(dn$ksi_upper_pct)*1.2),main="Severity among recorded accidents")
arrows(b,dn$ksi_lower_pct,b,dn$ksi_upper_pct,angle=90,code=3,length=.04)
text(b,dn$ksi_upper_pct,labels=paste0(formatC(dn$ksi_pct,digits=2,format="f"),"%"),pos=3)
dev.off()
start_png("Figure_6_4_3_weekday_counts")
b<-barplot(dow$mean_per_calendar_day,names.arg=substr(dow$group,1,3),col=c(rep(navy,5),rep(teal,2)),border=NA,ylab="Mean accident records per calendar day",ylim=c(0,max(dow$mean_per_calendar_day)*1.15),main="Calendar-adjusted accident counts by weekday")
text(b,dow$mean_per_calendar_day,round(dow$mean_per_calendar_day,1),pos=3,cex=.85)
dev.off()
start_png("Figure_6_4_4_hourly_severity")
plot(0:23,hourly$ksi_pct,type="b",pch=16,col=teal,xlab="Hour of recorded accident",ylab="Fatal or serious accidents (%)",xaxt="n",ylim=range(c(hourly$ksi_lower_pct,hourly$ksi_upper_pct)),main="Severity share by hour of day")
axis(1,at=0:23,labels=sprintf("%02d",0:23),cex.axis=.8)
arrows(0:23,hourly$ksi_lower_pct,0:23,hourly$ksi_upper_pct,angle=90,code=3,length=.025,col=teal)
dev.off()
daily$ym <- format(daily$date,"%Y-%m")
monthly <- aggregate(count~ym,daily,mean)
monthly$date <- as.Date(paste0(monthly$ym,"-15"))
write_tab(monthly,"monthly_mean_daily_counts")
start_png("Figure_6_4_5_monthly_trend")
plot(range(daily$date),range(monthly$count)*c(.95,1.05),type="n",xlab="Calendar year",ylab="Mean daily accident records within month",main="Monthly trend across supplied years")
for (run in list(2005:2007,2009:2014)) {q<-monthly[as.integer(substr(monthly$ym,1,4))%in%run,];lines(q$date,q$count,col=navy,lwd=2)}
rect(as.Date("2008-01-01"),par("usr")[3],as.Date("2008-12-31"),par("usr")[4],col=grey,border=NA)
text(as.Date("2008-07-01"),mean(par("usr")[3:4]),"2008\nnot supplied",cex=.9)
dev.off()
annual <- aggregate(count~year,daily,mean); annual$days<-as.integer(table(daily$year)); annual$total<-as.integer(tapply(daily$count,daily$year,sum))
write_tab(annual,"annual_counts")
start_png("Figure_6_4_6_vehicle_involvement",1800,950)
par(mfrow=c(1,2),mar=c(4.5,4.5,2.5,1))
b<-barplot(vehicles$records/1000,names.arg=vehicles$group,col=navy,border=NA,xlab="Vehicles involved in the accident",ylab="Accident records (thousands)",main="Vehicle involvement counts",ylim=c(0,max(vehicles$records/1000)*1.18))
text(b,vehicles$records/1000,labels=paste0(formatC(vehicles$share_pct,digits=1,format="f"),"%"),pos=3)
b<-barplot(vehicles$ksi_pct,names.arg=vehicles$group,col=teal,border=NA,xlab="Vehicles involved in the accident",ylab="Fatal or serious accidents (%)",main="Severity by vehicle involvement",ylim=c(0,max(vehicles$ksi_upper_pct)*1.18))
arrows(b,vehicles$ksi_lower_pct,b,vehicles$ksi_upper_pct,angle=90,code=3,length=.04)
text(b,vehicles$ksi_upper_pct,labels=paste0(formatC(vehicles$ksi_pct,digits=1,format="f"),"%"),pos=3)
dev.off()

capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
capture.output(list(flow=flow,daynight=dn,week=week,period=period,weekday=dow,tests=do.call(rbind,tests),odds=do.call(rbind,contrasts),poisson_dispersion=disp,nb=summary(nb),hac=nb_hac(nb,daily)),file=file.path(out,"analysis_log.txt"))
cat("Complete:",format(nrow(d),big.mark=","),"records;",nrow(daily),"days;",sum(!valid_time),"invalid times.\n")
cat("Output folder:",normalizePath(out),"\n")
