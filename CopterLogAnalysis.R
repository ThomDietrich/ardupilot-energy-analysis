#!/usr/bin/Rscript

# Clear workspace
remove(list = ls())

#####################################################################
# Logfiles to analyze

logfiles = data.frame()

# CustomCopter
###logfiles <- rbind(logfiles, c("2016-02-28 15-56-34.log", "2016-02-28 15-56-34 CC 'Quadrat'", "", ""), stringsAsFactors = FALSE)
###logfiles <- rbind(logfiles, c("2016-02-28 17-27-54.log", "2016-02-28 17-27-54 CC '23' ersterTest", "", ""), stringsAsFactors = FALSE)
#
logfiles <- rbind(logfiles, c("2016-03-11 11-51-26.log", "2016-03-11 11-51-26 CC '23'", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-03-11 12-01-04.log", "2016-03-11 12-01-04 CC '23'", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-03-11 12-09-31.log", "2016-03-11 12-09-31 CC '23' grAkku", 1, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-03-11 15-24-14.log", "2016-03-11 15-24-14 CC '23' grGeschw", 1, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-03-11 15-32-57.log", "2016-03-11 15-32-57 CC '23' klGeschw", 1, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-03-11 15-46-32.log", "2016-03-11 15-46-32 CC '23' grGewicht", 1, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-04-30 17-45-40.log", "2016-04-30 17-45-40 CC 'Steig'", 2, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-04-30 18-03-59.log", "2016-04-30 18-03-59 CC 'Steig'", 2, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-05-01 14-12-24.log", "2016-05-01 14-12-24 CC 'Steig' Wind", 3, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 13-44-35.log", "2016-10-23 13-44-35 CC 'Lw-Quadrat'", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 13-57-28.log", "2016-10-23 13-57-28 CC 'Lw-Benchmark'", 4, ""), stringsAsFactors = FALSE)

# SoloCopter
logfiles <- rbind(logfiles, c("2016-10-23 14-12-01.log", "2016-10-23 14-12-01 Solo 'Lw-Quadrat'", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 14-32-45.log", "2016-10-23 14-32-45 Solo 'Lw-Benchmark' Akku a", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 14-41-58.log", "2016-10-23 14-41-58 Solo 'Lw-Benchmark' Akku b(Flug1)", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 14-49-50.log", "2016-10-23 14-49-50 Solo 'Lw-Benchmark' Akku b(Flug2)", 13, ""), stringsAsFactors = FALSE)

logfiles <- rbind(logfiles, c("2016-10-27 17-19-01.log", "2016-10-27 17-19-01 Solo 'Lw-Quadrat' 3/4 Akku", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-25-25.log", "2016-10-27 17-25-25 Solo 'Lw-Quadrat' voller Akku", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-29-25.log", "2016-10-27 17-29-25 Solo 'Lw-Quadrat' 80%", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-33-10.log", "2016-10-27 17-33-10 Solo 'Lw-Benchmark' 67%", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-40-31.log", "2016-10-27 17-40-31 Solo 'Lw-Benchmark2' voller Akku", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-48-39.log", "2016-10-27 17-48-39 Solo 'Lw-Benchmark' bei 57%", 13, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-54-43.log", "2016-10-27 17-54-43 Solo 'Lw-Benchmark2' voller Akku", 13, ""), stringsAsFactors = FALSE)

names(logfiles) <- (c("filename", "description", "class", "weather"))

# class: https://plot.ly/r/reference/#scatter-marker

# Solo Hover/LOITER_UNLIM
hoverlog <- "2016-11-06 16-41-01-hover.log"

#####################################################################
# Settings

minSamplesForMean <- 30
transitionSampleSeconds <- 3

# Normalization
batteryVoltageCustomCopter   <- 3 * 3.7 # 11.1
batteryVoltageSoloCopter <- 4 * 3.7 # 14.8
referenceVoltage         <- 4 * 3.7

plotlyUsername <- ""
plotlyApiKey <- ""

#####################################################################
# Install packages, load libraries

list.of.packages <- c("plotly", "Hmisc", "geosphere")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
remove(list = c("list.of.packages", "new.packages"))

library(plotly)
library(Hmisc)
library(geosphere)
packageVersion('plotly')
packageVersion('Hmisc')
packageVersion('geosphere')

Sys.setenv("plotly_username" = plotlyUsername)
Sys.setenv("plotly_api_key" = plotlyApiKey)
remove(plotlyUsername, plotlyApiKey)

#####################################################################
# Helper functions

#
Rad2deg <- function(rad) {(rad * 180) / (pi)}
Deg2rad <- function(deg) {(deg * pi) / (180)}

# convert ArduCopter command ID to human readable type
CId2type <- function(cid) {
  # reference: https://pixhawk.ethz.ch/mavlink
  cmd.CId.as.String <- as.list(c("WAYPOINT", "LOITER_UNLIM", "LOITER_TIME", "RETURN_TO_LAUNCH", "LAND", "TAKEOFF"))
  names(cmd.CId.as.String) <-  c(        16,             17,            19,                 20,     21,        22)
  if (! exists(cid, where=cmd.CId.as.String)){
    warning(paste("CId", cid, "not in the list of known ArduPilot copter mission commands. Please add..."))
    return("-unknown-")
  }
  return(cmd.CId.as.String[cid])
}

#####################################################################
# Function: Parse ArduCopter Logfile, return extracted data frames

ParseLogdata <- function(data) {
  copter.model <- if('STRT' %in% data$V1) "CustomCopter" else "3DRSolo"

  if (copter.model == "CustomCopter") {
    data.starttime <- data[data$V1 %in% c('STRT'), 2]
    time.factor <- 1 / 1000 / 1000 # 'Time' in [µs[]
    battery.voltage <- batteryVoltageCustomCopter
  } else if (copter.model == "3DRSolo")  {
    data.starttime <- 0
    time.factor <- 1 / 1000 # 'Time' in [ms]
    battery.voltage <- batteryVoltageSoloCopter
  } else warning("unexpected model")

  # Select CURR
  temp <- data[data$V1 %in% c('CURR'), ]
  curr <- data.frame(
    type =     temp$V1,
    Time =     as.integer(temp$V2),
    TimeRelS = as.numeric((as.integer(temp$V2) - data.starttime) * time.factor),
    #Throttle = as.integer(temp$V3),
    Volt =     as.integer(temp$V4) / 100, # in [V]
    Curr =     as.integer(temp$V5) / 100, # in [A]
    Power =    as.numeric(as.integer(temp$V5)/100 * as.integer(temp$V4)/100 / battery.voltage * referenceVoltage), # in [W] normalized to 4s
    #Vcc =      as.integer(temp$V6),
    CurrTot =  as.numeric(temp$V7),       # in [mAh]
    PowerTot = as.numeric(as.numeric(temp$V7)/1000 * as.integer(temp$V4)/100 / battery.voltage * referenceVoltage), # in [Wh] normalized to 4s
    #Volt2 =    as.integer(temp$V8),
    stringsAsFactors = FALSE
  )

  # Select MODE
  temp <- data[data$V1 %in% c('MODE'), ]
  mode <- data.frame(
    type =     temp$V1,
    Time =     as.integer(temp$V2),
    TimeRelS = as.numeric((as.integer(temp$V2) - data.starttime) * time.factor),
    Mode =     as.character(temp$V3),
    ModeNum =  as.integer(temp$V4),
    stringsAsFactors = FALSE
  )

  # Select CMD
  temp <- data[data$V1 %in% c('CMD'), ]
  cmd <- data.frame(
    type =     temp$V1,
    Time =     as.integer(temp$V2),
    TimeRelS = as.numeric((as.integer(temp$V2) - data.starttime) * time.factor),
    CTot =     as.integer(temp$V3),
    CNum =     as.integer(temp$V4),
    CId =      as.integer(temp$V5),
    CName =    rapply(as.list(temp$V5), CId2type),
    Prm1 =     as.integer(temp$V6),
    Prm2 =     as.integer(temp$V7),
    Prm3 =     as.integer(temp$V8),
    Prm4 =     as.integer(temp$V9),
    Lat =      as.numeric(temp$V10),
    Lng =      as.numeric(temp$V11),
    Alt =      as.numeric(temp$V12),
    stringsAsFactors = FALSE
  )
  # filter duplicate commands from "MSG [...] new mission" summary
  cmd <- cmd[cmd$TimeRelS >= 0, ]
  # merge mode and command logs
  mode.cmd <- merge(mode, cmd, by = c("Time", "TimeRelS", "type"), sort = TRUE, all = TRUE)

  # Select GPS
  temp <- data[data$V1 %in% c('GPS'), ]
  if (copter.model == "CustomCopter") {
    gps <- data.frame(
      Time =     as.integer(temp$V2),
      TimeRelS = as.numeric((as.integer(temp$V2) - data.starttime) * time.factor),
      Status =   as.integer(temp$V3),
      GMS =      as.integer(temp$V4),
      GWk =      as.integer(temp$V5),
      NSats =    as.integer(temp$V6),
      HDop =     as.numeric(temp$V7),
      Lat =      as.numeric(temp$V8),
      Lng =      as.numeric(temp$V9),
      RelAlt =   as.numeric(temp$V10),
      Alt =      as.numeric(temp$V11),
      Speed =    as.numeric(temp$V12),
      GCrs =     as.numeric(temp$V13),
      VZ =       as.numeric(temp$V14),
      U =        as.integer(temp$V15),
      stringsAsFactors = FALSE
    )
  } else if (copter.model == "3DRSolo") {
    gps <- data.frame(
      type =     temp$V1,
      Status =   as.integer(temp$V2),
      TimeGPS =  as.integer(temp$V3),
      Time =     as.integer(temp$V14),
      TimeRelS = as.numeric((as.integer(temp$V14) - data.starttime) * time.factor),
      Week =     as.integer(temp$V4),
      NSats =    as.integer(temp$V5),
      HDop =     as.numeric(temp$V6),
      Lat =      as.numeric(temp$V7),
      Lng =      as.numeric(temp$V8),
      RelAlt =   as.numeric(temp$V9),
      Alt =      as.numeric(temp$V10),
      Speed =    as.numeric(temp$V11),
      GCrs =     as.numeric(temp$V12),
      VZ =       as.numeric(temp$V13),
      stringsAsFactors = FALSE
    )
  }

  return(list("copter.model" = copter.model, "curr" = curr, "mode" = mode, "cmd" = cmd, "mode.cmd" = mode.cmd, "gps" = gps))
}

#####################################################################
# Function: Cut into seperate commands

GenerateSections <- function(mode.cmd, class, descr) {
  sections <- data.frame()
  in.auto.mode = FALSE
  for (i in 1:nrow(mode.cmd)) {
    if (mode.cmd[i, "type"] == "MODE") {
      in.auto.mode = if(mode.cmd[i, "Mode"] == "Auto") TRUE else FALSE
      if (exists("temp.section")) remove(temp.section)
    }
    if ((mode.cmd[i, "type"] == "CMD") & in.auto.mode) {
      if (! exists("temp.section")) {
        temp.section <- data.frame()
      } else {
        if (mode.cmd[i, "CName"] == "LOITER_TIME")
          temp.section[1, "timestamp_end"] <- mode.cmd[i, "TimeRelS"] - mode.cmd[i, "Prm1"]
        else
          temp.section[1, "timestamp_end"] <- mode.cmd[i, "TimeRelS"]
        sections <- rbind(sections, temp.section)
      }
      temp.section[1, "flight"] <- descr
      temp.section[1, "class"] <- class
      temp.section[1, "cmd_name"] <- mode.cmd[i,"CName"]
      temp.section[1, "timestamp_start"] <- mode.cmd[i,"TimeRelS"]
      temp.section[1, "cmd_id_prev"] <- if (exists("cmd.id.prev")) cmd.id.prev else NA
      temp.section[1, "cmd_id_this"] <- mode.cmd[i,"CNum"]
      cmd.id.prev <- mode.cmd[i,"CNum"]
    }
  }
  #remove(i, in.auto.mode, temp.section, cmd.id.prev)
  return(sections)
}

#####################################################################
# Function: Calculate movement distance

AddMovementData <- function(sections, mode.cmd, gps) {
  for (i in 1:nrow(sections)) {
    if(is.na(sections[i,"cmd_id_prev"]) | sections[i,"cmd_name"] == "RETURN_TO_LAUNCH") {
      sections[i, "cmd_angle"] <- NA
      sections[i, "gps_angle"] <- NA
      sections[i, "cmd_distance"] <- NA
      sections[i, "gps_distance"] <- NA
      sections[i, "gps_speed"] <- NA
      sections[i, "interpretation"] <- "ung. Daten"
      next()
    }

    # based on programed commands
    cmd.location.start <- mode.cmd[mode.cmd$CNum %in% sections[i, "cmd_id_prev"], c("TimeRelS", "Lat", "Lng", "Alt")]
    cmd.location.end <- mode.cmd[mode.cmd$CNum %in% sections[i, "cmd_id_this"], c("TimeRelS", "Lat", "Lng", "Alt")]
    cmd.time <- cmd.location.end$TimeRelS - cmd.location.start$TimeRelS
    cmd.distance.on.Ground <- distGeo(c(cmd.location.start$Lng, cmd.location.start$Lat), c(cmd.location.end$Lng, cmd.location.end$Lat))
    cmd.altitude.diff <- cmd.location.end$Alt - cmd.location.start$Alt
    cmd.angle <- atan2(cmd.altitude.diff, cmd.distance.on.Ground)
    cmd.angle <- Rad2deg(cmd.angle)
    cmd.angle <- round(cmd.angle, digits = 1)
    cmd.distance <- sqrt(cmd.distance.on.Ground ^ 2 + abs(cmd.altitude.diff) ^ 2)

    # based on GPS data
    gps.location.start <- gps[which.min(abs(gps$TimeRelS - sections[i, "timestamp_start"])), ]
    gps.location.end <- gps[which.min(abs(gps$TimeRelS - sections[i, "timestamp_end"])), ]
    gps.time <- gps.location.end$TimeRelS - gps.location.start$TimeRelS
    gps.distance.on.Ground <- distGeo(c(gps.location.start$Lng, gps.location.start$Lat), c(gps.location.end$Lng, gps.location.end$Lat))
    gps.altitude.diff <- gps.location.end$RelAlt - gps.location.start$RelAlt
    gps.angle <- atan2(gps.altitude.diff, gps.distance.on.Ground)
    gps.angle <- Rad2deg(gps.angle)
    gps.angle <- round(gps.angle, digits = 1)
    gps.distance <- sqrt(gps.distance.on.Ground ^ 2 + abs(gps.altitude.diff) ^ 2)
    gps.speed <- gps.distance / gps.time

    sections[i, "cmd_angle"] <- cmd.angle         # in [°]
    sections[i, "gps_angle"] <- gps.angle         # in [°]
    sections[i, "cmd_distance"] <- cmd.distance   # in [m]
    sections[i, "gps_distance"] <- gps.distance   # in [m]
    sections[i, "gps_speed"] <- gps.speed         # in [m/s]


    sections[i, "interpretation"] <-
      if(cmd.angle == 0) {
        if(cmd.distance == 0) "Hover" else "Straight Forward"
      }
      else if(cmd.angle == +90) "Straight Up"
      else if(cmd.angle == -90) "Straight Down"
      else if(cmd.angle > 0) "Rise at Angle"
      else if(cmd.angle < 0) "Decline at Angle"
      else warning()
  }
  return(sections)
}

#####################################################################
# Function: Calculate statistical consumption values

AddConsumptionData <- function(sections, curr) {
  for (i in 1:nrow(sections)) {
    if(is.na(sections[i,"cmd_id_prev"]) | sections[i,"cmd_name"] == "RETURN_TO_LAUNCH") {
      sections[i, "n_rows"] <- NA
      sections[i, "mean"] <- NA
      sections[i, "stddeviation"] <- NA
      sections[i, "energy"] <- NA
      next()
    }
    # full sample data
    delayed.start = sections[i, "timestamp_start"]
    delayed.end = sections[i, "timestamp_end"]
    temp.curr <- curr[(curr$TimeRelS >= delayed.start & curr$TimeRelS <= delayed.end), ]
    if (nrow(temp.curr) == 0) warning("temp.curr empty.")
    energy <- tail(temp.curr, 1)["CurrTot"]  - head(temp.curr, 1)["CurrTot"]

    # sample data without transition
    delayed.start = sections[i, "timestamp_start"] + transitionSampleSeconds
    delayed.end = sections[i, "timestamp_end"] - transitionSampleSeconds
    temp.curr <- curr[(curr$TimeRelS >= delayed.start & curr$TimeRelS <= delayed.end), ]
    sections[i, "n_rows"] <- nrow(temp.curr)
    if(nrow(temp.curr) <= minSamplesForMean) {
      sections[i, "mean"] <- NA
      sections[i, "stddeviation"] <- NA
      sections[i, "energy"] <- NA
      next()
    }
    times = c()
    samples = c()
    for (j in 1:nrow(temp.curr)) {
      if (j == 1) next()
      step.time <- temp.curr[j, "TimeRelS"] - temp.curr[j-1, "TimeRels"]
      step.sample <- temp.curr[j-1, "Power"]
      #print(paste("Time", step.time, "mean", step.mean))
      times <- append(times, step.time)
      samples <- append(samples, step.sample)
    }

    # resulting data
    sections[i, "mean"] <- wtd.mean(samples, weights=times, na.rm=FALSE) # in [W]
    sections[i, "stddeviation"] <- sqrt(wtd.var(samples, weights=times, na.rm=FALSE)) # in [W]
    sections[i, "energy"] <- energy # in [mAh]

    #print(wtd.quantile(samples, weights=times, probs=0.75, na.rm=FALSE))

    #plot_ly(temp.curr, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]") %>% add_lines(y=sections[i,"mean"], name="mean")
    #plot_ly(x = samples, type = "histogram")
  }
  return(sections)
}

#####################################################################
# Function: Calculate combined data
# https://stats.stackexchange.com/questions/25848/how-to-sum-a-standard-deviation

GetCombinedAngleData <- function(sections) {
  combined <- data.frame()
  angles <- levels(factor(sections$cmd_angle))
  for (i in 1:length(angles)) {
    angle.sections <- sections[sections$cmd_angle == angles[i], ]
    count <- nrow(angle.sections)
    speed.mean <- sum(angle.sections$gps_speed) / count
    mean.mean <- sum(angle.sections$mean) / count
    sum.stddeviation <- sqrt(sum(angle.sections$stddeviation ^2))
    combined[i, "cmd_angle"] <- angles[i]
    combined[i, "sect_count"] <- count
    combined[i, "speed_mean"] <- speed.mean
    combined[i, "power_mean"] <- mean.mean
    combined[i, "power_stddev"] <- sum.stddeviation

    # energy_pm in [mAh/m]
    combined[i, "energy_pm_mean"] <- wtd.mean(angle.sections$energy / angle.sections$gps_distance, weights=angle.sections$gps_distance, na.rm=FALSE)
    combined[i, "energy_pm_stddev"] <- sqrt(wtd.var(angle.sections$energy / angle.sections$gps_distance, weights=angle.sections$gps_distance, na.rm=FALSE))
  }
  return(combined)
}


#####################################################################
#    HOVER - Static Analysis of one specific flight.
#####################################################################

if(exists("hoverlog")) {
  print("Hoverlog ... ")
  filename <- hoverlog
  timestamp.start <- 195.776
  timestamp.end <- 723.918

  cat("Load File ... ")
  logdata = read.csv(file = filename, sep = ",", header = FALSE, strip.white = TRUE, stringsAsFactors = FALSE)

  cat("Parse File ... ")
  .data <- ParseLogdata(logdata)
  logdata.curr <- .data$curr
  logdata.mode <- .data$mode
  logdata.cmd <- .data$cmd
  logdata.mode.cmd <- .data$mode.cmd
  logdata.gps <- .data$gps

  cat("Generate Sections (fixed timestamps) ... ")
  sections.file <- data.frame(
    flight = "Hover 10min",
    class = 1,
    cmd_name = "LOITER_UNLIM",
    timestamp_start = timestamp.start,
    cmd_id_prev = 1,
    cmd_id_this = 2,
    timestamp_end = timestamp.end,
    stringsAsFactors = FALSE
  )

  cat("AddConsumption ... ")
  sections.file <- AddConsumptionData(sections.file, logdata.curr)
  cat("\n")

  logdata.curr.hover <- logdata.curr[(logdata.curr$TimeRelS >= timestamp.start & logdata.curr$TimeRelS <= timestamp.end), ]

  p <- plot_ly(logdata.curr, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]")
  p <- add_markers(p, x = logdata.mode$TimeRelS, y = 5, text = paste('MODE', logdata.mode$Mode, sep=' '), name = "modes")
  p <- add_markers(p, x = logdata.cmd$TimeRelS, y = logdata.cmd$CId / 5, text = paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
  p <- add_trace(p, type="bar", x = logdata.cmd$TimeRelS, y = 25, opacity = 0.1, text=paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
  p
  plotly_POST(p, filename = "Hover")
  p2 <- plot_ly(logdata.curr.hover, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]")
  p2
  plotly_POST(p2, filename = "HoverSteady")
  p3 <- plot_ly(logdata.curr.hover, x = ~Power, type = "histogram")
  p3
  plotly_POST(p3, filename = "HoverSteadyHistogram")
  hover.mean <- mean(logdata.curr.hover$Power)
  hover.stddeviation <- sqrt(var(logdata.curr.hover$Power))
  
  hover.curr.mean1 <- mean(logdata.curr.hover$Curr)
  hover.curr.mean2 <- 60*60 / 1000 * (tail(logdata.curr.hover,1)[,"CurrTot"] - head(logdata.curr.hover,1)[,"CurrTot"]) / ((tail(logdata.curr.hover,1)[,"TimeRelS"] - head(logdata.curr.hover,1)[,"TimeRelS"]))

  # H O V E R   R E S U L T S
  # TimeRelS_start = 195.776
  # TimeRelS_end = 723.918
  # Time Diff = 528,142 sec (= 8min 48 sec)
  # power mean = 262.7 W
  # power stddev = 5.2 W
  # current mean = 18.09 A
  # current consumption total 3365-683 = 2682 mAh
  # current consumption total by seconds 2682/528 * (60*60/1000) = 18.29 mA
  # power consumption total 47.71570-10.25866 = 37.45704 Wh

  remove(filename, timestamp.start, timestamp.end)
  remove(hover.curr.mean1, hover.curr.mean2)
  remove(p, p2, p3)
}


#####################################################################
#    MAIN
#####################################################################

sections.all <- data.frame()

for (i in 1:nrow(logfiles)) {
  filename = logfiles[i, "filename"]
  description = logfiles[i, "description"]
  class = as.integer(logfiles[i, "class"])

  if(i == 1) cat("\n\n")

  cat(paste("Analyzing", filename, "- "))

  cat("Load File ... ")
  logdata = read.csv(file = filename, sep = ",", header = FALSE, strip.white = TRUE, stringsAsFactors = FALSE)

  cat("Parse File ... ")
  .data <- ParseLogdata(logdata)
  logdata.curr <- .data$curr
  logdata.mode <- .data$mode
  logdata.cmd <- .data$cmd
  logdata.mode.cmd <- .data$mode.cmd
  logdata.gps <- .data$gps

  cat("GenerateSections ... ")
  sections.file <- GenerateSections(logdata.mode.cmd, class, description)
  cat("AddMovement ... ")
  sections.file <- AddMovementData(sections.file, logdata.mode.cmd, logdata.gps)
  cat("AddConsumption ... ")
  sections.file <- AddConsumptionData(sections.file, logdata.curr)
  sections.all <- rbind(sections.all, sections.file)
  cat("\n")
}
remove(i, filename, description, class, sections.file)

sections.all.filtered <- sections.all[! is.na(sections.all$mean), ]

combined <- GetCombinedAngleData(sections.all.filtered)


#####################################################################
# Plot first test graph

x <- list(title = "Time [s]")
y <- list(title = "Power [W]")
p <- plot_ly(logdata.curr, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]") %>% layout(xaxis = x, yaxis = y)
#p <- add_markers(p, x = logdata.mode$TimeRelS, y = 5, text = paste('MODE', logdata.mode$Mode, sep=' '), name = "modes")
#p <- add_markers(p, x = logdata.cmd$TimeRelS, y = logdata.cmd$CId / 5, text = paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
#p <- add_trace(p, type="bar", x = logdata.cmd$TimeRelS, y = 420, opacity = 0.1, text=paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
p
#remove(p)
plotly_POST(p, filename = "LastFlightPower")

# quick visual inspection
plot(sections.all.filtered$mean, x = sections.all.filtered$cmd_angle)
plot(combined$power_mean, x = combined$cmd_angle)

plot(sections.all.filtered$gps_speed, x = sections.all.filtered$cmd_angle)
plot(combined$speed_mean, x = combined$cmd_angle)

#####################################################################
# Graphing

size.factor = 5

sizes <- c(1 * min(sections.all.filtered$stddeviation), 1 * max(sections.all.filtered$stddeviation))
#sizes <- c(size.factor * 2, 2* max(sections.all.filtered$n_rows) / min(sections.all.filtered$n_rows))

colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
#class <- 'circle'
class <- sections.all.filtered$class

plotlyGraph <- plot_ly(sections.all.filtered,
    x = ~cmd_angle, y = sections.all.filtered$mean, color = ~flight, size = ~stddeviation, colors = colors,
    type = 'scatter', mode = 'markers', sizes = sizes, error_y = list(type = "data", array = ~stddeviation),
    marker = list(symbol = class, sizemode = 'diameter', opacity = 0.8, line = list(width = 2, color = '#FFFFFF')),
    text = ~paste(
      flight, '- CMD', cmd_id_this,
      '<br>', interpretation, '(', format(cmd_angle, digits = 3),
      '° )<br>Mean Power µ:', format(mean, digits = 4),
      'W<br>StdDeviation σ:', format(stddeviation, digits = 3),
      '<br>Speed:', format(gps_speed, digits = 3)
    )
) %>%
layout(title = 'Flight Analysis',
    xaxis = list(title = 'Climb Angle [Degree]',
        #gridcolor = 'rgb(255, 255, 255)',
        range = c(-120, +120),
        #type = 'log',
        #zerolinewidth = 1,
        #ticklen = 5,
        dtick = 22.5#,
        #gridwidth = 2
    ),
    yaxis = list(title = 'Power (normalized to 4s) [W]',
        #gridcolor = 'rgb(255, 255, 255)',
        range = c(0, 1.2 * max(sections.all.filtered$mean))#,
        #zerolinewidth = 1,
        #ticklen = 5,
        #gridwith = 2
    )#,
    #paper_bgcolor = 'rgb(243, 243, 243)',
    #plot_bgcolor = 'rgb(243, 243, 243)'
)
plotlyGraph

plotly_POST(plotlyGraph, filename = "energyProfileGraph")

remove(size.factor, sizes, colors, class)
