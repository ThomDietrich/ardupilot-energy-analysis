#!/usr/bin/Rscript

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

GenerateSections <- function(mode.cmd, class, model, descr) {
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
      temp.section[1, "model"] <- model
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
        if(cmd.distance == 0) "Hover" else "Straigth Forward"
      }
    else if(cmd.angle == +90) "Straigth Up"
    else if(cmd.angle == -90) "Straigth Down"
    else if(cmd.angle > 0) "Rise at Angle"
    else if(cmd.angle < 0) "Decline at Angle"
    else warning()
  }
  return(sections)
}

#####################################################################
# Function: Calculte statistical consumption values

AddConsumptionData <- function(sections, curr) {
  for (i in 1:nrow(sections)) {
    if(is.na(sections[i,"cmd_id_prev"]) | sections[i,"cmd_name"] == "RETURN_TO_LAUNCH") {
      sections[i, "n_rows"] <- NA
      sections[i, "power_mean"] <- NA
      sections[i, "power_stddeviation"] <- NA
      sections[i, "current"] <- NA
      sections[i, "current_stddev"] <- NA
      sections[i, "energy"] <- NA
      next()
    }
    # full sample data
    time.start = sections[i, "timestamp_start"]
    time.end = sections[i, "timestamp_end"]
    temp.curr <- curr[(curr$TimeRelS >= time.start & curr$TimeRelS <= time.end), ]
    if (nrow(temp.curr) == 0) warning("temp.curr empty.")
    energy.time <- time.end - time.start
    energy <- tail(temp.curr, 1)["CurrTot"]  - head(temp.curr, 1)["CurrTot"]
    remove(time.start, time.end)
    
    # sample data without transition
    delayed.start = sections[i, "timestamp_start"] + transitionSampleSeconds
    delayed.end = sections[i, "timestamp_end"] - transitionSampleSeconds
    temp.curr <- curr[(curr$TimeRelS >= delayed.start & curr$TimeRelS <= delayed.end), ]
    sections[i, "n_rows"] <- nrow(temp.curr)
    if(nrow(temp.curr) <= minSamplesForMean) {
      sections[i, "power_mean"] <- NA
      sections[i, "power_stddeviation"] <- NA
      sections[i, "current"] <- NA
      sections[i, "current_stddev"] <- NA
      sections[i, "energy"] <- NA
      sections[i, "energy_time"] <- NA
      sections[i, "energy_current"] <- NA
      next()
    }
    times = c()
    powers = c()
    currents = c()
    for (j in 1:nrow(temp.curr)) {
      if (j == 1) next()
      step.time <- temp.curr[j, "TimeRelS"] - temp.curr[j-1, "TimeRels"]
      step.power <- temp.curr[j-1, "Power"]
      step.current <- temp.curr[j-1, "Curr"]
      #print(paste("Time", step.time, "mean", step.mean))
      times <- append(times, step.time)
      powers <- append(powers, step.power)
      currents <- append(currents, step.current)
    }
    
    # resulting data
    sections[i, "power_mean"] <- wtd.mean(powers, weights=times, na.rm=FALSE)              # mean power over flight time, in [W]
    sections[i, "power_stddeviation"] <- sqrt(wtd.var(powers, weights=times, na.rm=FALSE)) # in [W]
    sections[i, "current"] <- wtd.mean(currents, weights=times, na.rm=FALSE)               # mean current, in [A]
    sections[i, "current_stddev"] <- sqrt(wtd.var(currents, weights=times, na.rm=FALSE))   # current stddev, in [A]
    sections[i, "energy"] <- energy                                                        # difference CurrTot before and after from full sample data, in [mAh]
    sections[i, "energy_time"] <- energy.time                                              # full sample data time, in [s]
    sections[i, "energy_current"] <- energy / energy.time * 60 * 60 / 1000                 # full sample data energy by time, in [A]
    
    #print(wtd.quantile(powers, weights=times, probs=0.75, na.rm=FALSE))
    
    #plot_ly(temp.curr, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]") %>% add_lines(y=sections[i,"mean"], name="mean")
    #plot_ly(x = powers, type = "histogram")
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
    mean.mean <- sum(angle.sections$power_mean) / count
    current.mean <- sum(angle.sections$current) / count
    current.stddev <- sum(angle.sections$current_stddev) / count
    sum.stddeviation <- sqrt(sum(angle.sections$power_stddeviation ^2))
    combined[i, "cmd_angle"] <- angles[i]            # angle for combination groups
    combined[i, "sect_count"] <- count               # number of sections combined for mean
    combined[i, "speed_mean"] <- speed.mean          # in [m/s]
    combined[i, "power_mean"] <- mean.mean           # in [W]
    combined[i, "power_stddev"] <- sum.stddeviation  # in [W]
    combined[i, "current_mean"] <- current.mean      # in [A]
    combined[i, "current_stddev"] <- current.stddev  # in [A]
    
    # energy used per meter, in [mAh/m]
    combined[i, "energy_pm_mean"] <- wtd.mean(angle.sections$energy / angle.sections$gps_distance, weights=angle.sections$gps_distance, na.rm=FALSE)
    combined[i, "energy_pm_stddev"] <- sqrt(wtd.var(angle.sections$energy / angle.sections$gps_distance, weights=angle.sections$gps_distance, na.rm=FALSE))
    
    # current calculated out of energy, in [A]
    combined[i, "energy_current_mean"] <- wtd.mean(angle.sections$energy_current, weights=angle.sections$gps_distance, na.rm=FALSE)
    combined[i, "energy_current_stddev"] <- sqrt(wtd.var(angle.sections$energy_current, weights=angle.sections$gps_distance, na.rm=FALSE))
  }
  return(combined)
}
