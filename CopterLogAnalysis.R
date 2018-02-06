#!/usr/bin/Rscript
#
# CopterLogAnalysis.R - Statistical multicopter energy consumption
# analysis based on Ardupilot dataflash log files.
# Copyright (C) 2017  Thomas Dietrich
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Clear workspace
remove(list = ls())

#####################################################################
# Logfiles to analyze ###############################################

# Solo Hover/LOITER_UNLIM
hoverlog <- "2016-11-06 16-41-01-hover.log"

#
logfiles = data.frame()

# CustomCopter
###logfiles <- rbind(logfiles, c("2016-02-28 15-56-34.log", "2016-02-28 15-56-34 CC 'Quadrat'", "", ""), stringsAsFactors = FALSE)
###logfiles <- rbind(logfiles, c("2016-02-28 17-27-54.log", "2016-02-28 17-27-54 CC '23' ersterTest", "", ""), stringsAsFactors = FALSE)
#
#logfiles <- rbind(logfiles, c("2016-03-11 11-51-26.log", "2016-03-11 11-51-26 CC '23'", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c("2016-03-11 12-01-04.log", "2016-03-11 12-01-04 CC '23'", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c("2016-03-11 12-09-31.log", "2016-03-11 12-09-31 CC '23' grAkku", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-03-11 15-24-14.log", "2016-03-11 15-24-14 CC '23' grGeschw", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-03-11 15-32-57.log", "2016-03-11 15-32-57 CC '23' klGeschw", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-03-11 15-46-32.log", "2016-03-11 15-46-32 CC '23' grGewicht", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c("2016-04-30 17-45-40.log", "2016-04-30 17-45-40 CC 'Steig'", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c("2016-04-30 18-03-59.log", "2016-04-30 18-03-59 CC 'Steig'", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c("2016-05-01 14-12-24.log", "2016-05-01 14-12-24 CC 'Steig' Wind", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c("2016-10-23 13-44-35.log", "2016-10-23 13-44-35 CC 'Lw-Quadrat'", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c("2016-10-23 13-57-28.log", "2016-10-23 13-57-28 CC 'Lw-Benchmark'", 0, ""), stringsAsFactors = FALSE)

# SoloCopter
logfiles <- rbind(logfiles, c("2016-10-23 14-12-01.log", "2016-10-23 14-12-01 Solo 'Lw-Quadrat'", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 14-32-45.log", "2016-10-23 14-32-45 Solo 'Lw-Benchmark' Akku a", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 14-41-58.log", "2016-10-23 14-41-58 Solo 'Lw-Benchmark' Akku b(Flug1)", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-23 14-49-50.log", "2016-10-23 14-49-50 Solo 'Lw-Benchmark' Akku b(Flug2)", 4, ""), stringsAsFactors = FALSE)

logfiles <- rbind(logfiles, c("2016-10-27 17-19-01.log", "2016-10-27 17-19-01 Solo 'Lw-Quadrat' 3/4 Akku", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-25-25.log", "2016-10-27 17-25-25 Solo 'Lw-Quadrat' voller Akku", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-29-25.log", "2016-10-27 17-29-25 Solo 'Lw-Quadrat' 80%", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-33-10.log", "2016-10-27 17-33-10 Solo 'Lw-Benchmark' 67%", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-40-31.log", "2016-10-27 17-40-31 Solo 'Lw-Benchmark2' voller Akku", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-48-39.log", "2016-10-27 17-48-39 Solo 'Lw-Benchmark' bei 57%", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c("2016-10-27 17-54-43.log", "2016-10-27 17-54-43 Solo 'Lw-Benchmark2' voller Akku", 4, ""), stringsAsFactors = FALSE)

names(logfiles) <- (c("filename", "description", "class", "weather"))

# class: https://plot.ly/r/reference/#scatter-marker

#####################################################################
# Settings ##########################################################

minSamplesForMean <- 30
transitionSampleSeconds <- 3

# Normalization
batteryVoltageCustomCopter   <- 3 * 3.7 # 11.1
batteryVoltageSoloCopter <- 4 * 3.7 # 14.8
referenceVoltage         <- 4 * 3.7

#####################################################################
# Install packages, load libraries ##################################

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

# plotlyUsername <- "user"
# plotlyApiKey <- "key"
source("plotlyCredentials.R")

Sys.setenv("plotly_username" = plotlyUsername)
Sys.setenv("plotly_api_key" = plotlyApiKey)
remove(plotlyUsername, plotlyApiKey)

#####################################################################
# Load functions files ##############################################
source("CopterLogAnalysisFunctions.R")

#####################################################################
#    HOVER - Static Analysis of one specific flight.   ##############
#####################################################################

if(exists("hoverlog")) {
  print("Hoverlog ... ")
  filename <- hoverlog
  timestamp.start <- 195.776
  timestamp.end <- 723.918

  cat("Load File ... ")
  logdata = read.csv(file = filename, sep = ",", header = FALSE, strip.white = TRUE, stringsAsFactors = FALSE)

  cat("Parse File ... ")
  raw.data <- ParseLogdata(logdata)
  logdata.curr <- raw.data$curr
  logdata.mode <- raw.data$mode
  logdata.cmd <- raw.data$cmd
  logdata.mode.cmd <- raw.data$mode.cmd
  logdata.gps <- raw.data$gps
  remove(raw.data)

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
  api_create(p, filename = "Hover")
  p2 <- plot_ly(logdata.curr.hover, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]")
  p2
  api_create(p2, filename = "HoverSteady")
  p3 <- plot_ly(logdata.curr.hover, x = ~Power, type = "histogram")
  p3
  api_create(p3, filename = "HoverSteadyHistogram")
  hover.mean <- mean(logdata.curr.hover$Power)
  hover.stddeviation <- sqrt(var(logdata.curr.hover$Power))
  
  hover.curr.mean <- mean(logdata.curr.hover$Curr)
  hover.curr.stddev <- sqrt(var(logdata.curr.hover$Curr))
  
  # for comparison and evaluation
  hover.curr.meanTot <- 60*60 / 1000 * (tail(logdata.curr.hover,1)[,"CurrTot"] - head(logdata.curr.hover,1)[,"CurrTot"]) / ((tail(logdata.curr.hover,1)[,"TimeRelS"] - head(logdata.curr.hover,1)[,"TimeRelS"]))

  # H O V E R   R E S U L T S
  # TimeRelS_start = 195.776
  # TimeRelS_end = 723.918
  # Time Diff = 528,142 sec (= 8min 48 sec)
  # power mean = 262.7 W
  # power stddev = 5.2 W
  # current mean = 18.09 A
  # current stddev = 0.36 A
  # current consumption total 3365-683 = 2682 mAh
  # current consumption total by seconds 2682/528 * (60*60/1000) = 18.29 A
  # power consumption total 47.71570-10.25866 = 37.45704 Wh

  #remove(filename, timestamp.start, timestamp.end)
  #remove(hover.curr.mean1, hover.curr.mean2)
  #remove(p, p2, p3)
}


#####################################################################
#    MAIN Analysis   ################################################
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
  raw.data <- ParseLogdata(logdata)
  logdata.curr <- raw.data$curr
  logdata.mode <- raw.data$mode
  logdata.cmd <- raw.data$cmd
  logdata.mode.cmd <- raw.data$mode.cmd
  logdata.gps <- raw.data$gps
  remove(raw.data)

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
# Plot first test graph #############################################

x <- list(title = "Time [s]")
y <- list(title = "Power [W]")
p <- plot_ly(logdata.curr, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]") %>% layout(xaxis = x, yaxis = y)
#p <- add_markers(p, x = logdata.mode$TimeRelS, y = 5, text = paste('MODE', logdata.mode$Mode, sep=' '), name = "modes")
#p <- add_markers(p, x = logdata.cmd$TimeRelS, y = logdata.cmd$CId / 5, text = paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
#p <- add_trace(p, type="bar", x = logdata.cmd$TimeRelS, y = 420, opacity = 0.1, text=paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
p
#remove(p)
api_create(p, filename = "LastFlightPower")

# quick visual inspection
plot(sections.all.filtered$mean, x = sections.all.filtered$cmd_angle)
plot(combined$power_mean, x = combined$cmd_angle)

plot(sections.all.filtered$gps_speed, x = sections.all.filtered$cmd_angle)
plot(combined$speed_mean, x = combined$cmd_angle)

#####################################################################
# Graphing ##########################################################

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
    autosize = F, width = 1500, height = 600,
    #autosize = F, width = 1000, height = 400,
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

api_create(plotlyGraph, filename = "energyProfileGraph")

remove(size.factor, sizes, colors, class)
