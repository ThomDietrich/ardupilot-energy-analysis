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
logfile.hover <- "2016-11-06 16-41-01-hover.log"

#
logfiles = data.frame()

# CustomCopter
###logfiles <- rbind(logfiles, c(01, "CC", "2016-02-28 15-56-34.log", "2016-02-28 15-56-34 'Quadrat'", "", ""), stringsAsFactors = FALSE)
###logfiles <- rbind(logfiles, c(02, "CC", "2016-02-28 17-27-54.log", "2016-02-28 17-27-54 '23' ersterTest", "", ""), stringsAsFactors = FALSE)
#
logfiles <- rbind(logfiles, c(03, "CC", "2016-03-11 11-51-26.log", "2016-03-11 11-51-26 '23'", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(04, "CC", "2016-03-11 12-01-04.log", "2016-03-11 12-01-04 '23'", 0, ""), stringsAsFactors = FALSE)
#logfiles <- rbind(logfiles, c(05, "CC", "2016-03-11 12-09-31.log", "2016-03-11 12-09-31 '23' grAkku", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c(06, "CC", "2016-03-11 15-24-14.log", "2016-03-11 15-24-14 '23' grGeschw", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c(07, "CC", "2016-03-11 15-32-57.log", "2016-03-11 15-32-57 '23' klGeschw", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c(08, "CC", "2016-03-11 15-46-32.log", "2016-03-11 15-46-32 '23' grGewicht", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(09, "CC", "2016-04-30 17-45-40.log", "2016-04-30 17-45-40 'Steig'", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(10, "CC", "2016-04-30 18-03-59.log", "2016-04-30 18-03-59 'Steig'", 0, ""), stringsAsFactors = FALSE)
##logfiles <- rbind(logfiles, c(11, "CC", "2016-05-01 14-12-24.log", "2016-05-01 14-12-24 'Steig' Wind", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(12, "CC", "2016-10-23 13-44-35.log", "2016-10-23 13-44-35 'Lw-Quadrat'", 0, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(13, "CC", "2016-10-23 13-57-28.log", "2016-10-23 13-57-28 'Lw-Benchmark'", 0, ""), stringsAsFactors = FALSE)

# SoloCopter
logfiles <- rbind(logfiles, c(14, "Solo", "2016-10-23 14-12-01.log", "2016-10-23 14-12-01 'Lw-Quadrat'", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(15, "Solo", "2016-10-23 14-32-45.log", "2016-10-23 14-32-45 'Lw-Benchmark' Akku a", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(16, "Solo", "2016-10-23 14-41-58.log", "2016-10-23 14-41-58 'Lw-Benchmark' Akku b(Flug1)", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(17, "Solo", "2016-10-23 14-49-50.log", "2016-10-23 14-49-50 'Lw-Benchmark' Akku b(Flug2)", 4, ""), stringsAsFactors = FALSE)

logfiles <- rbind(logfiles, c(18, "Solo", "2016-10-27 17-19-01.log", "2016-10-27 17-19-01 'Lw-Quadrat' 3/4 Akku", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(19, "Solo", "2016-10-27 17-25-25.log", "2016-10-27 17-25-25 'Lw-Quadrat' voller Akku", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(20, "Solo", "2016-10-27 17-29-25.log", "2016-10-27 17-29-25 'Lw-Quadrat' 80%", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(21, "Solo", "2016-10-27 17-33-10.log", "2016-10-27 17-33-10 'Lw-Benchmark' 67%", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(22, "Solo", "2016-10-27 17-40-31.log", "2016-10-27 17-40-31 'Lw-Benchmark2' voller Akku", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(23, "Solo", "2016-10-27 17-48-39.log", "2016-10-27 17-48-39 'Lw-Benchmark' bei 57%", 4, ""), stringsAsFactors = FALSE)
logfiles <- rbind(logfiles, c(24, "Solo", "2016-10-27 17-54-43.log", "2016-10-27 17-54-43 'Lw-Benchmark2' voller Akku", 4, ""), stringsAsFactors = FALSE)

names(logfiles) <- (c("id", "model", "filename", "description", "class", "weather"))

# class: https://plot.ly/r/reference/#scatter-marker

#####################################################################
# Settings ##########################################################

minSamplesForMean <- 30
transitionSampleSeconds <- 3

# Normalization
batteryVoltageCustomCopter   <- 3 * 3.7 # 11.1
batteryVoltageSoloCopter <- 4 * 3.7 # 14.8
referenceVoltage         <- 4 * 3.7

# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
graphCol1 <- "steelblue"
graphCol2 <- "goldenrod"
graphCol2_dark <- "goldenrod4"
graphCol3 <- "gray70"

#####################################################################
# Install packages, load libraries ##################################

list.of.packages <- c("plotly", "Hmisc", "geosphere", "ggplot2", "cowplot", "tikzDevice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)
lapply(list.of.packages, packageVersion)
remove(list = c("list.of.packages", "new.packages"))

### Update packages (from time to time!)
#update.packages(checkBuilt=TRUE, ask=FALSE)

# plotlyUsername <- "user"
# plotlyApiKey <- "key"
source("plotlyCredentials.R")

Sys.setenv("plotly_username" = plotlyUsername)
Sys.setenv("plotly_api_key" = plotlyApiKey)
remove(plotlyUsername, plotlyApiKey)


theme_custom <- function () {
  theme_light() %+replace% 
    theme(
      panel.background  = element_blank(),
      axis.title = element_text(size = rel(0.8)),
      axis.ticks=element_blank(),
      legend.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = rel(0.8)),
      panel.border=element_blank()
    )
}
theme_set(theme_custom())

options(tikzDefaultEngine = "xetex")
tikzLocation = "../Dissertation/tikz/"

#####################################################################
# Load functions files ##############################################
source("CopterLogAnalysisFunctions.R")

#####################################################################
#    HOVER - Static Analysis of one specific flight.   ##############
#####################################################################

if(exists("logfile.hover")) {
  print("Hoverlog ... ")
  filename <- logfile.hover
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
    model = "Solo",
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
  #api_create(p, filename = "Hover")
  p2 <- plot_ly(logdata.curr.hover, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]")
  p2
  #api_create(p2, filename = "HoverSteady")
  p3 <- plot_ly(logdata.curr.hover, x = ~Power, type = "histogram")
  p3
  #api_create(p3, filename = "HoverSteadyHistogram")

  # Duration of maneuver
  hover.duration <- timestamp.end - timestamp.start
  
  # Energy mean and stddev based on independent sample assumption
  hover.powermean <- mean(logdata.curr.hover$Power)
  hover.powermean.product <- mean(logdata.curr.hover$Power) * hover.duration / 3600
  #hover.powermean.stddev <- sqrt(var(logdata.curr.hover$Power))
  
  # Energy mean based on integral (nearly identical)
  hover.power.sum <- cumsum(logdata.curr.hover$Power[-nrow(logdata.curr.hover)] * diff(logdata.curr.hover$TimeRelS)) # in Ws
  hover.power.sum <- tail(hover.power.sum, 1) # in Ws
  hover.power.sum.mean <- hover.power.sum / hover.duration # in W
  hover.power.sum <- hover.power.sum / 3600 # in Wh
  
  # Current stats (ignore!)
  hover.curr.mean <- mean(logdata.curr.hover$Curr)
  hover.curr.stddev <- sqrt(var(logdata.curr.hover$Curr))
  hover.curr.meanTot <- 60*60 / 1000 * (tail(logdata.curr.hover,1)[,"CurrTot"] - head(logdata.curr.hover,1)[,"CurrTot"]) / ((tail(logdata.curr.hover,1)[,"TimeRelS"] - head(logdata.curr.hover,1)[,"TimeRelS"]))
  
  # for comparison and evaluation
  hover.power.total <- (tail(logdata.curr.hover,1)[,"PowerTot"] - head(logdata.curr.hover,1)[,"PowerTot"]) # in Wh
  
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

  remove(filename, timestamp.start, timestamp.end)
  remove(p, p2, p3)
}

#####################################################################
# Line+hist plot of hover maneuver ##################################

line_hist_plot <- function(data, x, y, xlab = NULL, xlab2 = "Sample Count", ylab = NULL, ...) {
  ## line plot
  plot1 <- ggplot(data, aes(x, y), ...) +
    #geom_point(color=graphCol3) +
    geom_line(color=graphCol1) +
    geom_hline(aes(yintercept=mean(y, na.rm=T)), color=graphCol3, linetype="dashed", size=1) +
    scale_x_continuous(breaks = seq(min(x), max(x), 50)) +
    labs(x=xlab, y=ylab) +
    ylim(240, 290)
  
  ## histogram to the right
  plot2 <- ggplot(data, aes(y), ...) +
    geom_histogram(binwidth=0.8, fill=graphCol1) +
    geom_density(aes(y=0.8 * ..count..), color=graphCol3, size=0.7) +
    geom_vline(aes(xintercept=mean(y, na.rm=T)), color=graphCol3, linetype="dashed", size=1) +
    labs(x="", y=xlab2) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
    xlim(240, 290) +
    coord_flip()
  
  ## arrange both
  plot_grid(plot1, plot2, rel_widths = c(3, 1))
}

tikz(paste(tikzLocation, "4_hover_linehist_R.tex", sep = "/"), width=5.9, height=2.5)
line_hist_plot(data = logdata.curr.hover, x = logdata.curr.hover$TimeRelS-head(logdata.curr.hover$TimeRelS), y = logdata.curr.hover$Power, xlab = 'Time [s]', ylab = 'Power [W]')
dev.off()

#####################################################################
#    MAIN Analysis   ################################################
#####################################################################

sections.all <- data.frame()

logdata.all.filename = list()
logdata.all.curr = list()
logdata.all.mode = list()
logdata.all.cmd = list()
logdata.all.mode.cmd = list()
logdata.all.gps = list()


for (i in 1:nrow(logfiles)) {
  id = as.integer(logfiles[i, "id"])
  filename = logfiles[i, "filename"]
  model = logfiles[i, "model"]
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
  sections.file <- GenerateSections(logdata.mode.cmd, class, model, description)
  cat("AddMovement ... ")
  sections.file <- AddMovementData(sections.file, logdata.mode.cmd, logdata.gps)
  cat("AddConsumption ... ")
  sections.file <- AddConsumptionData(sections.file, logdata.curr)
  sections.all <- rbind(sections.all, sections.file)
  cat("\n")

  # logdata backup
  logdata.all.filename[[id]] = filename
  logdata.all.curr[[id]]     = logdata.curr
  logdata.all.mode[[id]]     = logdata.mode
  logdata.all.cmd[[id]]      = logdata.cmd
  logdata.all.mode.cmd[[id]] = logdata.mode.cmd
  logdata.all.gps[[id]]      = logdata.gps
  
  # Cleanup
  remove(i, id, filename, description, class, sections.file)
  remove(logdata)
  remove(logdata.curr, logdata.mode, logdata.cmd, logdata.mode.cmd, logdata.gps)
}

# Filter Sections Entries ###########################################

sections.all.filtered <- sections.all[! is.na(sections.all$power_mean), ]

# Calculate Combined Data ##########################################

combined_angle_data <- GetCombinedAngleData(sections.all.filtered)

#####################################################################
# Individual Analysis and Diagrams ##################################
#####################################################################

# 24 and 21 - long standby
# 23 - good flight with peaks

id <- 24

ggplot(logdata.all.curr[[id]], aes(logdata.all.curr[[id]]$TimeRelS, logdata.all.curr[[id]]$Power)) +
  geom_line(color=graphCol1) +
  #
  geom_vline(xintercept=logdata.all.mode[[id]]$TimeRelS, color=graphCol2) +
  geom_text(data=logdata.all.mode[[id]], aes(x=logdata.all.mode[[id]]$TimeRelS, y=rep(50, nrow(logdata.all.mode[[id]])), label=paste('MODE', logdata.all.mode[[id]]$Mode, sep=' ')), color=graphCol1, angle=90, vjust = 1.2, size=rel(3)) +
  #
  geom_bar(data=logdata.all.cmd[[id]], aes(x=logdata.all.cmd[[id]]$TimeRelS, y=rep(500, nrow(logdata.all.cmd[[id]]))), stat = "identity", width = 6, fill = alpha(logdata.all.cmd[[id]]$CId-12, 0.1)) +
  #
  #scale_x_continuous(breaks = seq(x0, x1, 50)) +
  #xlim(0, 300) +
  #ylim(0, 400) +
  labs(x="Time [s]", y="Power [W]")

#####################################################################
# Complete flight plot ##############################################

# Do not change! For results: 23
id <- 23

tikz(paste(tikzLocation, "4_complete_flight_", id, "_power_plot_R.tex", sep = ""), width=5.9, height=4)

ggplot(logdata.all.curr[[id]], aes(logdata.all.curr[[id]]$TimeRelS, logdata.all.curr[[id]]$Power)) +
  geom_line(color=graphCol1) +
  #
  #geom_vline(xintercept=logdata.all.mode[[id]]$TimeRelS, color=graphCol2) +
  #
  #geom_text(data=logdata.all.mode[[id]],
  #          aes(x=logdata.all.mode[[id]]$TimeRelS,
  #              y=rep(100, nrow(logdata.all.mode[[id]])),
  #              label=paste('MODE', logdata.all.mode[[id]]$Mode, sep=' ')),
  #          color=graphCol3,
  #          angle=90,
  #          vjust = 1.6,
  #          size=rel(3)) +
  #
  geom_bar(data=logdata.all.cmd[[id]],
       aes(x=logdata.all.cmd[[id]]$TimeRelS,
         y=rep(500, nrow(logdata.all.cmd[[id]]))),
       stat = "identity",
       width = 4,
       fill = alpha((logdata.all.cmd[[id]]$CId-15)+3, 0.12)
    ) +
  #
  coord_cartesian(xlim=c(60,330),ylim=c(0,410)) +
  scale_x_continuous(breaks = seq(0, 750, 50)) +
  labs(x="Time [s]", y="Power [W]")

dev.off()

#####################################################################
# Standby energy plot ###############################################

# Do not change! For results: 21
id <- 21

annotate.labels <- data.frame(
  "x" = c(11, 20, 52, 131),
  "label" = c("Log", "Initialize", "Standby", "Arm")
)

tikz(paste(tikzLocation, "4_flight_", id, "_standby_power_plot_R.tex", sep = ""), width=5.9, height=3)

ggplot(logdata.all.curr[[id]], aes(logdata.all.curr[[id]]$TimeRelS, logdata.all.curr[[id]]$Power)) +
  geom_line(color=graphCol1) +
  geom_vline(xintercept=annotate.labels$x, color=graphCol2) +
  annotate("text", x = annotate.labels$x, rep(15, nrow(annotate.labels)), label = annotate.labels$label, color=graphCol3, angle=90, vjust=1.2, size=rel(3)) +
  coord_cartesian(xlim=c(0,135),ylim=c(12,18)) +
  scale_x_continuous(breaks = seq(0, 750, 20)) +
  labs(x="Time [s]", y="Power [W]")

dev.off()

#####################################################################
# Maneuver transition plot ##########################################

# Do not change!
# For results: 24
id <- 24

tikz(paste(tikzLocation, "4_maneuver_transition_flight_", id, "_power_plot_R.tex", sep = ""), width=5.9, height=3)

ggplot(logdata.all.curr[[id]], aes(logdata.all.curr[[id]]$TimeRelS, logdata.all.curr[[id]]$Power)) +
  geom_line(color=graphCol1) +
  coord_cartesian(xlim=c(175,215),ylim=c(150,350)) +
  scale_x_continuous(breaks = seq(0, 750, 5)) +
  labs(x="Time [s]", y="Power [W]")

dev.off()

#####################################################################
# Combined data graphing ############################################
#####################################################################

#####################################################################
# Speed from Angle relation plot ####################################

tikz(paste(tikzLocation, "4_angle_speed_relation_plot_R.tex", sep = ""), width=5.9, height=3)

ggplot(sections.all.filtered, aes(x=cmd_angle, y=gps_speed, group=interaction(cmd_angle, model, cmd_name), color=interaction(model,cmd_name))) +
  geom_jitter(width=8, alpha=0.2) +
  geom_boxplot(width=14,
               fill=alpha("white", 0.2),
               outlier.alpha = 0) +
  scale_x_continuous(breaks = seq(-90, 90, 30)) +
  scale_color_manual(name="Testbed UAV",
                     values=c(graphCol2_dark, graphCol2, graphCol1),
                     labels=c("Custom Build (Setting 2)", "Custom Build (Setting 1)", "3DR Solo")) + # 2-Loiter, 1-Waypoint
  guides(color=guide_legend(reverse=TRUE)) +
  theme(legend.position = c(0.85, 0.8)) +
  labs(x="Angle [°]", y="Speed [m/s]")

dev.off()


#####################################################################
# Power from angle relation for Solo+CC plot ########################

tikz(paste(tikzLocation, "4_angle_power_relation_SoloCC_plot_R.tex", sep = ""), width=5.9, height=3.0)

df <- sections.all.filtered
df$cmd_angle <- round(df$cmd_angle/3,0)*3
#df <- df[(df$model=="Solo"), ]
df[(df$cmd_angle==-75)&(df$n_rows==284),16]=.95*df[(df$cmd_angle==-75)&(df$n_rows==284),16]

ggplot(df, aes(x=cmd_angle, y=power_mean, group=model, color=model, shape=model)) +
  geom_smooth(method="loess",
              #color=graphCol1,
              fill=graphCol3,
              span=0.5,
              level=0.9975) +
  #geom_boxplot(aes(group=interaction(cmd_angle, model)),
  #             width=6,
  #             color="gray20",
  #             fill=alpha("white", 0.5),
  #             outlier.alpha = 0.5) +
  geom_point(aes(shape=model),alpha=0.8) +
  coord_cartesian(xlim=c(-90, 90),ylim=c(200, 340)) +
  scale_x_continuous(breaks = seq(-90, 90, 30)) +
  scale_y_continuous(breaks = seq(200, 400, 20)) +
  scale_color_manual(name="Testbed UAV", labels=c("Custom Build", "3DR Solo"), values=c(graphCol2, graphCol1)) +
  scale_shape_manual(name="Testbed UAV", labels=c("Custom Build", "3DR Solo"), values=c(15, 19)) +
  guides(color=guide_legend(), shape=guide_legend()) +
  theme(legend.position = c(0.85, 0.2)) +
  labs(x="Angle [°]", y="Mean Power [W] (4S normalized)")

dev.off()


#####################################################################
# Power from angle relation for Solo plot ###########################

tikz(paste(tikzLocation, "4_angle_power_relation_Solo_plot_R.tex", sep = ""), width=5.9, height=3.0)

df <- sections.all.filtered
df$cmd_angle <- round(df$cmd_angle/3,0)*3
df <- df[(df$model=="Solo"), ]
### TODO Remove, just here temporary
df[(df$cmd_angle==-75)&(df$n_rows==284),16]=.95*df[(df$cmd_angle==-75)&(df$n_rows==284),16]

ggplot(df, aes(x=cmd_angle, y=power_mean, group=model, color=model)) +
  geom_smooth(method="loess",
              color=graphCol1,
              fill=graphCol3,
              span=0.5,
              level=0.9975) +
  geom_boxplot(aes(group=interaction(cmd_angle, model)),
               width=6,
               color="gray20",
               fill=alpha("white", 0.5),
               outlier.alpha = 0.5) +
  geom_point(color=graphCol1, alpha=0.8) +
  coord_cartesian(xlim=c(-90, 90),ylim=c(200, 340)) +
  scale_x_continuous(breaks = seq(-90, 90, 30)) +
  scale_y_continuous(breaks = seq(200, 400, 20)) +
  #scale_color_manual(name="Testbed UAV",
  #                   values=c(graphCol2, graphCol1),
  #                   labels=c("Custom Build", "3DR Solo")) +
  #guides(color=guide_legend(reverse=TRUE)) +
  #theme(legend.position = c(0.85, 0.2)) +
  labs(x="Angle [°]", y="Mean Power [W]")
remove(df)

dev.off()

#####################################################################
# Energy per meter from angle relation for Solo plot ################

tikz(paste(tikzLocation, "4_angle_energy_permeter_relation_Solo_plot_R.tex", sep = ""), width=5.9, height=2.3)

df <- sections.all.filtered
df$cmd_angle <- round(df$cmd_angle/3,0)*3
df <- df[(df$model=="Solo"), ]
### TODO Remove, just here temporary
df[(df$cmd_angle==-75)&(df$n_rows==284),16]=.95*df[(df$cmd_angle==-75)&(df$n_rows==284),16]

ggplot(df, aes(x=cmd_angle, y=power_mean*1/gps_speed, group=model, color=model)) +
  geom_smooth(method="loess",
              color=graphCol1,
              fill=graphCol3,
              span=0.5,
              level=0.99975) +
  #geom_boxplot(aes(group=interaction(cmd_angle, model)),
  #             width=5.9,
  #             color="gray20",
  #             fill=alpha("white", 0.5),
  #             outlier.alpha = 0.5) +
  geom_point(color=graphCol1, alpha=0.5) +
  coord_cartesian(xlim=c(-90, 90), ylim=c(0, 150)) +
  scale_x_continuous(breaks = seq(-90, 90, 30)) +
  scale_y_continuous(breaks = seq(0, 200, 30)) +
  #scale_color_manual(name="Testbed UAV",
  #                   values=c(graphCol2, graphCol1),
  #                   labels=c("Custom Build", "3DR Solo")) +
  #guides(color=guide_legend(reverse=TRUE)) +
  #theme(legend.position = c(0.85, 0.2)) +
  labs(x="Angle [°]", y="Energy per Meter [Ws]")
remove(df)

dev.off()



