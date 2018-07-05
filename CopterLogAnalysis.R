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
# Close all devices still open from previous executions
graphics.off()

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
#graphCol1 <- "steelblue"
#graphCol2 <- "goldenrod"
#graphCol2_dark <- "goldenrod4"
#graphCol3 <- "gray70"

#graphCol1 <- "#353594"
graphCol1 <- "#4682B4" #steelblue
graphCol2 <- "#B8850A" #mygold
graphCol2_dark <- "goldenrod4"
graphCol3 <- "#B3B3B3" #mygray
graphCol3_dark <- "#737373"
graphCol3_darkdark <- "#434343"
graphCol4 <- "#780116" #myred

#####################################################################
# Install packages, load libraries ##################################

library(data.table)
library(xtable)
list.of.packages <- c("Hmisc", "geosphere", "ggplot2", "cowplot", "tikzDevice", "TTR", "xts", "forecast")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)
lapply(list.of.packages, packageVersion)
remove(list = c("list.of.packages", "new.packages"))

# Increase max number of warnings
options(nwarnings=200) 

### Update packages (from time to time!)
#update.packages(checkBuilt=TRUE, ask=FALSE)

# plotlyUsername <- "user"
# plotlyApiKey <- "key"
#source("plotlyCredentials.R")

#Sys.setenv("plotly_username" = plotlyUsername)
#Sys.setenv("plotly_api_key" = plotlyApiKey)
#remove(plotlyUsername, plotlyApiKey)


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
    flight_id = 0,
    flight = "Hover 10min",
    flightcmd_id = 0,
    model = "Solo",
    class = 1,
    cmd_name = "LOITER_UNLIM",
    timestamp_start = timestamp.start,
    cmd_id_prev = 1,
    cmd_id_this = 2,
    gps_angle = NA,
    timestamp_end = timestamp.end,
    stringsAsFactors = FALSE
  )
  
  cat("AddConsumption ... ")
  ret_list <- AddConsumptionData(sections.file, logdata.curr)
  sections.file <- ret_list$sections
  cat("\n")
  
  logdata.curr.hover <- logdata.curr[(logdata.curr$TimeRelS >= timestamp.start & logdata.curr$TimeRelS <= timestamp.end), ]
  
  #p <- plot_ly(logdata.curr, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]")
  #p <- add_markers(p, x = logdata.mode$TimeRelS, y = 5, text = paste('MODE', logdata.mode$Mode, sep=' '), name = "modes")
  #p <- add_markers(p, x = logdata.cmd$TimeRelS, y = logdata.cmd$CId / 5, text = paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
  #p <- add_trace(p, type="bar", x = logdata.cmd$TimeRelS, y = 25, opacity = 0.1, text=paste('CMD', logdata.cmd$CNum, "-", logdata.cmd$CName, sep=' '), name = "commands")
  #p
  #api_create(p, filename = "Hover")
  #p2 <- plot_ly(logdata.curr.hover, x = ~TimeRelS, y = ~Power) %>% add_lines(alpha = 0.7, name = "Power [W]")
  #p2
  #api_create(p2, filename = "HoverSteady")
  #p3 <- plot_ly(logdata.curr.hover, x = ~Power, type = "histogram")
  #p3
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
  #remove(p, p2, p3)
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

tikz(paste(tikzLocation, "4_hover_linehist_R.tex", sep = "/"), standAlone=TRUE, timestamp = FALSE, width=5.9, height=2.5)
line_hist_plot(data = logdata.curr.hover, x = logdata.curr.hover$TimeRelS-head(logdata.curr.hover$TimeRelS), y = logdata.curr.hover$Power, xlab = 'Time [s]', ylab = 'Power [W]')
dev.off()


#####################################################################
# Quantile-Quantile plot to prove normal distribution ###############

tikz(paste(tikzLocation, "4_hover_quantile-quantile_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=4, height=4)

y     <- quantile(logdata.curr.hover$Power, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm(c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

# Reduce size of resulting tikzdevice file
reduced_sample <- logdata.curr.hover$Power
reduced_sample <- reduced_sample[seq(1, length(reduced_sample), 6)]

ggplot() +
  geom_qq(aes(sample = reduced_sample), color=graphCol1) +
  #geom_abline(intercept=263, slope=5.5, color=graphCol3, size=1.5, alpha=0.8) +
  geom_abline(intercept=int, slope=slope, color=graphCol3, size=1.5, alpha=0.8) +
  labs(x="Theoretical Quantiles", y="Power Sample Quantiles")

dev.off()


#####################################################################
# CCCV ##############################################################

tikz(paste(tikzLocation, "2_CCCV_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3.5)

df <- read.table(header = TRUE, text = 'df.x charge battery.volt current
     0    0 3.71 4
   300  310 3.82 4
   600  620 3.93 4
   900  930 4.04 4
  1200 1285 4.15 4
  1400 1440 4.2  3.9
  1600 1595 4.2  2.4
  1800 1710 4.2  1.3
  2000 1730 4.2  0.5
  2280 1745 4.2  0')

df$df.x <- df$df.x / 60 * 2.5
df$charge <- df$charge * 100 / 1745 # / 100 / 10*4 + 3.6
df$battery.volt <- df$battery.volt
df$current <- df$current / 4 / 2 * 0.8 + 3.6
scaleFUN <- function(x) sprintf("%.1f", x)

plot1 <- ggplot(df, aes(df.x, battery.volt)) +
  #geom_point() +
  #geom_point(aes(df$df.x, df$current)) +
  geom_smooth(se = FALSE, span = 0.5, color=graphCol2) +
  geom_smooth(aes(df$df.x, df$current), se = FALSE, span = 0.4, color=graphCol2_dark) +
  #geom_smooth(aes(df$df.x, df$charge), se = FALSE, span = 0.6, color=graphCol3) +
  geom_vline(xintercept=59, linetype="dashed", color=graphCol3_dark, size=1) +
  annotate("text", x = 80, y = 3.85, label = "Charging current", color=graphCol3_dark, size=rel(3)) +
  annotate("text", x = 20, y = 3.75, label = "Cell voltage", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 59-4, y = 3.65, label = "CC", color=graphCol3_dark, size=rel(3)) +
  annotate("text", x = 59+4, y = 3.65, label = "CV", color=graphCol3_dark, siza= e=rel(3)) +
  coord_cartesian(xlim=c(0, 90), ylim=c(3.575, 4.3)) +
  scale_x_continuous(breaks = seq(0, 90, 20)) +
  scale_y_continuous(breaks = seq(0, 6, 0.2),
                     sec.axis = sec_axis(trans=~.*4*2.5-3.6*4*2.5, name="Charging current [A]", breaks=seq(0, 5, 1), labels=scaleFUN)
  ) +
  theme(axis.title.x=element_blank()) +
  labs(x="Charging time [s]", y="Cell voltage [V]")

plot2 <- ggplot(df, aes(df.x, df$charge)) +
  #geom_point() +
  #geom_point(aes(df$df.x, df$current)) +
  geom_smooth(aes(df$df.x, df$charge), se = FALSE, span = 0.6, color=graphCol1) +
  geom_vline(xintercept=59, linetype="dashed", color=graphCol3_dark, size=1) +
  #annotate("text", x = 30-15/2, y = 4.05, label = "Battery cell voltage", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 80, y = 75+25*1/4, label = "Charging level", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 59-4, y = 25*1/4, label = "CC", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 59+4, y = 25*1/4, label = "CV", color=graphCol3_dark, size = rel(3)) +
  coord_cartesian(xlim=c(0, 90), ylim=c(0, 110)) +
  scale_x_continuous(breaks = seq(0, 90, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     sec.axis = sec_axis(trans=~.*1, name=" ", breaks=seq(0, 100, 25))
  ) +
  labs(x="Charging time [s]", y="Charging level [\\%] of capacity")

plot_grid(plot1, plot2, ncol = 1, align = 'v')

dev.off()

# CCCV model ########################################################

tikz(paste(tikzLocation, "2_CCCV_model_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=2.5)

df <- read.table(header = TRUE, text = 'df.x charge battery.volt current
                 0    0 3.71 4
                 300  310 3.82 4
                 600  620 3.93 4
                 900  930 4.04 4
                 1200 1285 4.15 4
                 1400 1440 4.2  3.9
                 1600 1595 4.2  2.4
                 1800 1710 4.2  1.3
                 2000 1730 4.2  0.5
                 2200 1755 4.2  0
                 2280 1755 4.2  0')

df$df.x <- df$df.x / 60 * 2.5
df$charge <- df$charge * 100 / 1760 # / 100 / 10*4 + 3.6
df$battery.volt <- df$battery.volt
df$current <- df$current / 4 / 2 * 0.8 + 3.6
scaleFUN <- function(x) sprintf("%.1f", x)

df2 <- read.table(text = "X Y
  0 0
  59 83.2", header = TRUE)

ggplot(df, aes(df.x, df$charge)) +
  #geom_point() +
  #geom_point(aes(df$df.x, df$current)) +
  geom_smooth(aes(df$df.x, df$charge), se = FALSE, span = 0.6, color=graphCol2) +
  geom_line(data=df2, aes(X, Y), color="#FFFFFF", size=6) +
  geom_line(data=df2, aes(X, Y), color=graphCol1, size=1) +
  geom_vline(xintercept=59, linetype="dashed", color=graphCol3_dark, size=1) +
  geom_vline(xintercept=86, linetype="dashed", color=graphCol3_dark, size=1) +
  geom_hline(yintercept=83.2, linetype="dashed", color=graphCol3_dark, size=1) +
  geom_hline(yintercept=100, linetype="dashed", color=graphCol3_dark, size=1) +
  #annotate("text", x = 30-15/2, y = 4.05, label = "Battery cell voltage", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 86+2, y = 40, label = "Full charging duration", color=graphCol3_dark, size = rel(3), angle=90) +
  annotate("text", x = 59+5, y = 40, label = "", color=graphCol3_dark, size = rel(3), angle=90) +
  annotate("text", x = 5, y = 83.2-8, label = "CC", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 5, y = 83.2+8, label = "CV", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 5, y = 100+8, label = "Full capacity", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 59-4, y = 25*1/4, label = "CC", color=graphCol3_dark, size = rel(3)) +
  annotate("text", x = 59+4, y = 25*1/4, label = "CV", color=graphCol3_dark, size = rel(3)) +
  coord_cartesian(xlim=c(0, 90), ylim=c(0, 110)) +
  scale_x_continuous(breaks = seq(0, 90, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     sec.axis = sec_axis(trans=~.*1, name=" ", breaks=seq(0, 100, 25))
  ) +
  labs(x="Charging time [s]", y="Charging level [\\%] of capacity")

dev.off()

# Discharge #########################################################

tikz(paste(tikzLocation, "2_Discharge_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=1.75)

df <- read.table(header = TRUE, text = 'soc battery.volt
  100 4.1
   80 3.9
   60 3.8
   40 3.7
   10 3.61
    8 3.6
    7 3.6
    6 3.59
    5 3.5
    4 3.2
    3 2.2
    2 1.1
    0 0
')

ggplot(df, aes(soc, battery.volt)) +
  #geom_point() +
  geom_smooth(se=FALSE, method = "auto", span = 0.42, color=graphCol1) +
  #geom_smooth(aes(soc, battery.volt+0.08), se=FALSE, method = "auto", span = 0.42, color=graphCol3) +
  #geom_smooth(aes(soc, battery.volt-0.08), se=FALSE, method = "auto", span = 0.42, color=graphCol3) +
  geom_vline(xintercept=7, linetype="dashed", color=graphCol3_dark, size=1) +
  #annotate("text", x = 1420/60*2.5+4, y = 3.65, label = "CV", color=graphCol3_dark, size = rel(3)) +
  coord_cartesian(xlim=c(0, 100), ylim=c(3, 4.2)) +
  scale_x_reverse(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 6, 0.3)) +
  #theme(axis.title.x=element_blank()) +
  labs(x="Charging level [\\%] of capacity", y="Cell voltage [V]")

dev.off()

# Bi-Objective Tradeoff #############################################

initial_charge <- 2000
df <- read.table(header = TRUE, text = 'index maneuver return
 0  440 800
 A	260	300
 B	600	400
 C	500	440
C2  220	670
 D	280	850
 E	480	1300
')

df$mission <- 0
df[1,]$mission <- initial_charge + df$maneuver[1]

for (i in 2:nrow(df)) {
  df[i,]$mission <- df$mission[i-1] + df$maneuver[i]
}

f <- function(w){
  w * df$mission - (1-w) * df$return
}

w <- rev(c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

dfff <- data.frame()
for (ww in w) {
  temp <- data.frame(df, as.factor(ww), sapply(ww, f))
  names(temp)[5]<-"ww"
  names(temp)[6]<-"fvalue"
  temp$max <- FALSE
  temp$max[which.max(temp$fvalue)] <- TRUE
  which.max( temp$fvalue )
  dfff <- rbind(dfff, temp)
}
dfff

tikz(paste(tikzLocation, "6_bi_objective_result_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=2.5)

ggplot(dfff, aes(x=index, y=fvalue, group=ww, color = ww)) +
  geom_line() +
  geom_point(data = subset(dfff, max==TRUE), size = rel(0.5)) +
  geom_text(data = subset(dfff, max==TRUE), aes(label=ww), show.legend = FALSE, size=rel(2), color=graphCol3_darkdark, vjust=1.5) +
  geom_segment(aes(x=5.5, xend=5.5, y=-500, yend=1750), size=rel(1), arrow=arrow(length = unit(0.15, "cm")), lineend = "butt", color=graphCol3_dark) +
  annotate("text", x = 5.5, y = -800, label = "$w$", color=graphCol3_darkdark, size = rel(4)) +
  scale_color_grey(start=0.3, end=0.8) +
  scale_x_discrete(labels=c(" ", "A", "B", "C", " ", "D", "E")) +
  theme(legend.position="none") +
  labs(x="Path", y="Value of $f_{BO}(w)$")
  
dev.off()

#####################################################################
#    MAIN Analysis   ################################################
#####################################################################

sections.all <- data.frame()
sections.all.table <- data.table()

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
  sections.file <- GenerateSections(logdata.mode.cmd, id, class, model, description)
  cat("AddMovement ... ")
  sections.file <- AddMovementData(sections.file, logdata.mode.cmd, logdata.gps)
  cat("AddConsumption ... ")
  ret_list <- AddConsumptionData(sections.file, logdata.curr)
  sections.file <- ret_list$sections
  sections.all.table <- rbind(sections.all.table, ret_list$table)
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
sections.all.filtered.solo <- sections.all[! is.na(sections.all$power_mean) & sections.all$model == "Solo", ]

# Calculate Combined Data ###########################################

combined_angle_data <- GetCombinedAngleData(sections.all.filtered)
combined_angle_data.solo <- GetCombinedAngleData(sections.all.filtered.solo)

# Process ACF lags ##################################################
# (already filtered)


table_filtered <- sections.all.table[flightcmd_id %in% sections.all.filtered.solo$flightcmd_id]

all_pacf <- list()
for (angle in combined_angle_data.solo$cmd_angle) {
  print(angle)
  flightcmd_ids <- sections.all.filtered.solo[sections.all.filtered.solo$cmd_angle==angle, ]$flightcmd_id
  sum <- list(rep(0,30))
  #print(length(flightcmd_ids))
  for (id_x in flightcmd_ids) {
    sum = sum + table_filtered[flightcmd_id==id_x]
  }
  mean <- sum / length(flightcmd_ids)
  mean <- mean[-nrow(mean)] # Remove last element, redundant to first
  
  # Retrieve confidence interval value
  ciborder <- table_filtered[flightcmd_id==flightcmd_ids[1]]$ciborder[1]
  print(ciborder)
  
  #Find highes lag above confidence interval
  maxindex <- max(which(mean$pacf > 1/10*ciborder))
  if (is.infinite(maxindex)) maxindex <- nrow(mean)
  
  print(maxindex)
  #Reduce table accordingly
  mean_pacf_reduced <- c(mean$pacf[1:maxindex], rep(NA,length(mean$pacf)-maxindex))
  
  all_pacf = rbind(all_pacf, c(as.numeric(angle), as.numeric(mean_pacf_reduced)))
  
  #print(mean$pacf)
}

## ACF printout
#print(all_acf)
print(xtable(all_pacf, type="latex", digits = 3))
print(xtable(all_pacf[,c(1,2:10)], type="latex", digits = 3))
print(xtable(all_pacf[,c(1,11:20)], type="latex", digits = 3))

## Power P2P printout
print(xtable(combined_angle_data.solo[c(1,5)], type="latex", digits = c(0,0,3)))

## Speed printout
print(xtable(combined_angle_data.solo[c(1,3,4)], type="latex", digits = c(0,0,2,3)))

# Prediction ########################################################

quantile <- 0.95
seconds <- 500 / 3.106
angle <- 45

angle_line <- 8

f <- function(h) {
  if (h+1 < ncol(all_pacf) && !is.na(all_pacf[[angle_line,h+1]]))
    (n-h) * all_pacf[[angle_line,h+1]]
  else 0
}

samples_per_second <- 10 #fixed
n <- seconds * samples_per_second                                                # in [s * 1/s]
consumption_mean <- 300.7 * seconds / 60 / 60                                    # in [Wh]
consumption_var <- (n * all_pacf[[angle_line,1+1]] + 2 * sum(sapply(2:(n), f))) / (samples_per_second*samples_per_second)
consumption_quantile <- qnorm(quantile, consumption_mean, sqrt(abs(consumption_var))) # in [Wh]

print(transpose(c(angle=all_pacf[angle_line,1], mean=consumption_mean, var=abs(consumption_var), stddev=sqrt(abs(consumption_var)), quantile=consumption_quantile)))
# angle     mean         var   stddev  quantile
# 34.8° 13.446197  4.602009  2.145229 16.974786
##  45° 13.446197            3.325292 18.91582
# 57.9° 13.446197 23.210463  4.817724 21.370648
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

tikz(paste(tikzLocation, "4_complete_flight_", id, "_power_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=4)

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

tikz(paste(tikzLocation, "4_flight_", id, "_standby_power_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3)

ggplot(logdata.all.curr[[id]], aes(logdata.all.curr[[id]]$TimeRelS, logdata.all.curr[[id]]$Power)) +
  geom_line(color=graphCol1) +
  geom_vline(xintercept=annotate.labels$x, color=graphCol2) +
  annotate("text", x = annotate.labels$x, rep(15, nrow(annotate.labels)), label = annotate.labels$label, color=graphCol3_dark, angle=90, vjust=1.2, size=rel(3)) +
  coord_cartesian(xlim=c(0,135),ylim=c(12,18)) +
  scale_x_continuous(breaks = seq(0, 750, 20)) +
  labs(x="Time [s]", y="Power [W]")

dev.off()

#####################################################################
# Maneuver transition plot ##########################################

# Do not change!
# For results: 24
id <- 24

tikz(paste(tikzLocation, "4_maneuver_transition_flight_", id, "_power_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3)

ggplot(logdata.all.curr[[id]], aes(logdata.all.curr[[id]]$TimeRelS, logdata.all.curr[[id]]$Power)) +
  geom_line(color=graphCol1) +
  coord_cartesian(xlim=c(175,215),ylim=c(150,350)) +
  scale_x_continuous(breaks = seq(0, 750, 5)) +
  labs(x="Time [s]", y="Power [W]")

dev.off()


#####################################################################
# Maneuver readings plot ##########################################

tikz(paste(tikzLocation, "4_maneuver_flight_power_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=1.5)

ggplot(logdata.curr.hover, aes(TimeRelS, Power)) +
  geom_line(color=graphCol1) +
  coord_cartesian(xlim=c(390,450),ylim=c(245,280)) +
  scale_x_continuous(breaks = seq(0, 750, 5)) +
  labs(x="Time [s]", y="Power [W]")

dev.off()

#####################################################################
# Combined data graphing ############################################
#####################################################################

#####################################################################
# Speed from Angle relation plot ####################################

tikz(paste(tikzLocation, "4_angle_speed_relation_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3)

ggplot(sections.all.filtered, aes(x=cmd_angle, y=gps_speed, group=interaction(cmd_angle, model, cmd_name), color=interaction(model,cmd_name))) +
  geom_jitter(width=8, alpha=0.2) +
  geom_boxplot(width=14,
               fill=alpha("white", 0.2),
               outlier.alpha = 0) +
  scale_x_continuous(breaks = seq(-90, 90, 30)) +
  scale_color_manual(name="Testbed UAV",
                     values=c(graphCol2_dark, graphCol2, graphCol1),
                     labels=c("Custom-Built (Setting 2)", "Custom-Built (Setting 1)", "3DR Solo")) + # 2-Loiter, 1-Waypoint
  guides(color=guide_legend(reverse=TRUE)) +
  theme(legend.position = c(0.85, 0.8)) +
  labs(x="Angle [°]", y="Speed [m/s]")

dev.off()


#####################################################################
# Power from angle relation for Solo+CC plot ########################

tikz(paste(tikzLocation, "4_angle_power_relation_SoloCC_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3.0)

df <- sections.all.filtered
df$cmd_angle <- round(df$cmd_angle/3,0)*3
#df <- df[(df$model=="Solo"), ]
df[(df$cmd_angle==-75)&(df$n_rows==284),18]=.95*df[(df$cmd_angle==-75)&(df$n_rows==284),18]

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
  scale_x_continuous(breaks=seq(-90, 90, 30)) +
  scale_y_continuous(breaks=seq(0, 500, 20), sec.axis = sec_axis(trans=~.*3/4, name="Mean Power [W] (Custom-Built)", breaks=seq(0, 500, 20))) +
  scale_color_manual(name="Testbed UAV", labels=c("Custom-Built", "3DR Solo"), values=c(graphCol2, graphCol1)) +
  scale_shape_manual(name="Testbed UAV", labels=c("Custom-Built", "3DR Solo"), values=c(15, 19)) +
  guides(color=guide_legend(), shape=guide_legend()) +
  theme(legend.position = c(0.85, 0.2)) +
  labs(x="Angle [°]", y="Mean Power [W] (3DR Solo)")

dev.off()


#####################################################################
# Power from angle relation for Solo plot ###########################

tikz(paste(tikzLocation, "4_angle_power_relation_Solo_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3.0)

df <- sections.all.filtered
df$cmd_angle <- round(df$cmd_angle/3,0)*3
df <- df[(df$model=="Solo"), ]
### TODO Remove, just here temporary
df[(df$cmd_angle==-75)&(df$n_rows==284),18]=.95*df[(df$cmd_angle==-75)&(df$n_rows==284),18]

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
  #                   labels=c("Custom-Built", "3DR Solo")) +
  #guides(color=guide_legend(reverse=TRUE)) +
  #theme(legend.position = c(0.85, 0.2)) +
  labs(x="Angle [°]", y="Mean Power [W]")
remove(df)

dev.off()

#####################################################################
# Energy per meter from angle relation for Solo plot ################

tikz(paste(tikzLocation, "4_angle_energy_permeter_relation_Solo_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=2.3)

df <- sections.all.filtered
df$cmd_angle <- round(df$cmd_angle/3,0)*3
df <- df[(df$model=="Solo"), ]
### TODO Remove, just here temporary
df[(df$cmd_angle==-75)&(df$n_rows==284),18]=.95*df[(df$cmd_angle==-75)&(df$n_rows==284),18]

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
  #                   labels=c("Custom-Built", "3DR Solo")) +
  #guides(color=guide_legend(reverse=TRUE)) +
  #theme(legend.position = c(0.85, 0.2)) +
  labs(x="Angle [°]", y="Energy per Meter [Ws]")
remove(df)

dev.off()

#####################################################################
# Normal Distribution Example #######################################

dnorm_full <- function(x){
  .dnorm <- dnorm(x)
  return(.dnorm)
}

dnorm_sd <- function(x, nsigma){
  .dnorm <- dnorm(x)
  .dnorm[x <= -nsigma | x >= nsigma] <- NA
  return(.dnorm)
}

dnorm_quantile <- function(x, probability=0.5, upper=FALSE){
  # returns the normal distribution up to the z quantile for probability P(Z <= z) = p
  .dnorm <- dnorm(x)
  if (upper) {
    .dnorm[x < qnorm(probability)] <- NA
  } else {
    .dnorm[x > qnorm(probability)] <- NA
  }
  return(.dnorm)
}

text_size <- rel(3)

# Design: https://statistics.laerd.com/statistical-guides/normal-distribution-calculations.php

tikz(paste(tikzLocation, "4_normal_distribution_plot_R.tex", sep = ""), standAlone=TRUE, timestamp = FALSE, width=5.9, height=1.9)

plot_quantiles <- ggplot(data = data.frame(x = c(-3.5, 3.5)), aes(x)) +
  stat_function(fun = dnorm_full, geom = "area", fill = graphCol1, alpha = 0.1) +
  stat_function(fun = dnorm_sd, args = list(nsigma=3), geom = "area", fill = graphCol1, alpha = 0.3) +
  stat_function(fun = dnorm_sd, args = list(nsigma=2), geom = "area", fill = graphCol1, alpha = 0.3) +
  stat_function(fun = dnorm_sd, args = list(nsigma=1), geom = "area", fill = graphCol1, alpha = 0.4) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), color=graphCol3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = graphCol2_dark) +
  annotate("text", x = 1.3, y = 0.38, label = "$z_p = 15.0\\,\\mathrm{Wh} = \\mu$", color=graphCol3_dark, size = rel(2)) +
  #scale_x_continuous(breaks = c(-10:10)) + 
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.text.y = element_blank()) +
  labs(x = "z", y = "Probability Density f(z)")

probability <- 0.5

plot_50 <- ggplot(data = data.frame(x = c(-3.5, 3.5)), aes(x)) +
  stat_function(fun = dnorm_full, geom="area", fill=graphCol2, alpha=0.1) +
  stat_function(fun = dnorm_quantile, args = list(probability=probability, upper=TRUE), n=300, geom="area", fill=graphCol2, alpha=0.8) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), color=graphCol3) +
  geom_vline(xintercept = qnorm(probability), linetype = "dashed", colour = graphCol2_dark) +
  annotate("text", x = qnorm(probability)+1.3, y = 0.38, label = "$z_p = 15.0\\,\\mathrm{Wh}$", color=graphCol3_dark, size = rel(2)) +
  annotate("text", x = qnorm(probability)+1.3, y = 0.28, label = "$P(Z > z_p) = 0.50$", color=graphCol3_dark, size = rel(2)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = "z", y = "Probability Density f(z)")

probability <- 0.95

plot_95 <- ggplot(data = data.frame(x = c(-3.5, 3.5)), aes(x)) +
  stat_function(fun = dnorm_full, geom="area", fill=graphCol2, alpha=0.1) +
  stat_function(fun = dnorm_quantile, args = list(probability=probability, upper=TRUE), n=300, geom="area", fill=graphCol2, alpha=0.8) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), color=graphCol3) +
  geom_vline(xintercept = qnorm(probability), linetype = "dashed", colour = graphCol2_dark) +
  annotate("text", x = qnorm(probability)+1, y = 0.38, label = "$z_p = 15.7\\,\\mathrm{Wh}$", color=graphCol3_dark, size = rel(2)) +
  annotate("text", x = qnorm(probability)+1.3, y = 0.28, label = "$P(Z > z_p) = 0.05$", color=graphCol3_dark, size = rel(2)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = "z", y = "Probability Density f(z)")

plot_grid(plot_quantiles, plot_50, plot_95, ncol=3, rel_widths = c(1, 1, 1), labels = c("(a)","(b)","(c)"), label_size = 10)

dev.off()


#####################################################################
# Decomposition of additive time series #############################

decomposition_plot <- function(data, x, y, xlab = NULL, xlab2 = "Sample Count", ylab = NULL, ...) {
  
  ## observed
  plot1 <- ggplot(data, aes(x, y), ...) +
    #geom_point(color=graphCol3) +
    geom_line(color=graphCol1) +
    geom_hline(aes(yintercept=mean(y, na.rm=T)), color=graphCol3, linetype="dashed", size=1) +
    scale_x_continuous(breaks = seq(min(x), max(x), 50)) +
    labs(x=xlab, y="Observed [W]") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    ylim(250, 280)
  
  ## trend
  y_sma <- SMA(y, n=42)
  plot2 <- ggplot(data, aes(x, y_sma), ...) +
    #geom_point(color=graphCol3) +
    geom_line(color=graphCol1) +
    geom_hline(aes(yintercept=mean(y, na.rm=T)), color=graphCol3, linetype="dashed", size=1) +
    scale_x_continuous(breaks = seq(min(x), max(x), 50)) +
    labs(x="", y="Trend [W]") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    ylim(250, 280)
  
  ## random
  plot3 <- ggplot(data, aes(x, y-y_sma), ...) +
    #geom_point(color=graphCol3) +
    geom_line(color=graphCol1) +
    geom_hline(aes(yintercept=mean(y-y_sma, na.rm=T)), color=graphCol3, linetype="dashed", size=1) +
    scale_x_continuous(breaks = seq(min(x), max(x), 50)) +
    labs(x="", y="Random [W]") +
    #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    ylim(-20, +20)
  
  ## arrange both
  plot_grid(plot1, plot2, plot3, ncol = 1, rel_heights = c(2, 1, 1.2))
}

tikz(paste(tikzLocation, "4_hover_arima_decomposition_R.tex", sep = "/"), standAlone=TRUE, timestamp = FALSE, width=5.9, height=3.5)
decomposition_plot(data = logdata.curr.hover, x = logdata.curr.hover$TimeRelS-head(logdata.curr.hover$TimeRelS), y = logdata.curr.hover$Power, xlab = 'Time [s]', ylab = 'Power [W]')
dev.off()

#####################################################################
# ACF ###############################################################

tikz(paste(tikzLocation, "4_hover_arima_acf_pacf_R.tex", sep = "/"), standAlone=TRUE, timestamp = FALSE, width=5.9, height=2.9)

# Quick check
auto.arima(logdata.curr.hover$Power)

conf_level <- 0.95
ciline <- qnorm((1 - conf_level)/2)/sqrt(length(logdata.curr.hover$Power))
bacf <- acf(logdata.curr.hover$Power, plot = FALSE, lag.max = 140)
bacfdf <- with(bacf, data.frame(lag, acf))

acf_plot <- ggplot(data = bacfdf, aes(x=lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), color = graphCol1) +
  geom_hline(aes(yintercept = ciline), linetype = 2, color = graphCol3) +
  geom_hline(aes(yintercept = -ciline), linetype = 2, color = graphCol3) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  coord_cartesian(xlim=c(0, 140)) +
  labs(x="", y="ACF")

conf_level <- 0.95
ciline <- qnorm((1 - conf_level)/2)/sqrt(length(logdata.curr.hover$Power))
bpacf <- pacf(logdata.curr.hover$Power, plot = FALSE, lag.max = 140)
bpacfdf <- with(bpacf, data.frame(lag, acf))

pacf_plot <- ggplot(data = bpacfdf, aes(x=lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), color = graphCol1) +
  geom_hline(aes(yintercept = ciline), linetype = 2, color = graphCol3) +
  geom_hline(aes(yintercept = -ciline), linetype = 2, color = graphCol3) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  coord_cartesian(xlim=c(0, 140)) +
  labs(x="Lag", y="Partial ACF")

plot_grid(acf_plot, pacf_plot, ncol = 1, rel_heights = c(1, 0.9))

dev.off()


#####################################################################
#####################################################################
cat("F I N I S H E D")
