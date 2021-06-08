####################################
##### Cleaning the environment  ####

rm(list=ls())

####################################
##### Loading packages  ############

library(readxl)
library(dplyr)
library(stringr)
library(rstudioapi)

####################################
##### Setting working directory ####

setwd(dirname(getActiveDocumentContext()$path)) 


####################################
##### Creating data-set  ###########



filename <- list.files("Maps")

atlasmaps <- str_split_fixed(filename, ",", n=2) 
colnames(atlasmaps) <- c("citekey","title")
atlasmaps <- cbind(atlasmaps, done = 0, num.maps = 1, year = -99, hyear = -99,
                   lyear = -99, state1 = c(""), state2 = c(""), state3 = c(""),
                   state4 = c(""), state5 = c(""), state6 = (""),  state7 = c(""),
                   state8 = c(""), state9 = (""), state10 = c(""))

write.csv2(atlasmaps, file = "atlasmaps.csv", row.names=F)

