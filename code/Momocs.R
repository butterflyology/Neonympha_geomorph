# 2021-09-30
# Chris Hamm
# Momocs analysis of Neonympha wing shape


#Preliminaries
set.seed(1138)

library(tidyverse)
library(Momocs)
library(geomorph)


#Load data

### Pattern raw data
Pattern_raw <- read.csv("data/geomorph/Pattern_raw.csv", header = TRUE)
dim(Pattern_raw)


### Structure raw data
Structure_raw <- read.csv("data/geomorph/Structure_raw.csv", header = TRUE)
dim(Structure_raw)


### Covariate raw data
Covariate_data <- read.csv("data/geomorph/BO0_12_Dec.csv", header = TRUE)
dim(Covariate_data)


### Merging data
Landmark_data <- merge(x = Pattern_raw, y = Structure_raw, by = "Id")
dim(Landmark_data)

Combined_data <- merge(x = Landmark_data, y = Covariate_data, by = "Id")
dim(Combined_data)

# Prepare structure data for analysis
Structure_array <- arrayspecs(Combined_data[, 54:83], p = 15, k = 2)

Structure_gpa <- gpagen(Structure_array, PrinAxes = FALSE, Proj = TRUE, ProcD = TRUE, print.progress = FALSE)

dimnames(Structure_array)[[3]] <- Combined_data$Taxon
Structure_2d <- two.d.array(Structure_gpa$coords)
rownames(Structure_2d) <- Combined_data$Taxon


