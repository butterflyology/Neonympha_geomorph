# Fourier analysis on butterfly wing shape

# 2024-02-18
# Chris Hamm


# Preliminaries ----

## Load packages

set.seed(8762432)
library("geomorph")
library("Momocs")



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


# Procrustes superimposition ----
### Pattern data
Pattern_array <- arrayspecs(Combined_data[, 2:53], p = 26, k = 2)

Pattern_gpa <- gpagen(Pattern_array, PrinAxes = FALSE, Proj = TRUE, ProcD = TRUE, print.progress = FALSE)

dimnames(Pattern_array)[[3]] <- Combined_data$Taxon
Pattern_2d <- two.d.array(Pattern_gpa$coords)
rownames(Pattern_2d) <- Combined_data$Taxon


### Structure data
Structure_array <- arrayspecs(Combined_data[, 54:83], p = 15, k = 2)

Structure_gpa <- gpagen(Structure_array, PrinAxes = FALSE, Proj = TRUE, ProcD = TRUE, print.progress = FALSE)

dimnames(Structure_array)[[3]] <- Combined_data$Taxon
Structure_2d <- two.d.array(Structure_gpa$coords)
rownames(Structure_2d) <- Combined_data$Taxon

str(Structure_gpa)

# Need to create a subset of data that jsut contain the outline of the wing, no internal componants.
outer_landmarks <- Structure_2d[, c(1, 2, 3, 4, 5, 6, 9, 10, 13, 14, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30) ]

head(outer_landmarks)
str(outer_landmarks)
dim(outer_landmarks)

# Convert for fourier analysis

# Convert the landmark data to a dataframe
outer_landmarks_df <- as.data.frame(outer_landmarks)
head(outer_landmarks_df)

# Extracting coordinates and converting them into a format suitable for Momocs
outlines_list <- lapply(1:nrow(outer_landmarks_df), function(i) {
  as.matrix(outer_landmarks_df[i, ])
})

species_labels <- attr(outer_landmarks, "dimnames")[[1]]



