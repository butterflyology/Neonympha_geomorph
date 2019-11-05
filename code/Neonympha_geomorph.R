# Chris Hamm 2019-11-04
# R code for geometric morphometrics


# Preliminaries ----

## Load packages

set.seed(8762432)
library("geomorph")
library("MASS")
library("spaceMovie")

# save package descriptions as text file
# sink("misc/Neonympha_geommorph_packages.txt")
# sessionInfo()
# sink()

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


# Data Visualization ----

## Pattern data
N_ar_p_mean <- colMeans(Pattern_2d[Combined_data$Taxon == "N.ar", ])
N_he_p_mean <- colMeans(Pattern_2d[Combined_data$Taxon == "N.he", ])
N_fr_p_mean <- colMeans(Pattern_2d[Combined_data$Taxon == "N.fr", ])
N_mi_p_mean <- colMeans(Pattern_2d[Combined_data$Taxon == "N.mi", ])

plot(colMeans(Pattern_2d[c(seq(1, 51, 2)), ]), colMeans(Pattern_2d[c(seq(2, 52, 2)), ]), type = 'n', xlab = 'X Coordinate', ylab ='Y Coordinate', xlim = c(-0.35, 0.25), ylim = c(-0.35, 0.35), las = 1)
points(as.matrix(Pattern_2d[, c(seq(1, 51, 2)) ]), pch = 19, as.matrix(Pattern_2d[, c(seq(2, 52, 2)) ]), col = rgb(0.745098, 0.745098, 0.745098, 0.25))
points(colMeans(Pattern_2d[, c(seq(1, 51, 2)) ]), pch = 19, colMeans(Pattern_2d[, c(seq(2, 52, 2)) ]), col = "black")	
points(N_ar_p_mean[ c(seq(1, 51, 2))] , N_ar_p_mean[ c(seq(2, 52, 2))] , pch = 19, col = SW_palette("Main")[1]) 
points(N_he_p_mean[ c(seq(1, 51, 2))] , N_he_p_mean[ c(seq(2, 52, 2))] , pch = 19, col = SW_palette("Inquisitor")[1]) 
points(N_fr_p_mean[ c(seq(1, 51, 2))] , N_fr_p_mean[ c(seq(2, 52, 2))] , pch = 19, col = SW_palette("ROTJ")[3]) 
points(N_mi_p_mean[ c(seq(1, 51, 2))] , N_mi_p_mean[ c(seq(2, 52, 2))] , pch = 19, col = SW_palette("TESB")[8]) 
legend('bottomleft', legend = c('Individual landmarks','Grand mean', expression(paste(italic('N. areolatus'))), expression(paste(italic('N. helicta'))), expression(paste(italic('N. m. fransisci'))), expression(paste(italic('N. m. mitchellii')))), col = c('grey','black', SW_palette("Main")[1], SW_palette("Inquisitor")[1], SW_palette("ROTJ")[3], SW_palette("TESB")[8]), pch = 19, bty = 'n')


## Structure data
N_ar_s_mean <- colMeans(Structure_2d[Combined_data$Taxon == "N.ar", ])
N_he_s_mean <- colMeans(Structure_2d[Combined_data$Taxon == "N.he", ])
N_fr_s_mean <- colMeans(Structure_2d[Combined_data$Taxon == "N.fr", ])
N_mi_s_mean <- colMeans(Structure_2d[Combined_data$Taxon == "N.mi", ])

plot(colMeans(Structure_2d[c(seq(1, 29, 2)), ]), colMeans(Structure_2d[c(seq(2, 30, 2)), ]), type = 'n', xlab = 'X Coordinate', ylab ='Y Coordinate', xlim = c(-0.4, 0.3), ylim = c(-0.4, 0.4), las = 1)
points(as.matrix(Structure_2d[, c(seq(1, 29, 2)) ]), pch = 19, as.matrix(Structure_2d[, c(seq(2, 30, 2)) ]), col = rgb(0.745098, 0.745098, 0.745098, 0.25))
points(colMeans(Structure_2d[, c(seq(1, 29, 2)) ]), pch = 19, colMeans(Structure_2d[, c(seq(2, 30, 2)) ]), col = "black")	
points(N_ar_s_mean[ c(seq(1, 29, 2))] , N_ar_s_mean[ c(seq(2, 30, 2))] , pch = 19, col = SW_palette("Main")[1])
points(N_he_s_mean[ c(seq(1, 29, 2))] , N_he_s_mean[ c(seq(2, 30, 2))] , pch = 19, col = SW_palette("Inquisitor")[1]) 
points(N_fr_s_mean[ c(seq(1, 29, 2))] , N_fr_s_mean[ c(seq(2, 30, 2))] , pch = 19, col = SW_palette("ROTJ")[3]) 
points(N_mi_s_mean[ c(seq(1, 29, 2))] , N_mi_s_mean[ c(seq(2, 30, 2))] , pch = 19, col = SW_palette("TESB")[8]) 
legend('topleft', legend = c('Individual landmarks','Grand mean', expression(paste(italic('N. areolatus'))), expression(paste(italic('N. helicta'))), expression(paste(italic('N. m. fransisci'))), expression(paste(italic('N. m. mitchellii')))), col = c('grey','black', SW_palette("Main")[1], SW_palette("Inquisitor")[1], SW_palette("ROTJ")[3], SW_palette("TESB")[8]), pch = 19, bty = 'n')


# Analyses ----
## Linear Discriminant Analysis

### Pattern
# The first step is to conduct a Principle Componants Analysis (PCA)
Pattern_pca <- prcomp(Pattern_2d)
summary(Pattern_pca)$importance[, 1:4]

# Recall that the last four PC dimensions (49:52 in this case) are empty when dealing with morphometric data - we remove these to avoid deficient dimensions

Pattern_pca_shape <- cbind(Combined_data[, c(1, 84)], Pattern_pca$x[, 1:48])

Pattern_lda <- lda(as.matrix(Pattern_pca_shape[, 3:50]), Pattern_pca_shape$Taxon, method ="mle" )

Pattern_lda_scores <- as.matrix(Pattern_pca_shape[, 3:50]) %*% as.matrix(Pattern_lda$scaling) # Here we use matrix multiplication to multiply the Pattern pcs by the lda scaling. This gives us LDA positions for each individual.
Pattern_pca_lda <- cbind(Pattern_pca_shape, Pattern_lda_scores)

plot(Pattern_pca_lda$LD1, Pattern_pca_lda$LD2, col = c(SW_palette("Main")[1], SW_palette("ROTJ")[3], SW_palette("Inquisitor")[1], SW_palette("TESB")[8])[Pattern_pca_lda$Taxon], pch = 19, xlab= 'LD1 (75.5%)', ylab ='LD 2 (15.5%)', xlim =c(-5, 6), ylim = c(-4.5, 4), cex = 2, las = 1 ) # % of variance explained at bottom of output when calling Pattern_lda
legend('bottomleft', legend = c(expression(paste(italic('N. areolatus'))), expression(paste(italic('N. helicta'))), expression(paste(italic('N. m. fransisci'))), expression(paste(italic('N. m. mitchellii')))), bty = 'n', pch = 19, col = c(SW_palette("Main")[1], SW_palette("Inquisitor")[1], SW_palette("ROTJ")[3], SW_palette("TESB")[8]), pt.cex = 2)


### Structure 
Structure_pca <- prcomp(Structure_2d)
summary(Structure_pca)$importance[, 1:4]

Structure_pca_shape <- cbind(Combined_data[, c(1, 84)], Structure_pca$x[, 1:26])

Structure_lda <- lda(as.matrix(Structure_pca_shape[, 3:28]), Structure_pca_shape$Taxon, method ="mle" )

Structure_lda_scores <- as.matrix(Structure_pca_shape[, 3:28]) %*% as.matrix(Structure_lda$scaling)
Structure_pca_lda <- cbind(Structure_pca_shape, Structure_lda_scores)

plot(Structure_pca_lda$LD1, Structure_pca_lda$LD2, col = c(SW_palette("Main")[1], SW_palette("ROTJ")[3], SW_palette("Inquisitor")[1], SW_palette("TESB")[8])[Structure_pca_lda$Taxon], pch = 19, xlab= 'LD1 (74.7%)', ylab ='LD 2 (1.5%)', xlim =c(-5, 6), ylim = c(-4.5, 4), cex = 2, las = 1 )
legend('bottomleft', legend = c(expression(paste(italic('N. areolatus'))), expression(paste(italic('N. helicta'))), expression(paste(italic('N. m. fransisci'))), expression(paste(italic('N. m. mitchellii')))), bty = 'n', pch = 19, col = c(SW_palette("Main")[1], SW_palette("Inquisitor")[1], SW_palette("ROTJ")[3], SW_palette("TESB")[8]), pt.cex = 2)


### Covariate
# Not all individuals have border ocelli in cells 1 and 2, so we will omit those from this analysis. We will select the area, lenght and width of border ocelli from cells 3:6.

Cov_ALW <- Combined_data[, c(87:90, 95:102)]
Cov_ALW_pca <- prcomp(Cov_ALW)
summary(Cov_ALW_pca)$importance[, 1:4]

Cov_ALW_shape <- cbind(Combined_data[, c(1, 84)], Cov_ALW_pca$x) # These data are not morphometric, so we don't need to remove dimensions.

Cov_ALW_lda <- lda(as.matrix(Cov_ALW_shape[, 3:14]), Cov_ALW_shape$Taxon, method = "mle")

Cov_ALW_lda_scores <- as.matrix(Cov_ALW_shape[, 3:14]) %*% as.matrix(Cov_ALW_lda$scaling)

Cov_ALW_pca_lda <- cbind(Cov_ALW_shape, Cov_ALW_lda_scores)

plot(Cov_ALW_pca_lda$LD1, Cov_ALW_pca_lda$LD2, col = c(SW_palette("Main")[1], SW_palette("ROTJ")[3], SW_palette("Inquisitor")[1], SW_palette("TESB")[8])[Cov_ALW_pca_lda$Taxon], xlim= c(-4.5, 4), ylim = c(-4, 4),pch = 19, xlab= "LD1 (83.8%)", ylab ="LD 2 (12.8%)", cex = 2, las = 1 )
legend('topleft', legend = c(expression(paste(italic('N. areolatus'))), expression(paste(italic('N. helicta'))), expression(paste(italic('N. m. fransisci'))), expression(paste(italic('N. m. mitchellii')))), bty = 'n', pch = 19, col = c(SW_palette("Main")[1], SW_palette("Inquisitor")[1], SW_palette("ROTJ")[3], SW_palette("TESB")[8]), pt.cex = 2)


## Discriminant analysis

### Pattern
Pattern_DF <- predict(Pattern_lda)
Pattern_pca_DF <- cbind(Pattern_pca_lda, Pattern_DF$posterior)

N_ar_P_DF <- Pattern_pca_DF[Pattern_pca_DF$Taxon == "N.ar", ]
N_he_P_DF <- Pattern_pca_DF[Pattern_pca_DF$Taxon == "N.he", ]
N_fr_P_DF <- Pattern_pca_DF[Pattern_pca_DF$Taxon == "N.fr", ]
N_mi_P_DF <- Pattern_pca_DF[Pattern_pca_DF$Taxon == "N.mi", ]

par(mfrow = c(2, 2))
hist(N_ar_P_DF$N.ar, col = SW_palette("Main")[1], main = "", xlab = "", ylab = "Frequency", ylim = c(0, 100), las = 1, cex.lab = 1.5)
text(0.3, 80, expression(paste(bolditalic("N. areolatus"))), cex = 1.5)
abline(v = mean(N_ar_P_DF$N.ar), lty = 2, lwd = 2)

hist(N_he_P_DF$N.he, col = SW_palette("Inquisitor")[1], ylab = "", xlab = "", main = "", ylim = c(0, 20), las = 1)
text(0.25, 15, expression(paste(bolditalic("N. helicta"))), cex = 1.5)
abline(v = mean(N_he_P_DF$N.he), lty = 2, lwd = 2)

hist(N_fr_P_DF$N.fr, col = SW_palette("ROTJ")[3], main = "", xlab = "Posterior probability", ylab = "Frequency", ylim = c(0, 15), xlim = c(0, 1), las = 1, cex.lab = 1.5)
text(0.35, 12, expression(paste(bolditalic("N. m. francisci"))), cex = 1.5)
abline(v = mean(N_fr_P_DF$N.fr), lty = 2, lwd = 2)

hist(N_mi_P_DF$N.mi, col = SW_palette("TESB")[8], main = "", ylab = "", xlab = "Posterior probability", ylim = c(0, 60), xlim = c(0, 1), las = 1, cex.lab = 1.5)
text(0.35, 50, expression(paste(bolditalic("N. m. mitchellii"))), cex = 1.5)
abline(v = mean(N_mi_P_DF$N.mi), lty = 2, lwd = 2)
par(mfrow = c(1, 1))


### Structure
Structure_DF <- predict(Structure_lda)
Structure_pca_DF <- cbind(Structure_pca_lda, Structure_DF$posterior)

N_ar_S_DF <- Structure_pca_DF[Structure_pca_DF$Taxon == "N.ar", ]
N_he_S_DF <- Structure_pca_DF[Structure_pca_DF$Taxon == "N.he", ]
N_fr_S_DF <- Structure_pca_DF[Structure_pca_DF$Taxon == "N.fr", ]
N_mi_S_DF <- Structure_pca_DF[Structure_pca_DF$Taxon == "N.mi", ]

par(mfrow = c(2, 2))
hist(N_ar_S_DF$N.ar, col = SW_palette("Main")[1], main = "", xlab = "", ylab = "Frequency", ylim = c(0, 100), las = 1, cex.lab = 1.5)
text(0.3, 80, expression(paste(bolditalic("N. areolatus"))), cex = 1.5)
abline(v = mean(N_ar_S_DF$N.ar), lty = 2, lwd = 2)

hist(N_he_S_DF$N.he, col = SW_palette("Inquisitor")[1], ylab = "", xlab = "", main = "", ylim = c(0, 20), las = 1)
text(0.25, 15, expression(paste(bolditalic("N. helicta"))), cex = 1.5)
abline(v = mean(N_he_S_DF$N.he), lty = 2, lwd = 2)

hist(N_fr_S_DF$N.fr, col = SW_palette("ROTJ")[3], main = "", xlab = "Posterior probability", ylab = "Frequency", ylim = c(0, 15), xlim = c(0, 1), las = 1, cex.lab = 1.5)
text(0.35, 12, expression(paste(bolditalic("N. m. francisci"))), cex = 1.5)
abline(v = mean(N_fr_S_DF$N.fr), lty = 2, lwd = 2)

hist(N_mi_S_DF$N.mi, col = SW_palette("TESB")[8], main = "", ylab = "", xlab = "Posterior probability", ylim = c(0, 60), xlim = c(0, 1), las = 1, cex.lab = 1.5)
text(0.35, 50, expression(paste(bolditalic("N. m. mitchellii"))), cex = 1.5)
abline(v = mean(N_mi_S_DF$N.mi), lty = 2, lwd = 2)
par(mfrow = c(1, 1))


#### Covariate
Cov_ALW_DF <- predict(Cov_ALW_lda)
Cov_ALW_pca_DF <- cbind(Cov_ALW_pca_lda, Cov_ALW_DF$posterior)

N_ar_Cov_DF <- Cov_ALW_pca_DF[Cov_ALW_pca_DF$Taxon == "N.ar", ]
N_he_Cov_DF <- Cov_ALW_pca_DF[Cov_ALW_pca_DF$Taxon == "N.he", ]
N_fr_Cov_DF <- Cov_ALW_pca_DF[Cov_ALW_pca_DF$Taxon == "N.fr", ]
N_mi_Cov_DF <- Cov_ALW_pca_DF[Cov_ALW_pca_DF$Taxon == "N.mi", ]

par(mfrow = c(2, 2))
hist(N_ar_Cov_DF$N.ar, col = SW_palette("Main")[1], main = "", xlab = "", ylab = "Frequency", ylim = c(0, 100), las = 1, cex.lab = 1.5)
text(0.3, 80, expression(paste(bolditalic("N. areolatus"))), cex = 1.5)
abline(v = mean(N_ar_Cov_DF$N.ar), lty = 2, lwd = 2)

hist(N_he_Cov_DF$N.he, col = SW_palette("Inquisitor")[1], ylab = "", xlab = "", main = "", ylim = c(0, 20), las = 1)
text(0.25, 15, expression(paste(bolditalic("N. helicta"))), cex = 1.5)
abline(v = mean(N_he_Cov_DF$N.he), lty = 2, lwd = 2)

hist(N_fr_Cov_DF$N.fr, col = SW_palette("ROTJ")[3], main = "", xlab = "Posterior probability", ylab = "Frequency", ylim = c(0, 15), xlim = c(0, 1), las = 1, cex.lab = 1.5)
text(0.35, 12, expression(paste(bolditalic("N. m. francisci"))), cex = 1.5)
abline(v = mean(N_fr_Cov_DF$N.fr), lty = 2, lwd = 2)

hist(N_mi_Cov_DF$N.mi, col = SW_palette("TESB")[8], main = "", ylab = "", xlab = "Posterior probability", ylim = c(0, 60), xlim = c(0, 1), las = 1, cex.lab = 1.5)
text(0.35, 50, expression(paste(bolditalic("N. m. mitchellii"))), cex = 1.5)
abline(v = mean(N_mi_Cov_DF$N.mi), lty = 2, lwd = 2)
par(mfrow = c(1, 1))


## Procrustes ANOVA

### Create the geomorph.data.frame
Neonympha_data.frame <- geomorph.data.frame(Pattern = Pattern_gpa$coords, Pattern_logCS = log(Pattern_gpa$Csize), Structure = Structure_gpa$coords, Structure_logCS = log(Structure_gpa$Csize), Taxon = Combined_data$Taxon, BO_size = Combined_data[, 87:90], BO_LW = Combined_data[, 95:102], Lat = Combined_data$Lat)
attributes(Neonympha_data.frame)

### Pattern

#### Pattern and by taxon
# With recent update to geomorph now need to compare with anova and pairwise after two procD.lm's?
plm_1a <- procD.lm(Pattern ~ Taxon, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(plm_1a)
reveal.model.designs(plm_1a)

plm_1b <- procD.lm(Pattern ~ 1, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(plm_1b)

anova(plm_1b, plm_1a) # including pattern improves analysis

# order matters here, the null model is second by default
plm_1_pw <- pairwise(fit = plm_1a, fit.null = plm_1b, groups = Neonympha_data.frame$Taxon)
summary(plm_1_pw, test.type = "dist", confidence = 0.95) # really cool; sofa king cool


#### Pattern by taxon and size (factor interaction)
plm_2a <- procD.lm(Pattern ~ Taxon * Pattern_logCS, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(plm_2a) # Taxon sig, no interaction

plm_2b <- procD.lm(Pattern ~ Taxon + Pattern_logCS, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE) # Taxon significant
summary(plm_2b) # taxon sig, log pattern nothing

anova(plm_2b, plm_2a) # no difference between models, and both suck. No need to do further testing.


#### Pattern centroid size by taxon and latitude
plm_3a <- procD.lm(Pattern_logCS ~ Taxon * Lat, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(plm_3a) # Taxon interaction significant

plm_3b <- procD.lm(Pattern_logCS ~ Taxon + Lat, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(plm_3b) # additive not significant

anova(plm_3b, plm_3a) # interaction model sig better

plm_3_pw <- pairwise(fit = plm_3b, fit.null = plm_3b, groups = Neonympha_data.frame$Taxon)

summary(plm_3_pw, test.type = "dist", confidence = 0.95) # can't distinguish taxa



### Structure
#### Structure by taxon
# apd_2 <- advanced.procD.lm(f1 = Structure ~ 1, Structure ~ Taxon, data = Neonympha_data.frame, iter = 1e4, seed = 76234, print.progress = FALSE)
# summary(apd_2)

slm_1a <- procD.lm(Structure ~ Taxon, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(slm_1a)

slm_1b <- procD.lm(Structure ~ 1, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(slm_1b)

anova(slm_1b, slm_1a)

slm_1_pw <- pairwise(fit = slm_1a, fit.null = slm_1b, groups = Neonympha_data.frame$Taxon) # much better fit

summary(slm_1_pw, test.type = "dist", confidence = 0.95) # Structire can distinguish taxa, pattern can too. Interesting that structure distinguished Nar from Nhe. Some latitude effect? 

# Structure by latitude
slm_2a <- procD.lm(Structure ~ Taxon * Lat, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(slm_2a)

slm_2b <- procD.lm(Structure ~ Taxon + Lat, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(slm_2b)

anova(slm_2b, slm_2a) # Interaction model is better

slm_2_pw <- pairwise(fit = slm_2a, fit.null = slm_2b, groups = Neonympha_data.frame$Taxon) # 

summary(slm_2_pw, test.type = "dist", confidence = 0.95)

### Covariate size by taxon
clm_1a <- procD.lm(Combined_data[, 87:90] ~ Taxon, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(clm_1a)

clm_1b <- procD.lm(Combined_data[, 87:90] ~ 1, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(clm_1b)

anova(clm_1b, clm_1a)

clm_1_pw <- pairwise(fit = clm_1a, fit.null = clm_1b, groups = Neonympha_data.frame$Taxon) # 

summary(clm_1_pw, test.type = "dist", confidence = 0.95)


### Covariate

#### Covariate length and width by taxon
# apd_4 <- advanced.procD.lm(f1 = Combined_data[, 95:102] ~ 1, Combined_data[, 95:102] ~ Taxon, data = Neonympha_data.frame, iter = 1e4, seed = 76234, print.progress = FALSE)
# summary(apd_4)

clm_2a <- procD.lm(Combined_data[, 95:102] ~ Taxon, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(clm_2a)

clm_2b <- procD.lm(Combined_data[, 95:102] ~ 1, data = Neonympha_data.frame, effect.type = "F", SS.type = "I", iter = 1e4, seed = 76234, print.progress = TRUE)
summary(clm_2b)

anova(clm_2b, clm_2a)

clm_2_pw <- pairwise(fit = clm_2a, fit.null = clm_2b, groups = Neonympha_data.frame$Taxon) # 

summary(clm_2_pw, test.type = "dist", confidence = 0.95)

# Misc ----

## Warp grids

### Pattern
# Need to do calculate the mean shape for each taxon
Pp <- dim(Pattern_gpa$coords)[1]
Pk <- dim(Pattern_gpa$coords)[2]
P_group <- Combined_data$Taxon
PY <- array(NA, dim = c(Pp, Pk, length(levels(P_group))))
dimnames(PY)[[3]] <- levels(P_group)

for(i in 1:length(levels(P_group))){
  grp <- Pattern_2d[which(P_group == levels(P_group)[i]), ]
  foo <- arrayspecs(grp, Pp, Pk)
  PY[, , i] <- mshape(foo)
}
dim(PY)

P_reference <- mshape(Pattern_gpa$coords) # Mean shape for Neonympha pattern

Nar_PY <- PY[, , 1]
Nhe_PY <- PY[, , 3]
Nfr_PY <- PY[, , 2]
Nmi_PY <- PY[, , 4]

plot(Pattern_pca_lda$LD1, Pattern_pca_lda$LD2, col = c(SW_palette("Main")[1], SW_palette("ROTJ")[3], SW_palette("Inquisitor")[1], SW_palette("TESB")[8])[Pattern_pca_lda$Taxon], pch = 19, xlab= 'LD1 (75.5%)', ylab ='LD 2 (15.5%)', xlim =c(-5, 7), ylim = c(-5, 5), cex = 2, las = 1 )

# Warp grid for N. areolata
par(fig = c(0.08, 0.28, 0.58, 0.93), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
N_ar_pars <- gridPar(tar.pt.bg = SW_palette("Main")[1], tar.pt.size = 1.5)
plotRefToTarget(P_reference, Nar_PY, method = "TPS", mag = 3, gridPars = N_ar_pars)

# Warp grid for N. helicta
par(fig = c(0.08, 0.28, 0.1, 0.45), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
N_he_pars <- gridPar(tar.pt.bg = SW_palette("Inquisitor")[1], tar.pt.size = 1.5)
plotRefToTarget(P_reference, Nhe_PY, method = "TPS", mag = 3, gridPars = N_he_pars)

# Warp grid for N. m. francisi
par(fig = c(0.76, 0.96, 0.1, 0.45), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
N_fr_pars <- gridPar(tar.pt.bg = SW_palette("ROTJ")[3], tar.pt.size = 1.5)
plotRefToTarget(P_reference, Nfr_PY, method = "TPS", mag = 3, gridPars = N_fr_pars)

# Warp grid for N. m. mitchellii
par(fig = c(0.76, 0.96, 0.58, 0.93), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
N_mi_pars <- gridPar(tar.pt.bg = SW_palette("TESB")[8], tar.pt.size = 1.5)
plotRefToTarget(P_reference, Nmi_PY, method = "TPS", mag = 3, gridPars = N_mi_pars)

dev.off()


### Structure
Sp <- dim(Structure_gpa$coords)[1]
Sk <- dim(Structure_gpa$coords)[2]
S_group <- Combined_data$Taxon
SY <- array(NA, dim = c(Sp, Sk, length(levels(S_group))))
dimnames(SY)[[3]] <- levels(S_group)

for(i in 1:length(levels(S_group))){
  grp <- Structure_2d[which(S_group == levels(S_group)[i]), ]
  foo <- arrayspecs(grp, Sp, Sk)
  SY[, , i] <- mshape(foo)
}
dim(SY)

S_reference <- mshape(Structure_gpa$coords) # Mean shape for Neonympha Structure

Nar_SY <- SY[, , 1]
Nhe_SY <- SY[, , 3]
Nfr_SY <- SY[, , 2]
Nmi_SY <- SY[, , 4]

plot(Structure_pca_lda$LD1, Structure_pca_lda$LD2, col = c(SW_palette("Main")[1], SW_palette("ROTJ")[3], SW_palette("Inquisitor")[1], SW_palette("TESB")[8])[Structure_pca_lda$Taxon], pch = 19, xlab= 'LD1 (75.5%)', ylab ='LD 2 (15.5%)', xlim =c(-5, 7), ylim = c(-5, 5), cex = 2, las = 1 )

# Warp grid for N. areolata
par(fig = c(0.08, 0.28, 0.55, 0.9), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
plotRefToTarget(S_reference, Nar_SY, method = "TPS", mag = 3, gridPars = N_ar_pars)

# Warp grid for N. helicta
par(fig = c(0.08, 0.28, 0.15, 0.5), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
plotRefToTarget(S_reference, Nhe_SY, method = "TPS", mag = 3, gridPars = N_he_pars)

# Warp grid for N. m. francisi
par(fig = c(0.76, 0.96, 0.15, 0.5), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
plotRefToTarget(S_reference, Nfr_SY, method = "TPS", mag = 3, gridPars = N_fr_pars)

# Warp grid for N. m. mitchellii
par(fig = c(0.76, 0.96, 0.55, 0.9), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
plotRefToTarget(S_reference, Nmi_SY, method = "TPS", mag = 3, gridPars = N_mi_pars)

dev.off()
