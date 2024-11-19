### Classify ###

#Required packages
#lidR
if(!require(lidR)){install.packages("lidR")}
library("lidR")

if(!require(sf)){install.packages("sf")}
library("sf")

if(!require(ggplot2)){install.packages("ggplot2")}
library("ggplot2")

if(!require(raster)){install.packages("raster")}
library("raster")

if(!require(viridis)){install.packages("viridis")}
library("viridis")

if(!require(mapview)){install.packages("mapview")}
library("mapview")

if(!require(future)){install.packages("future")}
library("future")

if(!require(parallel)){install.packages("parallel")}
library("parallel")

if(!require(parallelly)){install.packages("parallelly")}
library("parallelly")

#usage of more cores (CPU)
cores <- availableCores()
cores <- as.integer(cores)
plan(multisession, workers = as.integer(cores))
set_lidr_threads(as.integer(cores))

las_check(ctg2)

opt_chunk_buffer(ctg2) <- 30
opt_chunk_size(ctg2) <- 400
opt_laz_compression(ctg2) <- TRUE
opt_restart(ctg2) <- 1

opt_output_files(ctg2) <- "C:/Masterarbeit_R/Daten/Classify/Noise/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

### COde_test ###
las <- readLAS(ctg2@data$filename[1])
plot(las)
las <- classify_noise(las, algorithm = sor())
plot(las, color = "Classification")
p1 <- c(mean(las$X), max(las$Y))
p2 <- c(mean(las$X), min(las$Y))
las_tr <- clip_transect(las, p1, p2, width = 5, xz = TRUE)
ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
las_denoise <- filter_poi(las, Classification != LASNOISE)
las_denoise <- classify_noise(las_denoise, algorithm = ivf())
plot(las_denoise, color = "Classification")
las_denoise <- filter_poi(las_denoise, Classification != LASNOISE)
las_denoise <- lidR::classify_ground(las_denoise, algorithm = csf())
plot(las_denoise, color = "Classification")
plot(filter_ground(las_denoise))
p1 <- c(mean(las_denoise$X), max(las_denoise$Y))
p2 <- c(mean(las_denoise$X), min(las_denoise$Y))
las_tr <- clip_transect(las_denoise, p1, p2, width = 5, xz = TRUE)
ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))


### COde_test ###
las <- readLAS(ctg2@data$filename[1000])
#plot(las)
las <- classify_noise(las, algorithm = sor())
#plot(las, color = "Classification")
p1 <- c(mean(las$X), max(las$Y))
p2 <- c(mean(las$X), min(las$Y))
las_tr <- clip_transect(las, p1, p2, width = 5, xz = TRUE)
ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
las_denoise <- filter_poi(las, Classification != LASNOISE)
las_denoise <- classify_noise(las_denoise, algorithm = ivf())
#plot(las_denoise, color = "Classification")
las_denoise <- filter_poi(las_denoise, Classification != LASNOISE)
las_denoise <- lidR::classify_ground(las_denoise, algorithm = csf())
#plot(las_denoise, color = "Classification")
#plot(filter_ground(las_denoise))
p1 <- c(mean(las_denoise$X), max(las_denoise$Y))
p2 <- c(mean(las_denoise$X), min(las_denoise$Y))
las_tr <- clip_transect(las_denoise, p1, p2, width = 5, xz = TRUE)
ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
plot(las)


### Noise entfernen ###
for (i in 1:length(ctg2$filename)) {
  las_denoise <- readLAS(ctg2@data$filename[i])
  las_denoise <- classify_noise(ctg2, algorithm = sor())
  las_denoise <- filter_poi(las, Classification != LASNOISE)
  las_denoise <- classify_noise(las_denoise, algorithm = ivf())
  las_denoise <- filter_poi(las_denoise, Classification != LASNOISE)
  # Extrahiere nur den Dateinamen 
  filename <- basename(ctg$filename[i]) 
  # Entferne die Dateiendung ".laz"
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- paste0("C:/Masterarbeit_R/Daten/Classify/Noise", filename, "_denoise.laz")
  writeLAS(las_denoise, neu_filename, index=TRUE)
  print(i)
}


### LasCatalog ohne Noise ###
ctg3 <- readLAScatalog("C:/Masterarbeit_R/Daten/Classify/Noise")
las_check(ctg3)
lidR:::catalog_laxindex(ctg3)
las_check(ctg3)

