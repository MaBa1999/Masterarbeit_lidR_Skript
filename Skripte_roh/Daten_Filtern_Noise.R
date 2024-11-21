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






### Daten Noise filtern ###

lidR::las_check(ctg2)

#Einstellen der Größe der Chunks und des Buffers + Startpunkt der Berechnung + Speicherart/-ortes
lidR::opt_chunk_buffer(ctg2) <- 30
lidR::opt_chunk_size(ctg2) <- 400
lidR::opt_restart(ctg2) <- 1
lidR::opt_laz_compression(ctg2) <- TRUE
lidR::opt_output_files(ctg2) <- "./Daten/Classify/Noise/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

## Noise entfernen gesamtes Gebiet
#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

##Exploitative Datenanalyse
# Code_test
las <- lidR::readLAS(ctg2@data$filename[1])
lidR::plot(las)
las <- lidR::classify_noise(las, algorithm = sor())
lidR::plot(las, color = "Classification")
p1 <- c(base::mean(las$X), base::max(las$Y))
p2 <- c(base::mean(las$X), base::min(las$Y))
las_tr <- lidR::clip_transect(las, p1, p2, width = 5, xz = TRUE)
ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
las_denoise <- lidR::filter_poi(las, Classification != LASNOISE)
las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
lidR::plot(las_denoise, color = "Classification")
las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
las_denoise <- lidR::classify_ground(las_denoise, algorithm = csf())
lidR::plot(las_denoise, color = "Classification")
lidR::plot(lidR::filter_ground(las_denoise))
p1 <- c(base::mean(las_denoise$X), base::max(las_denoise$Y))
p2 <- c(base::mean(las_denoise$X), base::min(las_denoise$Y))
las_tr <- lidR::clip_transect(las_denoise, p1, p2, width = 5, xz = TRUE)
ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))

# COde_test 2
las <- lidR::readLAS(ctg2@data$filename[1000])
#lidR::plot(las)
las <- lidR::classify_noise(las, algorithm = sor())
#lidR::plot(las, color = "Classification")
p1 <- c(base::mean(las$X), base::max(las$Y))
p2 <- c(base::mean(las$X), base::min(las$Y))
las_tr <- lidR::clip_transect(las, p1, p2, width = 5, xz = TRUE)
ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
las_denoise <- lidR::filter_poi(las, Classification != LASNOISE)
las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
#lidR::plot(las_denoise, color = "Classification")
las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
las_denoise <- lidR::classify_ground(las_denoise, algorithm = csf())
#lidR::plot(las_denoise, color = "Classification")
#lidR::plot(lidR::filter_ground(las_denoise))
p1 <- c(base::mean(las_denoise$X), base::max(las_denoise$Y))
p2 <- c(base::mean(las_denoise$X), base::min(las_denoise$Y))
las_tr <- lidR::clip_transect(las_denoise, p1, p2, width = 5, xz = TRUE)
ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
lidR::plot(las)

#Schlife Filterung Noise
for (i in 1:base::length(ctg2$filename)) {
  las_denoise <- lidR::readLAS(ctg2@data$filename[i])
  las_denoise <- lidR::classify_noise(las_denoise, algorithm = sor())
  las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
  las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
  las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
  # Extrahiere nur den Dateinamen 
  filename <- base::basename(ctg$filename[i])
  # Entferne die Dateiendung ".laz"
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/Schleife/", filename, "_denoise.laz")
  lidR::writeLAS(las_denoise, neu_filename, index=TRUE)
  base::print(i)
}

### Test
ctg_test <- lidR::readLAScatalog(c(ctg2@data$filename[1], ctg2@data$filename[2], ctg2@data$filename[3], ctg2@data$filename[4], ctg2@data$filename[5]))
lidR::opt_laz_compression(ctg_test) <- TRUE
lidR::las_check(ctg_test)
lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_test) <- 30
lidR::classify_noise(ctg_test, algorithm = sor())
ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_test) <- TRUE
lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::las_check(ctg_test)
lidR::opt_chunk_buffer(ctg_test) <- 0
for (i in 1:length(ctg_test@data$filename)) {
  test <- lidR::readLAS(ctg_test@data$filename[i])
  test <-  lidR::filter_poi(test, Classification != LASNOISE)
  filename <- base::basename(ctg_test$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/sor/", filename, ".laz")
  lidR::writeLAS(test, neu_filename, index=TRUE)
  base::print(i)
}
ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_test) <- TRUE
lidR::las_check(ctg_test)

ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_test) <- TRUE
lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/ivf/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_test) <- 30
lidR::classify_noise(ctg_test, algorithm = ivf())
ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/ivf/")
lidR::opt_laz_compression(ctg_test) <- TRUE
lidR::las_check(ctg_test)
lidR::opt_chunk_buffer(ctg_test) <- 0
for (i in 1:length(ctg_test@data$filename)) {
  test <- lidR::readLAS(ctg_test@data$filename[i])
  test <-  lidR::filter_poi(test, Classification != LASNOISE)
  filename <- base::basename(ctg_test$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/ivf/", filename, ".laz")
  lidR::writeLAS(test, neu_filename, index=TRUE)
  base::print(i)
}

ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/ivf/")
lidR::opt_laz_compression(ctg_test) <- TRUE
lidR::las_check(ctg_test)


## LasCatalog ohne Noise erstellen
ctg3 <- lidR::readLAScatalog("./Daten/Classify/Noise")
lidR::las_check(ctg3)
lidR:::catalog_laxindex(ctg3)
lidR::las_check(ctg3)





## Noise - Filterung <- teilweise multicore
#Filterung mit sor()
ctg_noise <- ctg2
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::las_check(ctg_noise)
lidR::opt_output_files(ctg_noise) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_noise) <- 30
lidR::classify_noise(ctg_noise, algorithm = sor())
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::opt_output_files(ctg_noise) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::las_check(ctg_noise)
lidR::opt_chunk_buffer(ctg_noise) <- 0
#hier werden die Daten Überschrieben
for (i in 1:length(ctg_noise@data$filename)) {
  Noise <- lidR::readLAS(ctg_noise@data$filename[i])
  Noise <-  lidR::filter_poi(Noise, Classification != LASNOISE)
  filename <- base::basename(ctg_noise$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/sor/", filename, ".laz")
  lidR::writeLAS(Noise, neu_filename, index=FALSE)
  base::print(i)
}
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR:::catalog_laxindex(ctg2)
lidR::las_check(ctg_noise)

#Filterung mit ivf()
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::opt_output_files(ctg_noise) <- "./Daten/Classify/Noise/ivf/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_noise) <- 30
lidR::classify_noise(ctg_noise, algorithm = ivf())
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/ivf/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::las_check(ctg_noise)
lidR::opt_chunk_buffer(ctg_noise) <- 0
#hier werden die Daten Überschrieben
for (i in 1:length(ctg_noise@data$filename)) {
  Noise <- lidR::readLAS(ctg_noise@data$filename[i])
  Noise <-  lidR::filter_poi(Noise, Classification != LASNOISE)
  filename <- base::basename(ctg_noise$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/ivf/", filename, ".laz")
  lidR::writeLAS(Noise, neu_filename, index=FALSE)
  base::print(i)
}

#Schleife Filterung & Löschen: Noise mit sor() & ivf() <- nicht Multicore
#for (i in 1:base::length(ctg2$filename)) {
  #las_denoise <- lidR::readLAS(ctg2@data$filename[i])
  #las_denoise <- lidR::classify_noise(las_denoise, algorithm = sor())
  #las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
  #las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
  #las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
  ## Extrahiere nur den Dateinamen 
  #filename <- base::basename(ctg$filename[i])
  ## Entferne die Dateiendung ".laz"
  #filename <- tools::file_path_sans_ext(filename)
  #neu_filename <- base::paste0("./Daten/Classify/Noise/Schleife/", filename, "_denoise.laz")
  #lidR::writeLAS(las_denoise, neu_filename, index=TRUE)
  #base::print(i)
#}

## Noise - Filterung <- teilweise multicore
#Filterung mit sor()
ctg_noise <- ctg2
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::las_check(ctg_noise)
lidR::opt_output_files(ctg_noise) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_noise) <- 30
lidR::classify_noise(ctg_noise, algorithm = sor())
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::opt_output_files(ctg_noise) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::las_check(ctg_noise)
lidR::opt_chunk_buffer(ctg_noise) <- 0
#hier werden die Daten Überschrieben
for (i in 1:length(ctg_noise@data$filename)) {
  Noise <- lidR::readLAS(ctg_noise@data$filename[i])
  Noise <-  lidR::filter_poi(Noise, Classification != LASNOISE)
  filename <- base::basename(ctg_noise$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/sor/", filename, ".laz")
  lidR::writeLAS(Noise, neu_filename, index=FALSE)
  base::print(i)
}
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR:::catalog_laxindex(ctg2)
lidR::las_check(ctg_noise)

#Filterung mit ivf()
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::opt_output_files(ctg_noise) <- "./Daten/Classify/Noise/ivf/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_noise) <- 30
lidR::classify_noise(ctg_noise, algorithm = ivf())
ctg_noise <- lidR::readLAScatalog("./Daten/Classify/Noise/ivf/")
lidR::opt_laz_compression(ctg_noise) <- TRUE
lidR::las_check(ctg_noise)
lidR::opt_chunk_buffer(ctg_noise) <- 0
#hier werden die Daten Überschrieben
for (i in 1:length(ctg_noise@data$filename)) {
  Noise <- lidR::readLAS(ctg_noise@data$filename[i])
  Noise <-  lidR::filter_poi(Noise, Classification != LASNOISE)
  filename <- base::basename(ctg_noise$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/ivf/", filename, ".laz")
  lidR::writeLAS(Noise, neu_filename, index=FALSE)
  base::print(i)
}

##Exploitative Datenanalyse
# Code_test
# las <- lidR::readLAS(ctg2@data$filename[1])
# lidR::plot(las)
# las <- lidR::classify_noise(las, algorithm = sor())
# lidR::plot(las, color = "Classification")
# p1 <- c(base::mean(las$X), base::max(las$Y))
# p2 <- c(base::mean(las$X), base::min(las$Y))
# las_tr <- lidR::clip_transect(las, p1, p2, width = 5, xz = TRUE)
# ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
#   geom_point(size = 0.5) +
#   coord_equal() +
#   theme_minimal() +
#   scale_color_gradientn(colours = height.colors(50))
# las_denoise <- lidR::filter_poi(las, Classification != LASNOISE)
# las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
# lidR::plot(las_denoise, color = "Classification")
# las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
# las_denoise <- lidR::classify_ground(las_denoise, algorithm = csf())
# lidR::plot(las_denoise, color = "Classification")
# lidR::plot(lidR::filter_ground(las_denoise))
# p1 <- c(base::mean(las_denoise$X), base::max(las_denoise$Y))
# p2 <- c(base::mean(las_denoise$X), base::min(las_denoise$Y))
# las_tr <- lidR::clip_transect(las_denoise, p1, p2, width = 5, xz = TRUE)
# ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
#   geom_point(size = 0.5) +
#   coord_equal() +
#   theme_minimal() +
#   scale_color_gradientn(colours = height.colors(50))

# COde_test 2
# las <- lidR::readLAS(ctg2@data$filename[1000])
# #lidR::plot(las)
# las <- lidR::classify_noise(las, algorithm = sor())
# #lidR::plot(las, color = "Classification")
# p1 <- c(base::mean(las$X), base::max(las$Y))
# p2 <- c(base::mean(las$X), base::min(las$Y))
# las_tr <- lidR::clip_transect(las, p1, p2, width = 5, xz = TRUE)
# ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
#   geom_point(size = 0.5) +
#   coord_equal() +
#   theme_minimal() +
#   scale_color_gradientn(colours = height.colors(50))
# las_denoise <- lidR::filter_poi(las, Classification != LASNOISE)
# las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
# #lidR::plot(las_denoise, color = "Classification")
# las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
# las_denoise <- lidR::classify_ground(las_denoise, algorithm = csf())
# #lidR::plot(las_denoise, color = "Classification")
# #lidR::plot(lidR::filter_ground(las_denoise))
# p1 <- c(base::mean(las_denoise$X), base::max(las_denoise$Y))
# p2 <- c(base::mean(las_denoise$X), base::min(las_denoise$Y))
# las_tr <- lidR::clip_transect(las_denoise, p1, p2, width = 5, xz = TRUE)
# ggplot2::ggplot(payload(las_tr), aes(X,Z, color = Classification)) +
#   geom_point(size = 0.5) +
#   coord_equal() +
#   theme_minimal() +
#   scale_color_gradientn(colours = height.colors(50))
# lidR::plot(las)
