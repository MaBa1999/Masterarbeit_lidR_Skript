### Masterarbeit ABA Lidar Sollingen ###
# Bearbeitung des Moduls "Masterarbeit"
# Dieses Skript soll mit Laserscanningdaten und teresstrischen Daten ein Model für Voratskarten des Solling ermöglichen.
# Autor: Marc B.
# Betreuer: Herr Prof. Dr. Magdon, Herr Starker
# Version 1.0

#__________________________________________________________________________________________________________________________________________#

### Vorbereitungen: ###

#Lade nötige R-Packages
required_packages <- c("lidR", "sf", "ggplot2", "raster", "viridis", "mapview", "future", "parallel", "parallelly", "RCSF", "lasR")

for (pkg in required_packages) {
  if (!require(pkg, character.only = T)) {
    if (pkg == 'lasR') {
      install.packages('lasR', repos = 'https://r-lidar.r-universe.dev')
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = T)
  }
}


## Parallel-Prozessing 
cores <- parallelly::availableCores()
#Hier je nach Leistung des PCs einstellen, hier werden nur 1/4 der Threads der CPU genutzt:
cores <- as.integer(cores/4)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)


## Arbeitsverzeichnis
#Hier das Arbeitsverzeichnis einstellen
base::setwd("C:/Masterarbeit_R")

#__________________________________________________________________________________________________________________________________________#

### Daten laden ###

#LidR catalog
ctg <- lidR::readLAScatalog("./Daten/Rohdaten/ALS-Daten_solling/04_classified")

## Daten Überprüfen

base::print(ctg)
lidR::las_check(ctg)
lidR::plot(ctg, mapview = TRUE)

## Daten Überlappung filtern

#Einstellung der Chunksgröße und des Buffers
lidR::opt_chunk_buffer(ctg) <- 200
lidR::opt_chunk_size(ctg) <- 400

#Filtern der Überlappungen
lidR::opt_laz_compression(ctg) <- TRUE
lidR::opt_output_files(ctg) <- "./Daten/Zwischendaten/filter_doppelte_punkte/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::filter_duplicates(ctg)
ctg2 <- lidR::readLAScatalog("./Daten/Zwischendaten/filter_doppelte_punkte")
lidR::las_check(ctg2)

#Index für schnelleres Berechnen
lidR:::catalog_laxindex(ctg2)
lidR::las_check(ctg2)

#__________________________________________________________________________________________________________________________________________#

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

#Test
ctg_test <- lidR::readLAScatalog(c(ctg2@data$filename[1], ctg2@data$filename[2], ctg2@data$filename[3], ctg2@data$filename[4], ctg2@data$filename[5]))
lidR::las_check(ctg_test)
lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_test) <- 30
lidR::classify_noise(ctg_test, algorithm = sor())
ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/sor/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::las_check(ctg_test)
lidR::opt_chunk_buffer(ctg_test) <- 0
for (i in 1:length(ctg_test@data$filename)) {
  test <- lidR::readLAS(ctg_test@data$filename[i])
  test <-  lidR::filter_poi(test, Classification != LASNOISE)
  filename <- base::basename(ctg$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/sor/", filename, ".laz")
  lidR::writeLAS(test, neu_filename, index=TRUE)
  base::print(i)
}
ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::las_check(ctg_test)

ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/sor/")
lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/ivf/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_chunk_buffer(ctg_test) <- 30
lidR::classify_noise(ctg_test, algorithm = ivf())
ctg_test <- lidR::readLAScatalog("./Daten/Classify/Noise/ivf/")
lidR::las_check(ctg_test)

lidR::opt_output_files(ctg_test) <- "./Daten/Classify/Noise/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

ctg_test <- lidR::filter_poi(ctg_test, Classification != LASNOISE)


for ( i in 1:base::length(base::list.files(path = "./Daten/Classify/Noise/sor/"))) {
  las_denoise <- lidR::readLAS(ctg2@data$filename[i])
  las_denoise <- lidR::filter_poi(ctg2@data$filename[i], Classification != LASNOISE)
  filename <- base::basename(ctg$filename[i])
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise/sor/", filename, "_denoise.laz")
  lidR::writeLAS(las_denoise, neu_filename, index=TRUE)
  base::print(i)
}
lidR::opt_output_files(ctg2) <- "./Daten/Classify/Noise/ivf/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::classify_noise(ctg2, algorithm = ivf())
ctg2 <- lidR::readLAScatalog("./Daten/Classify/Noise/ivf/")
for ( i in 1:base::length(base::list.files(path = "./Daten/Classify/Noise/ivf/"))) {
  ctg2@data$filename[i] <- lidR::filter_poi(ctg2@data$filename[i], Classification != LASNOISE)
}



## LasCatalog ohne Noise erstellen
ctg3 <- lidR::readLAScatalog("./Daten/Classify/Noise")
lidR::las_check(ctg3)
lidR:::catalog_laxindex(ctg3)
lidR::las_check(ctg3)

#__________________________________________________________________________________________________________________________________________#

