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


## Noise entfernen gesamtes Gebiet
for (i in 1:base::length(ctg2$filename)) {
  las_denoise <- lidR::readLAS(ctg2@data$filename[i])
  las_denoise <- lidR::classify_noise(ctg2, algorithm = sor())
  las_denoise <- lidR::filter_poi(las, Classification != LASNOISE)
  las_denoise <- lidR::classify_noise(las_denoise, algorithm = ivf())
  las_denoise <- lidR::filter_poi(las_denoise, Classification != LASNOISE)
  # Extrahiere nur den Dateinamen 
  filename <- base::basename(ctg$filename[i]) 
  # Entferne die Dateiendung ".laz"
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- base::paste0("./Daten/Classify/Noise", filename, "_denoise.laz")
  lidR::writeLAS(las_denoise, neu_filename, index=TRUE)
  base::print(i)
}

## LasCatalog ohne Noise erstellen
ctg3 <- lidR::readLAScatalog("./Daten/Classify/Noise")
lidR::las_check(ctg3)
lidR:::catalog_laxindex(ctg3)
lidR::las_check(ctg3)

#__________________________________________________________________________________________________________________________________________#

