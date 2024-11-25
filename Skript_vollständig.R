### Masterarbeit ABA Lidar Sollingen ###
# Bearbeitung des Moduls "Masterarbeit"
# Dieses Skript soll mit Laserscanningdaten und teresstrischen Daten ein Model für Voratskarten des Solling ermöglichen.
# Autor: Marc B.
# Betreuer: Herr Prof. Dr. Magdon, Herr Starker
# Version 1.0

#__________________________________________________________________________________________________________________________________________#

### Vorbereitungen: ###

#Lade nötige R-Packages
required_packages <- c("lidR", "sf", "ggplot2", "raster", "viridis", "mapview", "future", "parallel", "parallelly", "RCSF", "lasR", "terra")

for (pkg in required_packages) {
  if (!require(pkg, character.only = T)) {
    if (pkg == 'lasR') {
      install.packages('lasR', repos = 'https://r-lidar.r-universe.dev', 'https://cloud.r-project.org')
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = T)
  }
}

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

## Parallel-Prozessing 
cores <- parallelly::availableCores()
#Hier je nach Leistung des PCs einstellen, hier werden nur 1/4 der Threads der CPU genutzt:
cores <- as.integer(cores/4)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

#Filtern der Überlappungen
lidR::opt_laz_compression(ctg) <- TRUE
lidR::opt_output_files(ctg) <- "./Daten/Zwischendaten/filter_doppelte_punkte/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::filter_duplicates(ctg)
ctg2 <- lidR::readLAScatalog("./Daten/Zwischendaten/filter_doppelte_punkte/")
lidR::las_check(ctg2)

#Index für schnelleres Berechnen
lidR:::catalog_laxindex(ctg2)
lidR::las_check(ctg2)

#__________________________________________________________________________________________________________________________________________#

### Daten Noise filtern ###

ctg2 <- lidR::readLAScatalog("./Daten/Zwischendaten/filter_doppelte_punkte/")
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

#Filterung mit sor()
ctg2 <- lidR::classify_noise(ctg2, algorithm = sor())
ctg2 <- lidR::readLAScatalog("./Daten/Classify/Noise/", filter = "-drop_class 18")
lidR::las_check(ctg2)
lidR:::catalog_laxindex(ctg2)

ctg3 <- lidR::readLAScatalog("./Daten/Classify/Noise/", filter = "-drop_class 18")

#Überprüfen der Daten 
lidR::opt_laz_compression(ctg3) <- TRUE
lidR::las_check(ctg3)

#__________________________________________________________________________________________________________________________________________#

### Daten Ground classify ###

#Parallelprozessing einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

ctg3 <- lidR::readLAScatalog("./Daten/Classify/Noise/", filter = "-drop_class 18")
lidR::las_check(ctg3)

lidR::opt_chunk_buffer(ctg3) <- 30
lidR::opt_chunk_size(ctg3) <- 550
lidR::opt_restart(ctg3) <- 1
lidR::opt_laz_compression(ctg3) <- TRUE
lidR::opt_output_files(ctg3) <- "./Daten/Classify/Ground/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

ctg3 <- lidR::classify_ground(ctg3, lidR::csf(sloop_smooth = FALSE,
                                        class_threshold = 0.5,
                                        cloth_resolution = 2))

#__________________________________________________________________________________________________________________________________________#

### Rasterrize ###

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores/4)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

ctg3 <- lidR::readLAScatalog("./Daten/Classify/Ground/", filter = "-drop_class 18")
#lidR:::catalog_laxindex(ctg3)
lidR::las_check(ctg3)

lidR::opt_chunk_buffer(ctg3) <- 50
lidR::opt_chunk_size(ctg3) <- 610
lidR::opt_restart(ctg3) <- 1
lidR::opt_laz_compression(ctg3) <- TRUE
lidR::opt_output_files(ctg3) <- "./Daten/Daten_Rasterize_Digital_terrain_model/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

dtm_tin <- lidR::rasterize_terrain(ctg3,
                                   res = 1,
                                   algorithm = lidR::tin(extrapolate =  knnidw(k = 10, p = 2, rmax = 50),
                                   use_class = c(2L, 9L)))

lidR::plot_dtm3d(dtm_tin, bg = "white")

#__________________________________________________________________________________________________________________________________________#

### Höhennormalisierung ###

ctg3 <- lidR::readLAScatalog("./Daten/Classify/Ground/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(ctg3) <- 50
lidR::opt_chunk_size(ctg3) <- 610
lidR::opt_restart(ctg3) <- 1
lidR::opt_laz_compression(ctg3) <- TRUE
lidR::opt_output_files(ctg3) <- "./Daten/Daten_Hoehennormalisierung/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores/3)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

nctg <- lidR::normalize_height(ctg3,
                               algorithm = lidR::tin(extrapolate =  knnidw(k = 10, p = 2, rmax = 50)),
                               use_class = c(2L, 9L))

# Überprüfen der Normalisierung
lidR::las_check(nctg)
lidR:::catalog_laxindex(nctg)

Test <- lidR::readLAS(nctg@data$filename[500])
Test <- lidR::filter_poi(Test,
                         Classification == 2L)
plot(Test, size = 3, bg = "white", color = "Classification")

Test <- lidR::readLAS(nctg@data$filename[1000])
hist(filter_ground(Test)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")


#__________________________________________________________________________________________________________________________________________#

### Plot ###

