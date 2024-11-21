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
lidR::opt_chunk_buffer(ctg2) <- 30
lidR::opt_chunk_size(ctg2) <- 400
lidR::opt_restart(ctg2) <- 1
lidR::opt_laz_compression(ctg2) <- TRUE
lidR::opt_output_files(ctg2) <- "./Daten/Classify/Noise/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Filterung mit ivf()
ctg2 <- lidR::classify_noise(ctg2, algorithm = ivf())
ctg3 <- lidR::readLAScatalog("./Daten/Classify/Noise/", filter = "-drop_class 18")

#Überprüfen der Daten 
lidR::opt_laz_compression(ctg3) <- TRUE
lidR::las_check(ctg3)

#__________________________________________________________________________________________________________________________________________#

