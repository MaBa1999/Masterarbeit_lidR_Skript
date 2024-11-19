#Daten_Vorbereitung

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
as.integer(cores)
plan(multisession, workers = as.integer(cores))
set_lidr_threads(as.integer(cores))

#LidR catalog
ctg <- readLAScatalog("D://ALS-Daten_solling/04_classified")
print(ctg)
las_check(ctg)
plot(ctg, mapview = TRUE)

#plot(ctg, mapview = FALSE)
opt_chunk_size(ctg) <- 320

#Filtern von Überlappungen
#opt_output_files(ctg) <- paste0(tempdir(), "/tree_coordinate_{XLEFT}_{YBOTTOM}")
opt_laz_compression(ctg) <- TRUE
opt_output_files(ctg) <- "D://Berechnung/filter_doppelte_punkte/Chunks_coordinate_{XLEFT}_{YBOTTOM}_{ID}"
filter_duplicates(ctg)
ctg2 <- readLAScatalog("D://Berechnung/filter_doppelte_punkte")
las_check(ctg2)

#Index -> schnelleres Berechnen
lidR:::catalog_laxindex(ctg2)
las_check(ctg2)