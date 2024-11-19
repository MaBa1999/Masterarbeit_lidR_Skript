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

if(!require(RCSF)){install.packages("RCSF")}
library("RCSF")

#usage of more cores (CPU)
cores <- availableCores()
cores <- as.integer(cores)
plan(multisession, workers = as.integer(cores))
set_lidr_threads(as.integer(cores))

ctg3 <- readLAScatalog("C:/Masterarbeit_R/Daten/Classify/Noise")
las_check(ctg3)

opt_chunk_buffer(ctg3) <- 30
opt_chunk_size(ctg3) <- 400
opt_laz_compression(ctg3) <- TRUE
opt_restart(ctg3) <- 1

opt_output_files(ctg3) <- "C:/Masterarbeit_R/Daten/Classify/Ground/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

classify_ground(ctg3, csf(sloop_smooth = FALSE, class_threshold = 0.5, cloth_resolution = 2))

ctg4 <- readLAScatalog("C:/Masterarbeit_R/Daten/Classify/Ground")
lidR:::catalog_laxindex(ctg4)
las_check(ctg4)
