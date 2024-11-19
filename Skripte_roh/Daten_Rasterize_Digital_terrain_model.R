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

ctg4 <- readLAScatalog("C:/Masterarbeit_R/Daten/Classify/Ground")
las_check(ctg4)

opt_chunk_buffer(ctg4) <- 30
opt_chunk_size(ctg4) <- 400
opt_laz_compression(ctg4) <- TRUE
opt_restart(ctg4) <- 1

opt_output_files(ctg4) <- "C:/Masterarbeit_R/Daten/Daten_Rasterize_Digital_terrain_model/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

dtm_tin <- rasterize_terrain(ctg4, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin, bg = "white")
