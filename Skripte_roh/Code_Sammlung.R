ctg <- readLAScatalog("D://ALS-Daten_solling/04_classified")

opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- "D://Berechnung/las_lax/{ORIGINALFILENAME}_converted"

options(lidR.progress = list(cores = 4))


#for (i in 1:length(ctg$filename)) {
  #laz <- readLAS(ctg$filename[i])
  #filename <- ctg$filename[i]
  #neu_filename <- paste0("D://Berechnung/las_lax/", filename, "_converted.las")
  #writeLAS(laz, neu_filename, index=TRUE)
#}


for (i in 1:length(ctg$filename)) {
  laz <- readLAS(ctg$filename[i])
  # Extrahiere nur den Dateinamen 
  filename <- basename(ctg$filename[i]) 
  # Entferne die Dateiendung ".laz"
  filename <- tools::file_path_sans_ext(filename)
  neu_filename <- paste0("D://Berechnung/las_lax/", filename, "_converted.las")
  writeLAS(laz, neu_filename, index=TRUE)
}


#Laz to LAS
#opt_chunk_buffer(ctg) <- 0
#opt_chunk_size(ctg) <- 0
#opt_output_files(ctg) <- "D://Berechnung/las_lax/{ORIGINALFILENAME}_converted"

#for (i in 1:length(ctg$filename)) {
#laz <- readLAS(ctg$filename[i])
# Extrahiere nur den Dateinamen 
#filename <- basename(ctg$filename[i]) 
# Entferne die Dateiendung ".laz"
#filename <- tools::file_path_sans_ext(filename)
#neu_filename <- paste0("D://Berechnung/las_lax/", filename, "_converted.las")
#writeLAS(laz, neu_filename, index=TRUE)
#}

#ctg <- readLAScatalog("D://Berechnung/las_lax/")
#print(ctg)
#las_check(ctg)
#plot(ctg, mapview = TRUE)


#opt_restart()
