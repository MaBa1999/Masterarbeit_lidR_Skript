### Masterarbeit ABA Lidar Sollingen ###
# Bearbeitung des Moduls "Masterarbeit"
# Dieses Skript soll mit Laserscanningdaten und teresstrischen Daten ein Model für Voratskarten des Solling ermöglichen.
# Autor: Marc B.
# Betreuer: Herr Prof. Dr. Magdon, Herr Starker
# Version 1.0

#__________________________________________________________________________________________________________________________________________#

### Vorbereitungen: ###

#Lade nötige R-Packages
required_packages <- c("lidR", #Lidar Daten bearbeiten
                       "sf", #Geodaten erstellen
                       "ggplot2", #Karten erstellen
                       "raster", #Geodaten
                       "viridis",
                       "mapview", #Karten erstellen
                       "future", #Multicore für lidR
                       "parallel", #Multicore
                       "parallelly", #Multicore
                       "RCSF", #Classify Bodenpunkte
                       "lasR", #Lidar Daten bearbeiten
                       "terra", #Geodaten
                       "progress", #Progressbar
                       "rminer",
                       "randomForest", #Modelling
                       "tidyr", #Tabelle Spalten zusammenführen
                       "dplyr",
                       "lidRviewer",
                       "caret"
                       )

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == 'lasR') {
      install.packages('lasR', repos = 'https://r-lidar.r-universe.dev', 'https://cloud.r-project.org')
      } else {
      if (pkg == 'lidRviewer') {
        install.packages('lidRviewer', repos = c('https://r-lidar.r-universe.dev'))
        } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  }
}

base::rm(list = ls())

## Arbeitsverzeichnis
#Hier das Arbeitsverzeichnis einstellen
base::setwd("C:/Masterarbeit_R")

base::save.image("./Workspace/Start.RData")

#__________________________________________________________________________________________________________________________________________#

### Daten laden ###

base::rm(list = ls())
base::load("./Workspace/Start.RData")

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

## Punkteindexierung erstellen ##
lidR:::catalog_laxindex(ctg)
lidR::las_check(ctg)

#Filtern der Überlappungen
lidR::opt_laz_compression(ctg) <- TRUE
lidR::opt_output_files(ctg) <- "./Daten/Zwischendaten/filter_doppelte_punkte/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::filter_duplicates(ctg)
ctg2 <- lidR::readLAScatalog("./Daten/Zwischendaten/filter_doppelte_punkte/")
lidR::las_check(ctg2)

#Index für schnelleres Berechnen
lidR:::catalog_laxindex(ctg2)
lidR::las_check(ctg2)

base::save.image("./Workspace/Daten_gefiltert.RData")
#__________________________________________________________________________________________________________________________________________#

### Daten Noise filtern ###

base::rm(list = ls())
base::load("./Workspace/Daten_gefiltert.RData")

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

base::save.image("./Workspace/Daten_Clas_Noise.RData")

#__________________________________________________________________________________________________________________________________________#

### Daten Ground classify ###

base::rm(list = ls())
base::load("./Workspace/Daten_Clas_Noise.RData")

#Parallelprozessing einstellen
cores <- parallelly::availableCores()
cores <- base::as.integer(cores)
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

base::save.image("./Workspace/Daten_Clas_Ground.RData")

#__________________________________________________________________________________________________________________________________________#

### Rasterrize ###

base::rm(list = ls())
base::load("./Workspace/Daten_Clas_Ground.RData")

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- base::as.integer((cores*3)/7)
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
                                   algorithm = lidR::tin(),
                                   use_class = c(2L, 9L))

terra::writeRaster(dtm_tin, filename = "./Daten/Daten_Rasterize_Digital_terrain_model/Terrainmodel/Terrainmodel.tif", overwrite=FALSE, progress = TRUE)

x <- "./Daten/Daten_Rasterize_Digital_terrain_model/Terrainmodel/Terrainmodel.tif"
dtm_tin2 <- terra::rast(x)

lidR::plot_dtm3d(dtm_tin2, bg = "white")

base::save.image("./Workspace/DTM.RData")


#__________________________________________________________________________________________________________________________________________#

### Höhennormalisierung ###

base::rm(list = ls())
base::load("./Workspace/DTM.RData")

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

#Berechnung der Normalisierung
nctg <- lidR::normalize_height(ctg3,
                               algorithm = lidR::tin(),
                               use_class = c(2L, 9L))

# Überprüfen der Normalisierung
lidR::las_check(nctg)
lidR:::catalog_laxindex(nctg)

graphics::par(mfrow = c(5,2), mar = c(2, 2, 2, 2))
x <- 1:base::length(nctg@data$filename)
j <- base::sample(x, size = 10, replace = FALSE, prob = NULL)
for ( i in 1:10) {
  Test <- lidR::readLAS(nctg@data$filename[j[i]])  
  graphics::hist(lidR::filter_ground(Test)$Z, breaks = seq(-1, 1, 0.01), main = "", xlab = "Elevation")
}

for ( i in 1:10) {
  Test <- lidR::readLAS(nctg@data$filename[j[i]])
  Test <- lidR::filter_poi(Test, Classification == 2L)
  lidR::plot(Test, size = 3, bg = "white", color = "Classification")
}
Test <- lidR::readLAS(nctg@data$filename[500])
Test <- lidR::filter_poi(Test,
                         Classification == 2L)
lidR::plot(Test, size = 3, bg = "white", color = "Classification")

graphics::par(mfrow = c(1,1), mar = c(5, 4, 4, 2))

rm("Test", "i", "j", "x")
base::save.image("./Workspace/Normalisiert.RData")


#__________________________________________________________________________________________________________________________________________#

### Canopy Height model ###

base::rm(list = ls())
load("./Workspace/Normalisiert.RData")

nctg <- lidR::readLAScatalog("./Daten/Daten_Hoehennormalisierung/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(nctg) <- 400
lidR::opt_chunk_size(nctg) <- 1500
lidR::opt_restart(nctg) <- 1
lidR::opt_laz_compression(nctg) <- TRUE
lidR::opt_output_files(nctg) <- "./Daten/CHM/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores/4)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

chm <- lidR::rasterize_canopy(nctg, res = 1, algorithm = lidR::dsmtin())

terra::writeRaster(chm, filename = "./Daten/CHM.tif", overwrite=FALSE, progress = TRUE)

x <- "./Daten/CHM.tif"
chm <- terra::rast(x)

lidR::plot_dtm3d(chm, bg = "white")
terra::plot(chm)

lidR::opt_output_files(nctg) <- "./Daten/CHM2/Chunks_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

chm2 <- lidR::rasterize_canopy(nctg, res = 1, algorithm = lidR::p2r(subcircle = 0, na.fill = lidR::knnidw(k = 10, p = 2, rmax = 50)))

terra::writeRaster(chm2, filename = "./Daten/CHM2.tif", overwrite=FALSE, progress = TRUE)

x <- "./Daten/CHM.tif"
chm2 <- terra::rast(x)

lidR::plot_dtm3d(chm2, bg = "white")
terra::plot(chm2)

base::save.image("./Workspace/chm.RData")

#__________________________________________________________________________________________________________________________________________#

### Plot ###

base::rm(list = ls())
load("./Workspace/Normalisiert.RData")

#Plots einlesen + Überprüfen
Plots_Vorrat <- utils::read.csv("./Daten/Rohdaten/Daten_FVA/vol_stp_092023.txt", header = TRUE, sep = ";")
utils::View(Plots_Vorrat)
base::table(Plots_Vorrat$key)
base::table(Plots_Vorrat$kspnr)
base::table(Plots_Vorrat$vol_ha)
base::table(Plots_Vorrat$hoe_mod_mean)
base::table(base::is.na.data.frame(Plots_Vorrat))
graphics::boxplot(Plots_Vorrat$vol_ha, ylab = "Volumen je Hektar (m³/ha)", main = "Boxplot NW-FVA Volumen je Messpunkt")
graphics::boxplot(Plots_Vorrat$hoe_mod_mean, ylab = "mittlere Höhe je Plot (m)", main = "Boxplot NW-FVA mittlere Höhe je Messpunkt")


#Tranformieren von Gaus-Krüger zu ...
Transform_Plots_Vorrat <- base::data.frame(lon = Plots_Vorrat$rw, lat = Plots_Vorrat$hw, vol_ha = Plots_Vorrat$vol_ha)
for (i in 1:base::length(Plots_Vorrat$hw)) {
  point <- sf::st_sfc(sf::st_point(x = c(Plots_Vorrat$rw[i], Plots_Vorrat$hw[i]), dim = XY), crs = 31467)
  cords <- sf::st_coordinates(sf::st_transform(point, src = 31467, crs = 25832))
  Transform_Plots_Vorrat$lon[i] <- cords[1]
  Transform_Plots_Vorrat$lat[i] <- cords[2]
  print(i)
  rm("cords", "point")
}

lidR::plot(Transform_Plots_Vorrat, add = TRUE, col = "red")

#Filtern der benötigten Datenpunkte
#Vektor für Catalog erstellen

sf::st_write(sf::st_as_sf(nctg, "sf"), "./Daten/Vektor/nctg.shp")
Ausdehnung_nctg <- sf::read_sf("./Daten/Vektor/nctg.shp")
sf::st_crs(Ausdehnung_nctg)
sf::st_crs(Ausdehnung_nctg) <- 25832

sf::st_write(sf::st_as_sf(Transform_Plots_Vorrat, coords = c("lon", "lat")), "./Daten/Vektor/Plot_Vorrat_sf.shp")
Plot_Vorrat_sf <- sf::read_sf("./Daten/Vektor/Plot_Vorrat_sf.shp")
sf::st_crs(Plot_Vorrat_sf)  <- 4326
sf::st_transform(Plot_Vorrat_sf, src = 4326, crs = 25832)
sf::st_crs(Plot_Vorrat_sf)
sf::st_crs(Plot_Vorrat_sf) <- 25832

Plot_Vorrat_nctg <- sf::st_filter(Plot_Vorrat_sf, Ausdehnung_nctg)

base::plot(Plot_Vorrat_nctg)

base::rm("Plot_Vorrat_sf", "Transform_Plots_Vorrat", "Plots_Vorrat", "Ausdehnung_nctg", "i")

base::save.image("./Workspace/Plots.RData")

#__________________________________________________________________________________________________________________________________________#

### Plot_Metriken ###

base::rm(list = ls())
base::load("./Workspace/Plots.RData")

graphics::par(mfrow = c(1,1))
lidR::plot(nctg)
lidR::plot(Plot_Vorrat_nctg, add = TRUE, col = "red")

nctg <- lidR::readLAScatalog("./Daten/Daten_Hoehennormalisierung/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(nctg) <- 50
lidR::opt_chunk_size(nctg) <- 530
lidR::opt_restart(nctg) <- 1
lidR::opt_laz_compression(nctg) <- TRUE
lidR::opt_output_files(nctg) <- "./Daten/Daten_plot_metrics/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

lidR::opt_filter(nctg) <- "-drop_z_below 0" # Ignore points with elevations below 0

lidR::plot(nctg)
lidR::plot(Plot_Vorrat_nctg[1066,], add = TRUE, col = 'red')

Plot_Vorrat_nctg <- Plot_Vorrat_nctg[- 1066,]

D <- lidR::plot_metrics(nctg, lidR::.stdmetrics_z, Plot_Vorrat_nctg, radius = 13)

m <- stats::lm(vol_ha ~ zq85, data = D)
base::summary(m)

plot(D$vol_ha, stats::predict(m), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

base::save.image("./Workspace/Plot_metriken.RData")

#__________________________________________________________________________________________________________________________________________#

### Hold-out ###

base::rm(list = ls())
base::load("./Workspace/Plot_metriken.RData")

#Betrachtung der Volumenverteilung
graphics::boxplot(Plot_Vorrat_nctg$vol_ha)
graphics::hist(Plot_Vorrat_nctg$vol_ha)
stats::median(Plot_Vorrat_nctg$vol_ha)
base::mean(Plot_Vorrat_nctg$vol_ha)

#Aufteilen der Daten in drei Klassen
Einteilung <- base::max(Plot_Vorrat_nctg$vol_ha) - base::min(Plot_Vorrat_nctg$vol_ha)
Einteilung <- Einteilung / 3

Klasse1 <- dplyr::filter(Plot_Vorrat_nctg, vol_ha < Einteilung) 
Klasse2 <- dplyr::filter(Plot_Vorrat_nctg, vol_ha > Einteilung & vol_ha < (Einteilung * 2)) 
Klasse3 <- dplyr::filter(Plot_Vorrat_nctg, vol_ha > (Einteilung * 2)) 

#Jeweilige Klasse eigenes Holdout und zusammenfügen
holdout_Klasse1 <- rminer::holdout(Klasse1$vol_ha, ratio = .8, mode = 'random')
holdout_Klasse2 <- rminer::holdout(Klasse2$vol_ha, ratio = .8, mode = 'random')
holdout_Klasse3 <- rminer::holdout(Klasse3$vol_ha, ratio = .8, mode = 'random')

Train1 <- Klasse1[holdout_Klasse1$tr,]
Test1 <- Klasse1[holdout_Klasse1$ts,]

Train2 <- Klasse2[holdout_Klasse2$tr,]
Test2 <- Klasse2[holdout_Klasse2$ts,]

Train3 <- Klasse3[holdout_Klasse3$tr,]
Test3 <- Klasse3[holdout_Klasse3$ts,]

Train_all <- base::rbind(Train1, Train2, Train3)
Test_all <- base::rbind(Test1, Test2, Test3)

#laden der Daten für die Plotmetriken
nctg <- lidR::readLAScatalog("./Daten/Daten_Hoehennormalisierung/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(nctg) <- 50
lidR::opt_chunk_size(nctg) <- 530
lidR::opt_restart(nctg) <- 1
lidR::opt_laz_compression(nctg) <- TRUE
lidR::opt_output_files(nctg) <- "./Daten/Train_Plot_metriken/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- base::as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

lidR::opt_filter(nctg) <- "-drop_z_below 0"

#Überprüfen auf Grenzplots
lidR::plot(nctg)
lidR::plot(Train_all, add = TRUE, col = 'red')

#Plotmetriken für Trainingsdaten
PlotMetriksTrain <- lidR::plot_metrics(nctg, lidR::.stdmetrics, Train_all, radius = 13)

#Plotmetriken für Testdaten
lidR::opt_output_files(nctg) <- "./Daten/Test_Plot_metriken/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Überpüfen und entfernen von Grenzplts im Testdatensatz
lidR::plot(nctg)
lidR::plot(Test_all, add = TRUE, col = 'red')
Test_all <- Test_all[- 202,]
lidR::plot(nctg)
lidR::plot(Test_all[202,], add = TRUE, col = 'red')

#Plotmetriken für Testdatensatz
PlotMetriksTest <- lidR::plot_metrics(nctg, lidR::.stdmetrics, Test_all, radius = 13)

base::save.image("./Workspace/Holdout.RData")

#__________________________________________________________________________________________________________________________________________#

### Random Forest ###
base::rm(list = ls())
load("C:/Masterarbeit_R/Workspace/Holdout.RData")

#Überprüfen und entfernen von leeren Werten
table(is.na.data.frame(PlotMetriksTest))
table(is.na.data.frame(PlotMetriksTrain))
PlotMetriksTest <- stats::na.omit(PlotMetriksTest)
PlotMetriksTrain <- stats::na.omit(PlotMetriksTrain)

#Entfernen von Geometry
data <- sf::st_drop_geometry(PlotMetriksTrain)
xtest <- sf::st_drop_geometry(PlotMetriksTest)
ytest <- PlotMetriksTest$vol_ha

# Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

#random Forest mit allen Variablen
Model <- randomForest::randomForest(vol_ha ~ ., data = data, xtest = subset(xtest, select = -vol_ha), ytest = ytest)
Model

#Testen des Random Forests
Model$mse[length(Model$mse)]
Model$rsq[length(Model$rsq)] * 100
Model$test$mse[length(Model$test$mse)]
Model$test$rsq[length(Model$test$rsq)] * 100

# Überpfrüfen des RandomForest
caret::varImp(Model) == Model$importance
randomForest::varImpPlot(Model)


#___________________________________________________________#

### Kreuzvalidierung ###

#laden der Daten
load("C:/Masterarbeit_R/Workspace/Holdout.RData")
base::rm("Train1", "Train2", "Train3", "Train_all", "Test1", "Test2", "Test3", "Test_all", "Klasse1", "Klasse2", "Klasse3", "holdout_Klasse1", "holdout_Klasse2", "holdout_Klasse3")
table(is.na.data.frame(PlotMetriksTest))
table(is.na.data.frame(PlotMetriksTrain))
PlotMetriksTest <- stats::na.omit(PlotMetriksTest)
PlotMetriksTrain <- stats::na.omit(PlotMetriksTrain)
PlotMetriksTest <- sf::st_drop_geometry(PlotMetriksTest)
PlotMetriksTrain <- sf::st_drop_geometry(PlotMetriksTrain)

#einfaches Model mit allen Variablen
train_control <- caret::trainControl(method = "cv", number = 10)
rf_model1 <- caret::train(vol_ha ~ ., data = PlotMetriksTrain, method = "rf", maximize = TRUE, trControl = train_control)
rf_model1; plot(caret::varImp(rf_model1))

rf_model1$finalModel$mse[length(rf_model1$finalModel$mse)]
rf_model1$finalModel$rsq[length(rf_model1$finalModel$rsq)] *100
predictions <- predict(rf_model1, newdata = PlotMetriksTest)
postResample(predictions, PlotMetriksTest$vol_ha)

## Aufbereiten des kreuzvalidierten random forest
#Aufräumen der Daten in hinscicht auf Korellation
data_corr <- subset(PlotMetriksTrain, select = -vol_ha)
data_preprozess <- caret::preProcess(data_corr, method = "corr", na.remove = TRUE)
data_preprozess
sel <- c(data_preprozess$method$remove)
data_train <- subset(PlotMetriksTrain, select = -which(names(PlotMetriksTrain) %in% sel))
train_control <- caret::trainControl(method = "cv", number = 10)
rf_model2 <- caret::train(vol_ha ~ ., data = data_train, method = "rf", maximize = TRUE, trControl = train_control)
rf_model2; plot(caret::varImp(rf_model2))
rf_model2$finalModel$mse[length(rf_model2$finalModel$mse)]
rf_model2$finalModel$rsq[length(rf_model2$finalModel$rsq)] *100
predictions <- predict(rf_model2, newdata = PlotMetriksTest)
postResample(predictions, PlotMetriksTest$vol_ha)

## Erweitertes Aufräumen Testen
#preprocess Train
data <- subset(PlotMetriksTrain, select = -vol_ha)
preProc <- caret::preProcess(data, method = c("center", "scale", "corr"), na.remove = TRUE)
preProc
data_transformed <- stats::predict(preProc, newdata = data)
data_transformed <- base::data.frame(data_transformed, vol_ha = PlotMetriksTrain$vol_ha)

#Modell 3
train_control <- caret::trainControl(method = "cv", number = 10)
rf_model3 <- caret::train(vol_ha ~ ., data = data_transformed, method = "rf", maximize = TRUE, trContol = train_control)
rf_model3; plot(caret::varImp(rf_model3))
rf_model3$finalModel$mse[length(rf_model3$finalModel$mse)]
rf_model3$finalModel$rsq[length(rf_model3$finalModel$rsq)] *100

#Test Modell 3
data_test <- PlotMetriksTest
preProcT <- subset(PlotMetriksTest, select = -vol_ha)
preProcT <- caret::preProcess(preProcT, method = c("center", "scale"), na.remove = TRUE)
data_test <- stats::predict(preProcT, newdata = data_test)
data_test <- data.frame(data_test, vol_ha = PlotMetriksTest$vol_ha)
predictions <- stats::predict(rf_model3, newdata = data_test)
caret::postResample(predictions, PlotMetriksTest$vol_ha)


#Streudiagramme
#Train
plot(x = PlotMetriksTrain$vol_ha, y = predict(Model, newdata = PlotMetriksTrain), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

plot(x = PlotMetriksTrain$vol_ha, y = predict(rf_model1, newdata = PlotMetriksTrain), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

plot(x = PlotMetriksTrain$vol_ha, y = predict(rf_model2, newdata = PlotMetriksTrain), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

plot(x = PlotMetriksTrain$vol_ha, y = predict(rf_model3, newdata = PlotMetriksTrain), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

#Test
plot(x = PlotMetriksTest$vol_ha, y = predict(Model, newdata = PlotMetriksTest), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

plot(x = PlotMetriksTest$vol_ha, y = predict(rf_model1, newdata = PlotMetriksTest), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

plot(x = PlotMetriksTest$vol_ha, y = predict(rf_model2, newdata = PlotMetriksTest), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)

plot(x = PlotMetriksTest$vol_ha, y = predict(rf_model3, newdata = PlotMetriksTest), xlab = "Gemessenes Volumen (m³/ha)", ylab = "geschätztes Volumen (m³/ha)")
graphics::abline(0,1)



base::save.image("./Workspace/RandomForest.RData")



#__________________________________________________________________________________________________________________________________________#

### W-to-W ###

base::rm(list = ls())
base::load("./Workspace/RandomForest.RData")

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores/2)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

nctg <- lidR::readLAScatalog("./Daten/Daten_Hoehennormalisierung/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(nctg) <- 30
lidR::opt_chunk_size(nctg) <- 530
lidR::opt_restart(nctg) <- 1
lidR::opt_laz_compression(nctg) <- TRUE
lidR::opt_output_files(nctg) <- "./Daten/Wall_to_wall/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"
lidR::opt_filter(nctg) <- "-drop_z_below 0" # Ignore points with elevations below 0

metrics_w2w <- lidR::pixel_metrics(nctg, lidR::.stdmetrics, res = 23, pkg = "terra")

terra::writeRaster(metrics_w2w, "./Daten/WalltoWall.tif", overwrite = TRUE, progress = TRUE)

metrics_w2w_100 <- lidR::pixel_metrics(nctg, lidR::.stdmetrics, res = 100, pkg = "terra")
terra::writeRaster(metrics_w2w_100, "./Daten/WalltoWall_100.tif", overwrite = TRUE, progress = TRUE)

Gitter <-  terra::rast("./Gitter/RT_T32UNC_A035741_20240109T103324_B03.tif")
metrics_w2w_sentinel <- lidR::pixel_metrics(nctg, lidR::.stdmetrics, res = Gitter, pkg = "terra")
terra::writeRaster(metrics_w2w_sentinel, "./Daten/WalltoWall_sentinel.tif", overwrite = TRUE, progress = TRUE)



base::save.image("./Workspace/WtoW.RData")

#__________________________________________________________________________________________________________________________________________#

### Karte ###

base::rm(list = ls())
base::load("./Workspace/WtoW.RData")

metrics_w2w <- terra::rast("./Daten/WalltoWall.tif")
#metrics_w2w <- as.data.frame(metrics_w2w)
#metrics_w2w <- sf::st_drop_geometry(metrics_w2w)
#metrics_w2w <- stats::na.omit(metrics_w2w)

Pred <- raster::predict(metrics_w2w, model = rf_model1, filname = "./Ergebnisse/Daten/Vorratskarte_geschaetzt.tif", progress = "text", na.rm = TRUE)
terra::writeRaster(Pred, "./Ergebnisse/Daten/Vorratskarte_geschaetzt.tif", overwrite = TRUE, progress = TRUE)

lidR::plot(Pred, col = viridis::viridis(1000), main = "Holzvorratskarte (Volumen je Hektar)", box = TRUE)

Pred2 <- raster::predict(metrics_w2w, model = rf_model2, filname = "./Ergebnisse/Daten/Vorratskarte2_geschaetzt.tif", progress = "text", na.rm = TRUE)
terra::writeRaster(Pred2, "./Ergebnisse/Daten/Vorratskarte2_geschaetzt.tif", overwrite = TRUE, progress = TRUE)

lidR::plot(Pred2, col = viridis::viridis(1000))

metrics_w2w_100 <- terra::rast("./Daten/WalltoWall_100.tif")

Pred <- raster::predict(metrics_w2w_100, model = rf_model1, filname = "./Ergebnisse/Daten/Vorratskarte_geschaetzt_100.tif", progress = "text", na.rm = TRUE)
terra::writeRaster(Pred, "./Ergebnisse/Daten/Vorratskarte_geschaetzt_100.tif", overwrite = TRUE, progress = TRUE)

lidR::plot(Pred, col = viridis::viridis(1000), main = "Holzvorratskarte (Volumen je Hektar)", box = TRUE)



Pred <- raster::predict(metrics_w2w_sentinel, model = rf_model1, filname = "./Ergebnisse/Daten/Vorratskarte_geschaetzt_sentinel.tif", progress = "text", na.rm = TRUE)
terra::writeRaster(Pred, "./Ergebnisse/Daten/Vorratskarte_geschaetzt_sentinel.tif", overwrite = TRUE, progress = TRUE)

lidR::plot(Pred, col = viridis::viridis(1000), main = "Holzvorratskarte (Volumen je Hektar)", box = TRUE)



base::save.image("./Workspace/Karte.RData")


#__________________________________________________#

### Karte Validieren ###

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

Pred <- terra::rast("./Ergebnisse/Daten/Vorratskarte_geschaetzt.tif")

graphics::par(mfrow = c(1,2))

Verteilung_Preq <- terra::as.data.frame(Pred)
graphics::hist(Verteilung_Preq$lyr1)

Plots_Vorrat <- utils::read.csv("./Daten/Rohdaten/Daten_FVA/vol_stp_092023.txt", header = TRUE, sep = ";")
Plots_Vorrat <- subset(Plots_Vorrat, select = vol_ha)
graphics::hist(Plots_Vorrat$vol_ha)

base::mean(Verteilung_Preq$lyr1)
base::mean(Plots_Vorrat$vol_ha)  
graphics::boxplot(Verteilung_Preq$lyr1)
graphics::boxplot(Plots_Vorrat$vol_ha)

Kombiniert <- data.frame(Vol_ha = Verteilung_Preq$lyr1)
for (i in 1:base::length(Plots_Vorrat$vol_ha)){
  Kombiniert[i,2] <- Plots_Vorrat$vol_ha[i]
}
table(is.na.data.frame(Kombiniert))

ggplot2::ggplot(data = Kombiniert) +
  ggplot2::geom_histogram(ggplot2::aes(x = Vol_ha, fill = "Verteilung_Preq"), binwidth = 10, alpha = 0.5, color = "orange") +
  ggplot2::geom_histogram(ggplot2::aes(x = V2, fill = "Plots_Vorrat"), binwidth = 10, alpha = 0.5, color = "blue") +
  ggplot2::scale_fill_manual(values = c("Verteilung_Preq" = "orange", "Plots_Vorrat" = "blue"), name = "Kombiniert") +
  ggplot2::labs(title = "Häufigkeitsverteilung Vergleich", x = "Volumen pro Hektar", y = "Häufigkeit") +
  ggplot2::theme_minimal()

ggplot2::ggplot(data = Kombiniert) +
  ggplot2::geom_density(ggplot2::aes(x = Vol_ha, color = "Verteilung_Preq"), size = 1) +
  ggplot2::geom_density(ggplot2::aes(x = V2, color = "Plots_Vorrat"), size = 1) +
  ggplot2::labs(title = "Dichteverteilungen Vergleich", x = "Volumen pro Hektar", y = "Dichte") +
  ggplot2::scale_color_manual(values = c("Verteilung_Preq" = "orange", "Plots_Vorrat" = "blue")) +
  ggplot2::theme_minimal()

ggplot2::ggplot(data = Kombiniert) +
  ggplot2::geom_boxplot(ggplot2::aes(x = Vol_ha, fill = "Verteilung_Preq"), binwidth = 10, alpha = 0.5, color = "orange") +
  ggplot2::geom_boxplot(ggplot2::aes(x = V2, fill = "Plots_Vorrat"), binwidth = 10, alpha = 0.5, color = "blue") +
  ggplot2::scale_fill_manual(values = c("Verteilung_Preq" = "orange", "Plots_Vorrat" = "blue"), name = "Kombiniert") +
  ggplot2::labs(title = "Häufigkeitsverteilung Vergleich", x = "Volumen pro Hektar", y = "Häufigkeit") +
  ggplot2::theme_minimal()

Pred <- terra::rast("./Ergebnisse/Daten/Vorratskarte_geschaetzt.tif")
graphics::par(mfrow = c(1,1))
lidR::plot(Pred, mapview = TRUE)
lidR::plot(nctg, add = TRUE)

for (i in 1:base::length(nctg@data$filename)) {
  las <- lidR::readLAS(nctg@data$filename[i])
  Zuschnitte <- terra::crop(Pred, las, snap = "near")
  terra::writeRaster(Zuschnitte, paste("./Ergebnisse/Daten/Vorratskarte_zugeschnitten_", tools::file_path_sans_ext(base::basename(nctg@data$filename[i])), "_", i, ".tif"))
  base::print(i)
}
lidR::plot(Zuschnitte)

las <- lidR::readLAS(nctg@data$filename[118])
lidR::plot(las)

x <- PlotMetriksTest$vol_ha
y <- predict(rf_model1, newdata = data_test)
plot(x, y, xlab = "Test: Vol/ha", ylab = "Geschätzt Model 1")
graphics::abline(0,1)


base::save.image("./Workspace/Karte_validiert.RData")

#__________________________________________________#

#### Verbesserungen durch entfernen der von potenziellen Ausreißer ### 
base::load("./Workspace/Karte_validiert.RData")

#Entfernen Ausreißer basierend auf Boxplot
data <- Plot_Vorrat_nctg

q1 <- stats::quantile(data$vol_ha, 0.25)
q3 <- stats::quantile(data$vol_ha, 0.75)
iqr <- q3 - q1

GrenzeU <- q1 - 1.5 * iqr
GrenzeO <- q3 + 1.5 * iqr

data_clean <- data[data$vol_ha >= GrenzeU & data$vol_ha <= GrenzeO,]
print(data_clean)


#Aufteilen der Daten in drei Klassen
Einteilung <- base::max(data_clean$vol_ha) - base::min(data_clean$vol_ha)
Einteilung <- Einteilung / 3

Klasse1 <- dplyr::filter(data_clean, vol_ha < Einteilung) 
Klasse2 <- dplyr::filter(data_clean, vol_ha > Einteilung & vol_ha < (Einteilung * 2)) 
Klasse3 <- dplyr::filter(data_clean, vol_ha > (Einteilung * 2)) 

#Jeweilige Klasse eigenes Holdout und zusammenfügen
holdout_Klasse1 <- rminer::holdout(Klasse1$vol_ha, ratio = .8, mode = 'random')
holdout_Klasse2 <- rminer::holdout(Klasse2$vol_ha, ratio = .8, mode = 'random')
holdout_Klasse3 <- rminer::holdout(Klasse3$vol_ha, ratio = .8, mode = 'random')

Train1 <- Klasse1[holdout_Klasse1$tr,]
Test1 <- Klasse1[holdout_Klasse1$ts,]

Train2 <- Klasse2[holdout_Klasse2$tr,]
Test2 <- Klasse2[holdout_Klasse2$ts,]

Train3 <- Klasse3[holdout_Klasse3$tr,]
Test3 <- Klasse3[holdout_Klasse3$ts,]

Train_all <- base::rbind(Train1, Train2, Train3)
Test_all <- base::rbind(Test1, Test2, Test3)



#laden der Daten für die Plotmetriken
nctg <- lidR::readLAScatalog("./Daten/Daten_Hoehennormalisierung/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(nctg) <- 50
lidR::opt_chunk_size(nctg) <- 530
lidR::opt_restart(nctg) <- 1
lidR::opt_laz_compression(nctg) <- TRUE
lidR::opt_output_files(nctg) <- "./Daten/VerbesserungenPlotMetriken/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- base::as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

lidR::opt_filter(nctg) <- "-drop_z_below 0"

#Überprüfen auf Grenzplots
lidR::plot(nctg)
lidR::plot(Train_all, add = TRUE, col = 'red')

#Plotmetriken für Trainingsdaten
PlotMetriksTrain <- lidR::plot_metrics(nctg, lidR::.stdmetrics, Train_all, radius = 13)

#Plotmetriken für Testdaten
lidR::opt_output_files(nctg) <- "./Daten/VerbesserungenPlotTestMetriken/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Überpüfen und entfernen von Grenzplts im Testdatensatz
lidR::plot(nctg)
lidR::plot(Test_all, add = TRUE, col = 'red')
Test_all <- Test_all[- 202,]
lidR::plot(nctg)
lidR::plot(Test_all[202,], add = TRUE, col = 'red')

#Plotmetriken für Testdatensatz
PlotMetriksTest <- lidR::plot_metrics(nctg, lidR::.stdmetrics, Test_all, radius = 13)


#Überprüfen und entfernen von leeren Werten
table(is.na.data.frame(PlotMetriksTest))
table(is.na.data.frame(PlotMetriksTrain))
PlotMetriksTest <- stats::na.omit(PlotMetriksTest)
PlotMetriksTrain <- stats::na.omit(PlotMetriksTrain)

#Entfernen von Geometry
data <- sf::st_drop_geometry(PlotMetriksTrain)
xtest <- sf::st_drop_geometry(PlotMetriksTest)
ytest <- PlotMetriksTest$vol_ha

# Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)

#random Forest mit allen Variablen
Model <- randomForest::randomForest(vol_ha ~ ., data = data, xtest = subset(xtest, select = -vol_ha), ytest = ytest)
Model

#Testen des Random Forests
Model$mse[length(Model$mse)]
Model$rsq[length(Model$rsq)] * 100
Model$test$mse[length(Model$test$mse)]
Model$test$rsq[length(Model$test$rsq)] * 100

# Überpfrüfen des RandomForest
caret::varImp(Model) == Model$importance
randomForest::varImpPlot(Model)

## Erweitertes Aufräumen Testen
#preprocess Train
data <- sf::st_drop_geometry(subset(PlotMetriksTrain, select = -vol_ha))
preProc <- caret::preProcess(data, method = c("center", "scale", "corr"), na.remove = TRUE)
preProc
data_transformed <- stats::predict(preProc, newdata = data)
data_transformed <- base::data.frame(data_transformed, vol_ha = PlotMetriksTrain$vol_ha)

#Modell 3
train_control <- caret::trainControl(method = "cv", number = 10)
rf_model3 <- caret::train(vol_ha ~ ., data = data_transformed, method = "rf", maximize = TRUE, trContol = train_control)
rf_model3; plot(caret::varImp(rf_model3))
rf_model3$finalModel$mse[length(rf_model3$finalModel$mse)]
rf_model3$finalModel$rsq[length(rf_model3$finalModel$rsq)] *100

#Test Modell 3
data_test <- sf::st_drop_geometry(PlotMetriksTest)
preProcT <- subset(data_test, select = -vol_ha)
preProcT <- caret::preProcess(preProcT, method = c("center", "scale"), na.remove = TRUE)
data_test <- stats::predict(preProcT, newdata = data_test)
data_test <- data.frame(data_test, vol_ha = PlotMetriksTest$vol_ha)
predictions <- stats::predict(rf_model3, newdata = data_test)
caret::postResample(predictions, PlotMetriksTest$vol_ha)


base::save.image("./Workspace/VerbessertVersuch.RData")


#__________________________________________________#

### Versuch der Extrakton von Kronengrößen ####

base::load("./workspace/VerbessertVersuch.RData")

#laden der Daten 
nctg <- lidR::readLAScatalog("./Daten/Daten_Hoehennormalisierung/", filter = "-drop_class 18")
lidR::opt_chunk_buffer(nctg) <- 50
lidR::opt_chunk_size(nctg) <- 600
lidR::opt_restart(nctg) <- 1
lidR::opt_laz_compression(nctg) <- TRUE
lidR::opt_output_files(nctg) <- "./Daten/indiTree/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- base::as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

lidR::opt_filter(nctg) <- "-drop_z_below 0"

#laden des Canopy-height-model
x <- "./Daten/CHM.tif"
chm <- terra::rast(x)

#Berechnung der Baumpositionen
ttops <- lidR::locate_trees(chm, lidR::lmf(ws = 5))

#Hochsetzen der maximalen Dateigrößen im parallelprozessing
options(future.globals.maxSize = 800 * 1024^2)

# Extraktion einzelner Bäume
TreeCrowns <- lidR::segment_trees(nctg ,lidR::silva2016(chm, ttops, ID = "treeID")) #Fehlerhaft, es werden keine Bäume extrahiert

