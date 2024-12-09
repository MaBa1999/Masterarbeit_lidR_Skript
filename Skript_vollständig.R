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

base::save.image("./Workspace/Daten_Clas_Ground.RData")

#__________________________________________________________________________________________________________________________________________#

### Rasterrize ###

base::rm(list = ls())
base::load("./Workspace/Daten_Clas_Ground.RData")

#Parallelprozessing Einstellen
cores <- parallelly::availableCores()
cores <- as.integer((cores*3)/7)
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
                                   algorithm = lidR::tin(extrapolate = lidR::knnidw(k = 10, p = 2, rmax = 50),
                                   use_class = c(2L, 9L)))

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
                               algorithm = lidR::tin(extrapolate =  knnidw(k = 10, p = 2, rmax = 50)),
                               use_class = c(2L, 9L))

# Überprüfen der Normalisierung
lidR::las_check(nctg)
lidR:::catalog_laxindex(nctg)

graphics::par(mfrow = c(5,2))
x <- 1:base::length(nctg@data$filename)
j <- base::sample(x, size = 10, replace = FALSE, prob = NULL)
for ( i in 1:10) {
  Test <- lidR::readLAS(nctg@data$filename[j[i]])  
  hist(filter_ground(Test)$Z, breaks = seq(-1, 1, 0.01), main = "", xlab = "Elevation")
}

for ( i in 1:10) {
  Test <- lidR::readLAS(nctg@data$filename[j[i]])
  Test <- lidR::filter_poi(Test, Classification == 2L)
  plot(Test, size = 3, bg = "white", color = "Classification")
}
Test <- lidR::readLAS(nctg@data$filename[500])
Test <- lidR::filter_poi(Test,
                         Classification == 2L)
plot(Test, size = 3, bg = "white", color = "Classification")

base::save.image("./Workspace/Normalisiert.RData")

#__________________________________________________________________________________________________________________________________________#

### Plot ###

base::rm(list = ls())
load("./Worksapce/Normalisiert.RData")

#Plots einlesen + Überprüfen
Plots_Vorrat <- utils::read.csv("./Daten/Rohdaten/Daten_FVA/vol_stp_092023.txt", header = TRUE, sep = ";")
utils::View(Plots_Vorrat)
base::View(Plots_Vorrat)
base::table(Plots_Vorrat$key)
base::table(Plots_Vorrat$kspnr)
base::table(Plots_Vorrat$vol_ha)
base::table(Plots_Vorrat$hoe_mod_mean)
base::table(base::is.na.data.frame(Plots_Vorrat))

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

base::rm("Plot_Vorrat_sf", "Transform_Plots_Vorrat", "Plots_Vorrat", "Ausdehnung_nctg")

base::save.image("./Workspace/Plots.RData")

#__________________________________________________________________________________________________________________________________________#

### Plot_Metriken ###

base::rm(list = ls())
base::load("./Workspace/Plots.RData")

graphics::par(mfrow = c(1,1))
plot(nctg)
plot(Plot_Vorrat_nctg, add = TRUE, col = "red")

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

opt_filter(nctg) <- "-drop_z_below 0" # Ignore points with elevations below 0

plot(nctg)
plot(Plot_Vorrat_nctg[1066,], add = TRUE, col = 'red')

Plot_Vorrat_nctg <- Plot_Vorrat_nctg[- 1066,]

D <- plot_metrics(nctg, .stdmetrics_z, Plot_Vorrat_nctg, radius = 13)

m <- lm(vol_ha ~ zq85, data = D)
summary(m)

plot(D$vol_ha, predict(m))
abline(0,1)

base::save.image("./Workspace/Plot_metriken.RData")

#__________________________________________________________________________________________________________________________________________#

### Hold-out ###

base::rm(list = ls())
base::load("./Workspace/Plot_metriken.RData")

#Betrachtung der Volumenverteilung
boxplot(Plot_Vorrat_nctg$vol_ha)
hist(Plot_Vorrat_nctg$vol_ha)
median(Plot_Vorrat_nctg$vol_ha)
mean(Plot_Vorrat_nctg$vol_ha)

#Aufteilen der Daten in drei Klassen
Einteilung <- max(Plot_Vorrat_nctg$vol_ha) - min(Plot_Vorrat_nctg$vol_ha)
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
cores <- as.integer(cores)
future::plan(future::multisession, workers = cores)
lidR::set_lidr_threads(cores)
data.table::setDTthreads(restore_after_fork = TRUE, percent = 100)

opt_filter(nctg) <- "-drop_z_below 0"

#Überprüfen auf Grenzplots
lidR::plot(nctg)
lidR::plot(Train_all, add = TRUE, col = 'red')

#Plotmetriken für Trainingsdaten
PlotMetriksTrain <- lidR::plot_metrics(nctg, .stdmetrics_z, Train_all, radius = 13)

#Plotmetriken für Testdaten
lidR::opt_output_files(nctg) <- "./Daten/Test_Plot_metriken/Plot_coordinate_{ID}_{XLEFT}_{YBOTTOM}"

#Überpüfen und entfernen von Grenzplts im Testdatensatz
lidR::plot(nctg)
lidR::plot(Test_all, add = TRUE, col = 'red')
Test_all <- Test_all[- 202,]
lidR::plot(nctg)
lidR::plot(Test_all[202,], add = TRUE, col = 'red')

#Plotmetriken für Testdatensatz
PlotMetriksTest <- lidR::plot_metrics(nctg, .stdmetrics_z, Test_all, radius = 13)

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
Model <- randomForest::randomForest(vol_ha ~ ., data = data, ntree = 1000, xtest = subset(xtest, select = -vol_ha), ytest = ytest)
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
### Model-Verbesserungen ###

## Verbesserzes RandomForesz
data <- sf::st_drop_geometry(PlotMetriksTrain)
xtest <- sf::st_drop_geometry(PlotMetriksTest)
ytest <- PlotMetriksTest$vol_ha
data <- subset(data, select = c(vol_ha, zmean, zsd, zskew, zkurt, zentropy, pzabovezmean, pzabove2, zq5, zq20, zq35, zq50, zq65, zq80, zq95, zpcum1, zpcum3, zpcum5, zpcum7, zpcum9))
xtest <- subset(xtest, select = c(vol_ha, zmean, zsd, zskew, zkurt, zentropy, pzabovezmean, pzabove2, zq5, zq20, zq35, zq50, zq65, zq80, zq95, zpcum1, zpcum3, zpcum5, zpcum7, zpcum9))
Model <- randomForest::randomForest(vol_ha ~ ., data = data, ntree = 1000, xtest = subset(xtest, select = -vol_ha), ytest = ytest)
Model
randomForest::varImpPlot(Model)

#Überpfrüfen der Daten im Random_forrest
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
#data <- rbind(PlotMetriksTest, PlotMetriksTrain)

#einfaches Model mit allen Variablen
train_control <- caret::trainControl(method = "cv", number = 10)
rf_model1 <- caret::train(vol_ha ~ ., data = PlotMetriksTrain, method = "rf", ntree = 1000, maximize = TRUE, trControl = train_control)
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
rf_model2 <- caret::train(vol_ha ~ ., data = data_train, method = "rf", ntree = 1000, maximize = TRUE, trControl = train_control)
rf_model2; plot(caret::varImp(rf_model2))
rf_model2$finalModel$mse[length(rf_model2$finalModel$mse)]
rf_model2$finalModel$rsq[length(rf_model2$finalModel$rsq)] *100
predictions <- predict(rf_model2, newdata = PlotMetriksTest)
postResample(predictions, PlotMetriksTest$vol_ha)

## Erweitertes Aufräumen Testen
data <- PlotMetriksTrain
preProc <- preProcess(data, method = c("center", "scale"))
data_transformed <- predict(preProc, newdata = data)
preProc <- preProcess(subset(data_transformed, select = -vol_ha), method = "corr")
data_transformed <- predict(preProc, newdata = data_transformed)
train_control <- caret::trainControl(method = "cv", number = 10)
rf_model3 <- caret::train(vol_ha ~ ., data = data_transformed, method = "rf", ntree = 1000, maximize = TRUE, trControl = train_control)
rf_model3
rf_model3$finalModel$mse[length(rf_model3$finalModel$mse)]
rf_model3$finalModel$rsq[length(rf_model3$finalModel$rsq)] *100
data_test <- PlotMetriksTest
preProcT <- preProcess(data_test, method = c("center", "scale"))
data_test <- predict(preProcT, newdata = data_test)
predictions <- predict(rf_model3, newdata = data_test)
postResample(predictions, PlotMetriksTest$vol_ha)

base::save.image("./Workspace/RandomForest.RData")

#__________________________________________________________________________________________________________________________________________#

### W-to-W ###

base::rm(list = ls())
base::load("./Workspace/RandomForest.RData")
base::rm("data", "data_corr", "data_preprozess", "data_test", "data_train", "data_transformed", "preProc", "preProcT", "Einteilung", "x", "sel")

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
opt_filter(nctg) <- "-drop_z_below 0" # Ignore points with elevations below 0

metrics_w2w <- pixel_metrics(nctg, .stdmetrics_z, res = 23, pkg = "terra")

terra::writeRaster(metrics_w2w, "./Daten/WalltoWall.tif", overwrite = TRUE, progress = TRUE)

plot(metrics_w2w$zq85)
plot(metrics_w2w$zmean)
plot(metrics_w2w$pzabovezmean)

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

lidR::plot(Pred)

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
