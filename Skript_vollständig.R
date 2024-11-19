### Masterarbeit ABA Lidar Sollingen ###
# Bearbeitung des Moduls "Masterarbeit"
# Dieses Skript soll mit Laserscanningdaten und teresstrischen Daten ein Model für Voratskarten des Solling ermöglichen.
# Autor: Marc B.
# Betreuer: Herr Prof. Dr. Magdon, Herr Starker
# Version 1.0

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

#Parallel-Prozessing 
