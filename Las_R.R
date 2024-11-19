### Versuch mit lasR ###

#### Nötige Packages ###
if (!require(lidR)){install.packages('lidR'); library(lidR)}
if (!require(lasR)){install.packages('lasR', repos = 'https://r-lidar.r-universe.dev'); library(lasR)}

folder = "C:/Masterarbeit_R/Daten/Zwischendaten/filter_doppelte_punkte/"
pipeline = reader_las() + dtm()
ans <- exec(pipeline, on = folder, ncores = 24, progress = TRUE)
ans  

plot(ans)
plot_dtm3d(ans)
