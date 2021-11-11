library(png)
source <- 'data/tumor/IDP 2020-02-27-002_mousemulti - cropping (1, 9451, 7028, 6223, 6386).png'
my_im <- readPNG(source, native = FALSE, info = FALSE)
library(geojsonR)
js_file_name <- "data/tumor/IDP 2020-02-27-002_mousemulti - cropping (1, 9451, 7028, 6223, 6386).geojson"
file_js = FROM_GeoJson(url_file_string = js_file_name)
plot(1:5,1:5)
library(geojsonio)
spdf <- geojson_read(js_file_name,  what = "sp")
library(sp)
par(mar=c(0,0,0,0))
plot(spdf, col="white")
spdf_fortified <- broom::tidy(spdf)
# Plot it
library(ggplot2)
#spdf_swap <- rbind(spdf_fortified[spdf_fortified$id == "3",],
#                   spdf_fortified[spdf_fortified$id == "2",],
#                   spdf_fortified[spdf_fortified$id == "1",])
#spdf_fortified$id <- factor(spdf_fortified$id, levels = c("3","2","1"))
spdf_fortified <- spdf_fortified[!(spdf_fortified$long>3000 & spdf_fortified$long<3800 & spdf_fortified$lat>3600 & spdf_fortified$lat<4100 & spdf_fortified$id == "3"),]
p <- ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = -lat, group = group, color = id), fill=NA)  + theme_classic()+ coord_fixed()
p 
plot_df <- data.frame(cbind(1:4,1:4))
colnames(plot_df) <- c('x','y')
p <- ggplot(plot_df, aes(x,y)) + geom_point()
p
