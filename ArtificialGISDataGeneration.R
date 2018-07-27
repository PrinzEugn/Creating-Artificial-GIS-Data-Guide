## Script for taking in dxf files exported from Inkscape and using them to 
## generate different map visualizations.
##  
##  Input files are the following:
##   - muscatatuck.png (base image)
##   - PointsTest.dxf (point data made in Inkscape)
##   - PointTest2.dxf (point data made in Inkscape)
##   - BordersTest.dxf (polygon border data made in Inkscape)
##   
##  Output files are the following:
##   -- Images --
##       - points_map.png
##       - points_2_map.png
##       - points_3_map.png
##       - borders_map.png
##       - features_map.png (reference)
##       - heatmap_alpha.png
##       - heatmap_bw.png
##       - heatmap_2_alpha.png
##       - heatmap_2_bw.png
##       - heatmap_3_alpha.png
##       - heatmap_3_bw.png
##       - points_1_choropleth.png
##       - points_2_choropleth.png
##       - points_3_choropleth.png
##       - Sector_IDs_Counts.png (reference)
##   -- Tables --
##       - points_1.csv
##       - points_2.csv
##       - points_3.csv
##       - Sector_Counts.csv
##       
## 
## Author: Mark Simpson, Pennsylvania State University
## Email: marksimpson@psu.edu
## Created June-July 2018

# Setup ------------------------------------------------------------
setwd("~/GitHub/Creating-Artificial-GIS-Data-Guide")

# Clear workspace
rm(list = ls())

#### Load packages ####

library(dplyr)
library(ggplot2)
#install.packages("dplyr", "ggplot")

library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(GISTools)

# needed to convert inkscape-created dxf lines to polygon
#install.packages("sf")
library(sf)

# remote sensing tools, for ggplot integration
#install.packages("RStoolbox")
library(RStoolbox)

# for theme_nothing
#install.packages("cowplot")
library(cowplot)

# Alternative plotting
#install.packages("png")
library(png)
library(grid)

# Base image loading and plotting -----------------------------------------

# Load base map
filename <- "muscatatuck.png"

# load as "RasterBrick" with multiple bands 
image <- brick(filename)

# convert to spatialgriddataframe for later processing
image.grid <- as(image, 'SpatialGridDataFrame')

# plot SpatialGridDataFrame first band, black and white
plot(image.grid, col = grey(seq(0, 1, length = 256)))
 
# plot with default plotting
#plotRGB(image)

# plot image to test
ggRGB(image) +
    theme(axis.title = element_blank()) +
    theme( panel.background = element_rect(fill = NA),
        panel.ontop = TRUE )

# NOTE: resolution is 1778*1452

# plot image with reference grid
ggRGB(image) +
    theme(axis.title = element_blank()) +
    theme( panel.background = element_rect(fill = NA),
           panel.grid.major = element_line( color = "white"),
           panel.grid.minor = element_line(color = "gray60"),
           panel.ontop = TRUE ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                       expand = c(0,0))

#### Test alternate image plotting ####

# works but is a little slower

# img <-  readPNG(filename)
# image2 <- rasterGrob(img, interpolate = TRUE)
# 
# ggplot(points, aes(x = lons, y = lats)) +
#     annotation_custom(image2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#     geom_point(color = "Yellow")




# Read Vector Data ------------------------------------------------------------
#### Points set 1 #####
# Read dxf created in Inkscape, will import as lines
points.dxf <- readOGR(dsn = "PointsTest.dxf")

# use rgeos to get centroid of triangles made in Inkscape
points <-  gCentroid(points.dxf, byid = TRUE)

# should be SpatialPoints
class(points)

#### Points set 2 #####
# Read dxf created in Inkscape, will import as lines
points.2.dxf <- readOGR(dsn = "PointsTest2.dxf")

# use rgeos to get centroid of triangles made in Inkscape
points.2 <-  gCentroid(points.2.dxf, byid = TRUE)

# should be SpatialPoints
class(points.2)

plot(points.2)

#### Generate point set 3 ####
set.seed(8675309)
x <-  runif(150, 0, 1778)
y <-  runif(150, 0, 1458)

# bind into new dataframe
points.3 <- as.data.frame(cbind(x, y))

points.3

plot(points.3)

# convert to spatial points for kde raster
points.3.spat <- SpatialPoints(points.3)

#### Borders Reading + Processing ####

# Read dxf created in Inkscape
borders.raw <- readOGR(dsn = "BordersTest.dxf")

# Loads as SpatialLinesDataFrame
class(borders.raw)

## use sf package to convert lines from dxf into polygons

# convert to 'simple feature'
borders.sf <- st_as_sf(borders.raw) 

# convert lines to polygon
borders.polygon <- st_polygonize(borders.sf)

# convert back to spatial polygon
borders <- as(borders.polygon, "Spatial")

# Should be SpatialPolygonsDataFrame
class(borders)

#create id field for merging later, needs to be 0-11
borders@data$id <- 0:(nrow(borders@data) - 1)

# reset plotting device ( may not be necessary)
dev.off()

#test Plot
plot(borders, col = "slategray")


#### Plotting Test #### 
#test plot
plot(borders, col = "gray20")
plot(points, pch = 16, col = "red", add = TRUE)
plot(points.2, pch = 15, col = "blue", add = TRUE)
plot(points.3.spat, pch = 17, col = "green", add = TRUE)

# Write out plots for vector features ------------------------------------------
#### Write plot of points 1 ####
# Create writeable plot
# 
points.map <- ggplot() +
    
    # call points as dataframe
    geom_point(data = as.data.frame(points), 
               aes(x = x, y = y), 
               color = "black", 
               fill = "red",
               pch = 21,
               size = 8) +

    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

# note point size is much smaller when written out
points.map

## begin plotting chunk

# designate filename, resolution
png( filename = "points_map.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
points.map

# write plot
dev.off()

#### Write plot of points 2 ####
# Create writeable plot
points.2.map <- ggplot() +
    # call points as dataframe
    geom_point(data = as.data.frame(points.2), 
               aes(x = x, y = y), 
               color = "black", 
               fill = "blue",
               pch = 21,
               size = 8) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

# note point size is much smaller when written out
points.2.map

## begin plotting chunk

# designate filename, resolution
png( filename = "points_2_map.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
points.2.map

# write plot
dev.off()

#### Write plot of points 3 ####
# Create writeable plot
points.3.map <- ggplot() +
    # call points as dataframe
    geom_point(data = as.data.frame(points.3), 
               aes(x = x, y = y), 
               color = "black", 
               fill = "green3",
               pch = 21,
               size = 8) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

# note point size is much smaller when written out
points.3.map

## begin plotting chunk

# designate filename, resolution
png( filename = "points_3_map.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
points.3.map

# write plot
dev.off()

#### Write plot of borders ####

# convert to ggplot friendly format
borders.fort <- fortify(borders)

# Create writeable plot
# 
borders.map <- ggplot() +
    
    # borders
    geom_polygon(data = borders.fort, 
                 aes( x = long, 
                      y = lat, 
                      group = group), 
                 color = "white",
                 fill = NA,
                 alpha = .1,
                 size = 1) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

# test with black background for visibility
borders.map + theme(panel.background = element_rect(fill = "gray10"))

## begin plotting chunk

# designate filename, resolution
png( filename = "borders_map.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
borders.map

# write plot
dev.off()


#### Write plot of borders and points ####
# Create writeable plot
features.map <- ggplot() +
    
    # borders
    geom_polygon(data = borders.fort, 
                 aes( x = long, 
                      y = lat, 
                      group = group), 
                 color = "white",
                 fill = NA,
                 alpha = .1,
                 size = 1) +
    
    # call points 1 as dataframe
    geom_point(data = as.data.frame(points), 
               aes(x = x, y = y), 
               color = "black", 
               fill = "red",
               pch = 21,
               size = 8) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

features.map

## begin plotting chunk

# designate filename, resolution
png( filename = "features_map.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
features.map

# write to disk
dev.off()


# Heatmap 1- kernel density raster ---------------------------------------------
#### create and convert kd raster####

# Create kernel density raster, manually change bandwidth
point.dens <- kde.points(points, h = 250, lims = image.grid)

plot(point.dens)

# Should be SpatialPixelsDataFrame
class(point.dens)

# Convert to dataframe for graphing
point.dens.df <- as.data.frame(point.dens)


#### Plot heatmap with color and alpha ####
# Plot it out
heatmap.alpha <- ggplot() + 
    
    geom_tile(data = point.dens.df, 
              aes(x = Var1, y = Var2, fill = "red", alpha = kde)) + 
    
    theme_nothing() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

heatmap.alpha

## begin plotting chunk

# designate filename, resolution
png( filename = "heatmap_alpha.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452, 
     bg = "transparent")

# call previous plot
heatmap.alpha

# write plot
dev.off()
##
#### Plot black + white heatmap ####

# Plot it out
heatmap.bw <- ggplot() + 
    
    geom_tile(data = point.dens.df, 
              aes(x = Var1, y = Var2, fill = kde),
              alpha = 1) + 
    scale_fill_gradient(low = "black", high = "white") +
    theme_nothing() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = NULL, y = NULL) 


## begin plotting chunk

# designate filename, resolution
png( filename = "heatmap_bw.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452 )

# call previous plot
heatmap.bw

# write plot
dev.off()

##
# Heatmap 2- kernel density raster ---------------------------------------------
#### create and convert kd raster####

# Create kernel density raster, manually change bandwidth
points.2.dens <- kde.points(points.2, h = 600, lims = image.grid)

plot(points.2.dens)

# Should be SpatialPixelsDataFrame
class(points.2.dens)

# Convert to dataframe for graphing
points.2.dens.df <- as.data.frame(points.2.dens)

head(points.2.dens.df)

#### Plot heatmap with color and alpha ####
# Plot it out
heatmap.2.alpha <- ggplot() + 
    
    geom_tile(data = points.2.dens.df, 
              aes(x = Var1, y = Var2, alpha = kde),  fill = "blue") + 
    
    theme_nothing() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

heatmap.2.alpha

## begin plotting chunk

# designate filename, resolution
png( filename = "heatmap_2_alpha.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452, 
     bg = "transparent")

# call previous plot
heatmap.2.alpha

# write plot
dev.off()
##
#### Plot black + white heatmap ####

# Plot it out
heatmap.2.bw <- ggplot() + 
    
    geom_tile(data = points.2.dens.df, 
              aes(x = Var1, y = Var2, fill = kde),
              alpha = 1) + 
    scale_fill_gradient(low = "black", high = "white") +
    theme_nothing() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = NULL, y = NULL) 


## begin plotting chunk

# designate filename, resolution
png( filename = "heatmap_2_bw.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452 )

# call previous plot
heatmap.2.bw

# write plot
dev.off()

##
# Heatmap 3- kernel density raster ---------------------------------------------
#### create and convert kd raster####

# note this uses points.3.spat, not points.3

# Create kernel density raster, manually change bandwidth
points.3.dens <- kde.points(points.3.spat, h = 300, lims = image.grid)

plot(points.3.dens)

# Should be SpatialPixelsDataFrame
class(points.3.dens)

# Convert to dataframe for graphing
points.3.dens.df <- as.data.frame(points.3.dens)

head(points.3.dens.df)

#### Plot heatmap with color and alpha ####
# Plot it out
heatmap.3.alpha <- ggplot() + 
    
    geom_tile(data = points.3.dens.df, 
              aes(x = Var1, y = Var2, alpha = kde), fill = "green3") + 
    
    theme_nothing() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

heatmap.3.alpha

## begin plotting chunk

# designate filename, resolution
png( filename = "heatmap_3_alpha.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452, 
     bg = "transparent")

# call previous plot
heatmap.3.alpha

# write plot
dev.off()
##
#### Plot black + white heatmap ####

# Plot it out
heatmap.3.bw <- ggplot() + 
    
    geom_tile(data = points.3.dens.df, 
              aes(x = Var1, y = Var2, fill = kde),
              alpha = 1) + 
    scale_fill_gradient(low = "black", high = "white") +
    theme_nothing() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = NULL, y = NULL) 

## begin plotting chunk

# designate filename, resolution
png( filename = "heatmap_3_bw.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452 )

# call previous plot
heatmap.3.bw

# write plot
dev.off()

##
# Export Points Tables ---------------------------------------------------------
write.csv(as.data.frame(points), file = "points_1.csv", row.names = FALSE)

write.csv(as.data.frame(points.2), file = "points_2.csv", row.names = FALSE)

write.csv(as.data.frame(points.3), file = "points_3.csv", row.names = FALSE)


# Choropleth Map ----------------------------------------------------------
#### Borders preprocessing
#### Processing points to polygons ####

# calculate point density for choropleth map (# points/area)
point.1.chor <- poly.counts(points, borders)/poly.areas(borders)
point.2.chor <- poly.counts(points.2, borders)/poly.areas(borders)
point.3.chor <- poly.counts(points.3.spat, borders)/poly.areas(borders)

#add point density to data slot
borders@data$point.1 <- point.1.chor
borders@data$point.2 <- point.2.chor
borders@data$point.3 <- point.3.chor

# fortify using id to code areas
borders.fort.2 <- fortify(borders, region = "id")

# convert id to numeric for merging
borders.fort.2$id <- as.numeric(borders.fort.2$id)

# left join to re-add data into fortified borders
borders.df <- left_join(borders.fort.2, borders@data, by = "id")

borders.df$point.1
borders.df$point.2
borders.df$point.3

#### Choropleth for points 1 ####
point.1.choro.map <- ggplot() +
    
    # borders
    geom_polygon(data = borders.df, 
                 aes( x = long, 
                      y = lat, 
                      group = group,
                      fill = point.1), 
                 color = "black",
                 size = 1) +

    scale_fill_distiller(type = "seq",
                         palette = "Reds",
                         direction = 1) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

## begin plotting chunk

# designate filename, resolution
png( filename = "points_1_choropleth.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
point.1.choro.map

# write plot
dev.off()

#### Choropleth for points 2
point.2.choro.map <- ggplot() +
    
    # borders
    geom_polygon(data = borders.df, 
                 aes( x = long, 
                      y = lat, 
                      group = group,
                      fill = point.2), 
                 color = "black",
                 size = 1) +
    
    scale_fill_distiller(type = "seq",
                         palette = "Blues",
                         direction = 1) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

## begin plotting chunk

# designate filename, resolution
png( filename = "points_2_choropleth.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
point.2.choro.map

# write plot
dev.off()

#### Choropleth for points 3 ####
point.3.choro.map <- ggplot() +
    
    # borders
    geom_polygon(data = borders.df, 
                 aes( x = long, 
                      y = lat, 
                      group = group,
                      fill = point.3), 
                 color = "black",
                 size = 1) +
    
    scale_fill_distiller(type = "seq",
                         palette = "Greens",
                         direction = 1) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

## begin plotting chunk

# designate filename, resolution
png( filename = "points_3_choropleth.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
point.3.choro.map

# write plot
dev.off()


# Plot Labeled Polygons ####

## Calculate Point Totals #

# Calculate raw point totals
point.1.count <- poly.counts(points, borders)
point.2.count <- poly.counts(points.2, borders)
point.3.count <- poly.counts(points.3.spat, borders)

#add point counts to data slot
borders@data$point.1.count <- point.1.count
borders@data$point.2.count <- point.2.count
borders@data$point.3.count <- point.3.count

# get centroids
# use rgeos to get centroid of polygons made in Inkscape
b.centroids <-  gCentroid(borders, byid = TRUE)

# should be SpatialPoints
class(b.centroids)

# convert to df
b.centroids.df <- as.data.frame(b.centroids)

#add ID field
b.centroids.df$id <- 0:11

# Join point fields from borders
b.centroids.df <- left_join(b.centroids.df, borders@data, by = "id")

head(b.centroids.df)
colnames(b.centroids.df)

# purge non-relevant data
b.centroids.df <- b.centroids.df[ , c("x", "y", "id", "point.1", "point.2", "point.3", 
                                      "point.1.count", "point.2.count", "point.3.count"  )]
# Create plot
centroids.map <- ggplot() +
    
    # borders
    geom_polygon(data = borders.fort, 
                 aes( x = long, 
                      y = lat, 
                      group = group), 
                 color = "gray50",
                 fill = NA,
                 alpha = .1,
                 size = 2) +
    
    # call points dataframe
    # geom_point(data = b.centroids.df, 
    #            aes(x = x, y = y), 
    #            color = "slategray3", 
    #            fill = "yellow",
    #            pch = 19,
    #            size = 2) +
    
    # Add labels
    
    # Sector labels
    geom_text(data = b.centroids.df, 
              aes(label = id, x = x, y = y),
              size = 30,
              nudge_y = 10) +
    
    # Point 1 Labels
    geom_text(data = b.centroids.df, 
              aes(label = point.1.count, x = x, y = y),
              color = "Red",
              size = 15,
              nudge_x = -40, nudge_y = -50) +
    
    # Point 2 Labels
    geom_text(data = b.centroids.df, 
              aes(label = point.2.count, x = x, y = y),
              color = "Blue",
              size = 15,
              nudge_x = 0, nudge_y = -50) +
    
    # Point 3 Labels
    geom_text(data = b.centroids.df, 
              aes(label = point.3.count, x = x, y = y),
              color = "green4",
              size = 15,
              nudge_x = 40, nudge_y = -50) +
    
    # removes all theme stuff, theme_nothing from cowplot package
    theme_nothing() +
    scale_x_continuous(limits = c(0, 1778), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1452), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    
    #trying to get actual alpha
    theme(
        panel.background = element_rect(fill = "transparent") # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , panel.grid.major = element_blank() # get rid of major grid
        , panel.grid.minor = element_blank() # get rid of minor grid
        , legend.background = element_rect(fill = "transparent") # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

centroids.map


## begin plotting chunk

# designate filename, resolution
png( filename = "Sectors_IDs_Counts.png", 
     type = "cairo", 
     units = "px", 
     width = 1778, 
     height = 1452,
     bg = "transparent")

# call previous plot
centroids.map

# write plot
dev.off()

#### Write counts and density table ####
write.csv(b.centroids.df, file = "Sector_Counts.csv", row.names = FALSE)

#### GISTools version for comparison ####
# point.1.shades <- auto.shading(borders@data$point.1,  cols = brewer.pal(5,'Reds')) 
# 
# png( filename = "points_1_choropleth2.png", 
#      type = "cairo", 
#      units = "px", 
#      width = 1778, 
#      height = 1452,
#      bg = "whi")
# 
# choropleth(borders, borders@data$point.1, point.1.shades, border = "gray10")
# dev.off()
