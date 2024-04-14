# List of Packages
library(lidR)
library(sf)
library(terra)

#### Section 1 Preprocessing ####

# Define file paths
las_folder <- "LAS_FOLDER_PATH"
aoi_shapefile <- "AOI_SHAPEFILE_PATH"
output_folder <- "OUTPUT_FOLDER_PATH"

# Read LAS catalog
ctg_aoi <- readLAScatalog(las_folder)

# Clip the area of interest
aoi_clip <- st_read(aoi_shapefile)
ctg_clip_1 <- catalog_intersect(ctg_aoi, aoi_1_clip)
plot(ctg_clip_1)

# Set output directory for filtered LAS layer
opt_output_files(ctg_clip_1) <- paste(file.path(output_folder, "filtered_{ID}"), sep = "")
ctg_filtered_1 <- filter_duplicates(ctg_clip_1)

# Create the DEM
dem_ctg_aoi <- rasterize_terrain(ctg_filtered_1, res = 0.5, algorithm = knnidw())
col_1 <- height.colors(50)
plot(dem_ctg_aoi, col = col_1)

# Normalizing height
opt_output_files(ctg_filtered_1) <- paste(file.path(output_folder, "filtered_normalized", "norm1_{ID}"), sep = "")
norm_aoi <- normalize_height(ctg_filtered_1, dem_ctg_aoi)


#### Section 2 Tree Segmentation ####
# Loading Normalised File
norm_aoi_list <- list.files(file.path(output_folder, "filtered_normalized"), full.names = TRUE) 

# Initialize data frame for merged crown metrics
crowns_merge <- data.frame()

for (i in 1:length(norm_aoi_list)) { # Loop through each normalized file
  tile_path <- norm_aoi_list[i] # Get the file path
  tile <- readLAS(tile_path) # Read the LAS file
  tile_z <- tile[tile@data$Z >= 0 & tile@data$Z < 75, ] # Filter points based on height
  tile_veg <- tile_z[tile_z@data$Classification != 6, ] # Filter out non-vegetation points
  chm <- rasterize_canopy(tile_veg, 1, p2r(0.2, na.fill = tin())) # Generate Canopy Height Model (CHM)
  # Export CHM to file (optional)
  output_chm <- file.path(output_folder, "chm", paste0("chm_", i, ".tif")) # Define output path for CHM
  writeRaster(chm, output_chm, overwrite = TRUE) # Write CHM to file
  
  # Detect treesin CHM
  detected_trees <- locate_trees(chm, algorithm = lmf(3))

  # Segment trees
  segmented_trees <- segment_trees(tile_veg, algorithm = dalponte2016(chm, detected_trees))
  
  # Delineating crowns
  crowns <- crown_metrics(las = segmented_trees, .stdtreemetrics, geom = "concave", concaveman = c(3, 0))
  crowns_merge <- rbind(crowns_merge, crowns)
  print(paste("done", i))
}

# Write merged crown metrics to shapefile
st_write(crowns_merge, file.path(output_folder, "crown_metrics", "crown_merge.shp"), append = FALSE)

#### Section 3 Tree Volume Index Calculation ####
# Read the crown metrics shapefile
fuel <- st_read(output_folder, "crown_merge.shp")
# Calculate fuel volume
fuel_volume <- fuel %>%
  mutate(volume = fuel$Z * fuel$cnvhll)

# Write fuel volume to a new shapefile in the output folder
st_write(fuel_volume, file.path(output_folder, "crown_TVI.shp"))

# The polygon file is now ready to be viewed in GIS software

#### Section 4 Surface Fuel Mapping ####
# Define function for pixel metrics where Z <= 2.5m and Z > 1m
surface_metrics_2 <- function(z) {
  count_total <- length(z)
  count_under_2m <- length(z[z <= 2.5 & z > 1])
  return(
    stratum2 = count_under_2m / count_total)
}

# Define function for pixel metrics where Z equals 0
surface_metrics_0 <- function(z) {
  count_total <- length(z)
  count_0 <- length(z[z == 0])
  return(
    stratum0 = count_0 / count_total)
}

# Initialize lists for storing raster data
rasters2 <- list()
rasters0 <- list()

# Iterate through each normalized file
for (i in 1:length(norm_aoi_list)) {
  tile_path <- norm_aoi_list[i]
  tile <- readLAS(tile_path)
  tile_z <- tile[tile@data$Z >= 0 & tile@data$Z < 75, ]
  tile_veg <- tile_z[tile_z@data$Classification != 6, ]
  
  # Calculate surface metrics for Z <= 2.5m and Z > 1m
  metrics_z_2 <- pixel_metrics(tile_veg, surface_metrics_2(Z), res = 2)
  
  # Calculate surface metrics for Z equals 0
  metrics_z_0 <- pixel_metrics(tile_veg, surface_metrics_0(Z), res = 2)
  
  # Append raster data to lists
  rasters2 <- c(rasters2, metrics_z_2)
  rasters0 <- c(rasters0, metrics_z_0)
  
  print(paste('done', i))
}

# Convert lists to raster objects
rast.list2 <- lapply(1:length(rasters2),
                     function(x) {
                       rast(rasters2[x])
                     })
rast.list0 <- lapply(1:length(rasters0),
                     function(x) {
                       rast(rasters0[x])
                     })

# Set function for combining rasters
rast.list2$fun <- mean
rast.list0$fun <- mean

# Mosaic rasters
rast.mosaic2 <- do.call(mosaic, rast.list2)
rast.mosaic0 <- do.call(mosaic, rast.list0)

# Define output file paths
output_surface_2 <- paste0(output_folder, 'surface_metrics_2.TIFF')
output_surface_0 <- paste0(output_folder, 'surface_metrics_0.TIFF')

# The raster files are now ready to be viewed in GIS software
