# LiDAR Data Processing in R

This repository contains R code for processing LiDAR data collected using airborne LiDAR (Light Detection and Ranging) technology. The LiDAR metrics derived could be used at forest structure mapping for wildfire risk assessment. Individual tree segmentation is utilized for crown fuel and ladder fuel mapping, followed by the calculation of relative fuel volume. Pixel metrics are employed for detecting the presence of surface fuel.

## Data Processing and Analysis

Data for this study was collected using airborne LiDAR technology. After collection, the LiDAR data underwent processing stages to prepare it for analysis:

- **Data Normalization**: This initial step involved creating a Digital Elevation Model (DEM) and normalizing the height data to correct for terrain variations.

The methodology employed for assessing wildfire risks through LiDAR data encompasses several advanced techniques:

- **Canopy Height Model (CHM)**: Creation of CHM to assess vegetation height across the study area.
- **Individual Tree Identification and Segmentation**: Algorithms were applied to identify and segment individual trees within the LiDAR data, which are crucial for detailed fuel modeling.
- **Self-Defined Functions and Pixel Metrics**: Custom functions were developed to refine the analysis and utilize pixel-level data for accurate fuel mapping.

## Sections

This repository's code is organized into several sections:

### 1. Preprocessing

This section includes code for preprocessing LiDAR data, including clipping, filtering, creating a digital elevation model (DEM), and normalizing heights.

### 2. Tree Segmentation

Here, we segment trees from the preprocessed LiDAR data and delineate their canopy height and crown area.

### 3. Calculating Fuel Volume

In this section, we calculate fuel volume using functions to find out the relative volume of each individual tree surface for further volume density analysis.

## Required R Packages

To run the code in this repository, you'll need the following R packages:

- `terra`
- `lidR`
- `sf`
- `tidyverse`

Install these packages using `install.packages("package_name")` if you haven't already.

## Usage

1. Clone the repository to your local machine.
2. Open R or RStudio.
3. Set the working directory to the cloned repository.
4. Run the R scripts in sequence, adjusting file paths as necessary.
5. View the results in GIS software.

## Contributors

- Kit Shan Cheung (Wendy Cheung) - Main contributor

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

