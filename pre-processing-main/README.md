# pre-processing

## Files

- **PP_light_correction.mlx**: This script is used to obtain the parameters for correcting the light intensity in the images. It should be run first to get the necessary coefficients for light correction.
- **PeGSDiskPrep_static.m**: This is the main script used for processing data. It handles cropping, illumination adjustment, and visualizes force chains among detected particles.
- **IMG_0078.JPG**: Testing image.
- **PeGSDiskFindH.m**: This function detects particles in the image using the Hough Transform method. It is a support file from version 1.0 of the project.
- **PeGSNeighbourFind.m**: This function identifies neighboring particles and their interactions. It is also a support file from version 1.0 of the project.

## Usage

**1.Run the Light Correction Script**: Execute `PP_light_correction.mlx` to obtain the light correction coefficients. Make sure to save the coefficients for later use.

2.**Process the Image**: Run `PeGSDiskPrep_static.m` to process the image. This script will perform cropping, adjust the illumination based on the coefficients obtained, and visualize the force chains among the detected particles.
