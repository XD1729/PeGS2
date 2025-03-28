# imagePreprocess Module

## Introduction
The `imagePreprocess` module is a new addition to the PeGS2 (PhotoElastic Grain Solver) package. It is designed to improve image quality by correcting lighting issues and enhancing contrast. This preprocessing step significantly improves the accuracy of particle detection and force calculation in subsequent modules.

## Functionality
The module performs the following operations:
1. **Lighting Correction**: Evens out lighting across the image based on the distance from the center
2. **Contrast Enhancement**: Adjusts image contrast for better visibility of features
3. **Auto-detection**: Automatically identifies image type and applies appropriate processing parameters
4. **Obstacle Masking**: Blanks out regions containing obstacles or unwanted artifacts
5. **Region-specific Processing**: Applies different contrast adjustments to specific regions of the image

## Image Types and Auto-detection
The module supports multiple image types through its auto-detection feature. Current supported types include:

### 1. test.JPG Type
Images similar to `test.JPG` typically feature:
- Strong illumination gradients requiring lighting correction
- Need for specific contrast enhancement settings
- Triangular obstacle region at the bottom
- Require region-specific contrast adjustments to handle varied lighting conditions

### 2. Step09 Series Type
Images in the Step09 series typically feature:
- More uniform lighting, requiring minimal to no lighting correction
- Different contrast enhancement settings
- No obstacle regions requiring masking
- Uniform processing across the entire image

### Auto-detection Process
The auto-detection feature analyzes images using:
1. Filename analysis (looking for keywords like "test" or "step")
2. Image characteristics (brightness distribution, RGB channel analysis)

Based on this analysis, appropriate processing parameters are automatically applied.

## Integration with PeGS2
The `imagePreprocess` module is designed to be inserted before the `particleDetect` module in the PeGS workflow:

```
PeGSModular -> imagePreprocess -> particleDetect -> particleTrack -> contactDetect -> diskSolve -> adjacencyMatrix
```

Both `particleDetect` and `contactDetect` modules have been updated to automatically check for and use the enhanced images when available.

## Output Files
For each input image, the module produces:
- `[imagename]_enhanced.jpg`: Image with even lighting and enhanced contrast
- `imagePreprocess_params.txt`: Record of parameters used

All processed images are stored in the `processedImgDir` directory (default: `processed_images`).

## Parameters
The following parameters can be set in the `ipParams` structure:

### Auto-detection Control
- `auto_detect`: Enable/disable automatic image type detection (boolean, default: true)

### Basic Parameters
- `light_correction_coefficients`: Coefficients for polynomial light correction [a, b, c]
- `I_low`, `I_high`: Intensity adjustment limits for contrast enhancement

### Advanced Processing Options
- `obstacle_mask`: Apply masking to obstacle regions (boolean, default: true for test.JPG type)
- `region_specific`: Apply region-specific contrast enhancement (boolean, default: true for test.JPG type)

## Sample Usage
In `RunscriptSample.m`:

```matlab
%% imagePreprocess parameters with auto-detection
ipParams.auto_detect = true; % automatically detect image type and apply appropriate parameters

% You can still override specific parameters if needed
ipParams.I_low = 0.08; % override the auto-detected low intensity limit
ipParams.obstacle_mask = false; % disable obstacle masking even if auto-detected
```

## Manual Configuration
If you prefer to manually set parameters for specific image types:

### For test.JPG Type Images:
```matlab
%% Manual settings for test.JPG images
ipParams.auto_detect = false; % disable auto-detection
ipParams.light_correction_coefficients = [-9.268e-09, 1.192e-06, 0.076]; 
ipParams.I_low = 0.09;
ipParams.I_high = 1;
ipParams.obstacle_mask = true; % apply triangle obstacle masking
ipParams.region_specific = true; % apply different processing to corners and edges
```

### For Step09 Series Images:
```matlab
%% Manual settings for Step09 series
ipParams.auto_detect = false; % disable auto-detection
ipParams.light_correction_coefficients = [0, 0, 0]; % No lighting correction
ipParams.I_low = 0.05;
ipParams.I_high = 0.95;
ipParams.obstacle_mask = false; % no obstacle masking needed
ipParams.region_specific = false; % uniform processing
```

## Deriving Light Correction Coefficients
The `light_correction_coefficients` can be derived using the following approach:
1. Take a photo of a uniformly reflective surface under the same lighting conditions
2. Plot intensity vs. distance from center
3. Fit a polynomial curve to obtain coefficients

The default coefficients were derived for a typical circular photoelastic experiment setup.

## Region-specific Processing
Region-specific processing divides the image into several regions:
- Top-left corner
- Top-right corner
- Left edge
- Right edge
- Main central region

Each region gets separate contrast enhancement, which is particularly useful for images with varying lighting conditions across different areas. This approach closely mirrors the PeGSDiskPrep_static implementation.

---

*This module was extracted from PeGSDiskPrep_static.m and adapted for the modular PeGS2 framework.* 