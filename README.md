# PeGS2

Welcome to the first release of PeGS2, the modular, community updated version of the PhotoElastic Grain Solver (originally developed by Johnathan Kollmer). This package takes images from photoelastic granular material and converts it into vector contact and forces. 


___
Basic functionality


This package contains the following modules
* `PeGSModular`: this is the main body of the script, it calls the modules in order, each module can be commented out if you do not need to run it
* `particleDetect`: finds the particles using a circle finding algorithm
* `particleTrack`: tracks particles in subsequent frames and assigns them a persistent id (optional, can comment out for non-sequential images)
* `contactDetect` : detects contacts between potential neighbours using gradients in image intensity
* `diskSolve` : 'solves' the forces by matching the experimental image to a theoretical fringe pattern based
* `adjacencyMatrix`: creates an adjacency list of the particles in contact with each other (optional, but a handy way of seeing the packing)

Each module has a verbose description of the function and explanation of the required input parameters at the top of each function.

Also included are 4 required functions to aid in diskSolve, they should mostly run under the hood without user input. If they do not, please contact us! These functions are: forceBalance.m, fringe\_pattern\_original.m, solver\_original.m, and stress\_engine\_original.m

There are two sample files, the first are 4 sample images located in `testdata/images`. The second is `RunscriptSample.m`, which has a sample of how you can change the input parameters depending on your experiment. There are a lot of input parameters, so we recommend working module by module when you first get started. To get started, you should be able to run either `RunscriptSample.m` or `PeGSModular.m` right out of the box on these images. If you don't specify your inputs, the package will run using defaults that are set at the bottom of each module.

___
Module Structure

Each module takes the basic structure of:

1. File handling (setting up output location, reading data, setting parameters if not set)
2. Main function (usually in a for-loop over each image)
3. Documentation and data saving
4. Other functions
5. Parameter default function (sets default parameters)


___
File structure

To run PeGSModular, the only requirement is that you have photoelastic images, as the user, you'll have to specify their location. We *highly* recommend using the following file formatting:

```
topDir
+--imageDir
	++- image1.jpg
```



