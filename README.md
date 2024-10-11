# PeGS2

Welcome to the first release of PeGS2, the modular, community updated version of the PhotoElastic Grain Solver (originally developed by Johnathan Kollmer). This package takes images from photoelastic granular material and converts it into vector contact and forces. 

We recommend using this version if you are doing photoelasticimetry because we've fixed some bugs, optimized computation, and made it easily adaptible to different user needs with plug-and-play modules.
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
└───imageDir
	└───image1.jpg
```

If you run all of the modules, at the end you will have the following outputs and datastructure

```
topDir
└───imageDir
|	└───image1.jpg
└───particleDir
|	└───image1_centers.txt
|	└───Centers_image1.jpg
|	└───trajectories.fig
|	└───particleDetect_params.txt
└───contactDir
|	└───image1_contacts.mat
|	└───contactDetect_params.txt
|	└───Contacts_image1.jpg
└───solvedDir
|	└───image1_solved.mat
|	└───image1_Synth.jpg
|	└───diskSolve_params.txt
└───adjacencyDir
|	└───image1_Adjacency.txt
|	└───adjacencyMatrix_params.txt
└───Adjacency_list.txt
└───particle_positions.txt
```
___

This README was written by Carmen Lee (cllee3@ncsu.edu) on 01/08/2024. Modules were developed and adapted by Kerstin Nordstrom, Lori McCabe, Abrar Nassar, Ben MacMillan, Ted Brzinski, Karen Daniels, and Carmen Lee.

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
