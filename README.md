# PolyMorph3D

## Overview
A flexible MATLAB-based tool for quantifying microglia morphology and multi-colocalization in three dimensions, PolyMorph3D is a script-based MATLAB program mediated by graphical user interface (GUI) apps; thus, coding experience is not necessary. The user starts with a 3D stack (which must be TIFF), thresholds and identifies valid structures within channels, then colocalizes structures from up to 4 channels with the identified microglia. After processing a single image manually with the Main Script, a parameters file is saved and may be used for batch processing many TIFF images.

## Data Output
The skeletonization step produces the following microglia morphology information: total branching length, number of branches, number of branch points, and number of endpoints. The cell volume and territorial volume of each microglia is also saved. In addition to these single-cell resolution outputs, the script produces image-level data for mean fluorescent intensity for all channels processed, total microglia volume, image volume, maximum microglia volume, average distance between microglia, microglia density, total number of objects in other processed channels, density of objects in other processed channels, and total and average volume of objects in other processed channels. 

## Other Resources
This program is based on 3DMorph, which can be found here: 
https://github.com/ElisaYork/3DMorph
