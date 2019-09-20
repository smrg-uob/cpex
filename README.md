cpex: Crystal Plasticity Data Extraction
===============================================

[![Build Status](https://travis-ci.org/casimp/cpex.svg?branch=master)](https://travis-ci.org/casimp/cpex) 
[![Build status](https://ci.appveyor.com/api/projects/status/9cc2aej45li1pm97?svg=true)](https://ci.appveyor.com/project/casimp/cpex/branch/master)

What is cpex?
-------------

cpex has been developed at the University of Bristol (SMRG) to optimise (a) the extraction and (b) the analysis of crystal plasticity finite element (CPFE) data created using Abaqus.

## Data Extraction From .odb File

The data extraction routine is intended for Abaqus odb files and is run on the command line, with all pertinent rotations, stresses and strains being extracted in a single pass. This process has been optimised by minimising the number of calls to the odb file. 
Scrape data is stored into a .npz file, which is a zipped archive of files named after the variables they contain. This allows for the data to be easily reloaded for the subsequent analysis steps.

## CPFE Data Analysis

The analysis modules focuses on interrogating/resolving the lattice plane specific strain response for the CPFE modelled grains. 
During the analysis process, the symmetry of the crystal is accounted for, with the full family or set of crystallographically equvalent planes - {hkl} - being assessed. 
For every grain the strain normal to each of these crystallographic planes calculated/resolved from the grain's orientation and fully described 3D strain state. 
This information is then stored alongside the associated angle that this normal makes with the loading axis. 

The resolved {hkl} specific CPFE data can then be 'caked' into azimuthal slices. 
This is consistent with the approach and results acquired during monochromatic X-ray diffraction and typically involves studying the mean strain (for a given {hkl}) in grains that have reflections or plane normals that are parallel or perpendicular to the loading axis. 
Further use is made of the planes with normals that are not parallel or perpendicular to the loading axis, with an in-plane strain tensor being fit to the azimuthally varying strains. 
This is consistent with the approach used for the pyXe strain analysis software (http://github.com/casimp/pyxe) and is covered in more detail there.

Example Usage
-------------


Requirements
------------

cpex is built on Python's scientific stack (numpy, scipy, matplotlib). Testing and development were carried out using the Anaconda (v 2019.03) package manager, which built with the following versions:

-	Python: version 3.7.3
-	numpy: version 1.16.2
-	scipy: version 1.2.1
-	matplotlib: version 3.0.3

Backward compatability to Python 3.5 is likely but not guaranteed. Scraping functionality requires Abaqus.

Installation
------------

Install from the distribution using the setup.py script. The source is stored in the GitHub repo, which can be browsed at:

https://github.com/casimp/cpex

Simply download and unpack, then navigate to the download directory and run the following from the command-line:

```
pip install .
```