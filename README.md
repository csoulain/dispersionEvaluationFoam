# dispersionEvaluationFoam


**************************************************************************
             dispersionEvaluationFoam for OpenFOAM  
**************************************************************************

## General Informations

- dispersionEvaluationFoam is an OpenFOAM package to compute the dispersion
  tensor of a porous sample


- The method relies on the volume averaging approach and the closure problems
  derived by Carbonell and Whitaker (1983)

- The package contains:
        1. The core solver
        2. The associated boundary condition
        3. a tutorial

## History

First development in 2014 by C.S. First release on May 25,2021.


## Installation instructions :

- This toolbox has been tested on OpenFOAM v7 only

- It only needs a standard OpenFOAM installation from www.openfoam.org

- Read the COPYING_OPENFOAM file for information about OpenFOAM and this
  toolbox Copyrights.

- First, source the OpenFOAM configuration file, i.e. (example for ubuntu
  version) :

  	  source /opt/openfoamv7/etc/bashrc

- then in the "porousMedia4Foam" directory, run :

      wmake

  to install the package.


- see the ReleaseNotes.txt file for detailed information about the toolbox.
