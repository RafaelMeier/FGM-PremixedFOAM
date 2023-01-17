# FGM-PremixedFOAM

The FGM-PremixedFOAM is an OpenFOAM-based application that stands for Flamelet-Generated Manifold for direct numerical simulations of premixed flames.

The Flamelet-Generated Manifold is a validated and widely used computational chemical reduction technique based on the assumption that most chemical time scales are very small compared to the flow scales. In a time-scale analysis, the chemical scales are assumed to be at steady-state compared to the convection and diffusion transport processes.

The ```fgmPremixedFoam``` app is built from the ```rhoReactingBuoyantFoam``` where a new library named ```combustionFGMModel``` is called by the OpenFOAM. The ```rhoReactingBuoyantFoam``` is a native OpenFOAM solver that comes as an optional platform in ```reactingFoam```. The ```rhoReactingBuoyantFoam``` is a solver for chemical reactions using a density-based thermodynamics package with enhanced buoyancy treatment.

### Version

OpenFOAM-v2112

### Solver compilation

The solver compilation can be performed by running the provided bash script in the fgmPremixedFoam/
folder

```
bash Allwmake
```

or can be done manually typing in terminal window

```
wmakeLnInclude -u combustionFGMModel
wmake combustionFGMModel
wmake fgmPremixedFoam
```

### Test case

It is provided a 2D Bunsen flame test case for a two-dimensional slab domain subject to a prescribed  Poiseuille inlet velocity for stochiometric methane/air mixture at atmospheric conditions.

To run the case just type in the terminal window the following command

```
blockMesh
fgmPremixedFoam >& log&
```

The steady-state solution of the proposed case should result in the contours according to Figures below.


<img src="/images/sourcePV.png" width="300"/> <img src="/images/Temperature.png" width="300"/> 
