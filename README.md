# HyATraX
## Description 
The Hydrogen Abstraction Transition State Explorer (**HyATraX**) is a workflow to find and optimize transition states of initial
hydrogen abstraction reactions by hydrogen atoms or hydroxyl radicals from polyaromatic hydrocarbons (PAHs) and is also capable of 
handling PAHs with heteroatomic functional groups. Presently, **HyATraX** can handle H-abstractions of 
 * C-**H** bonds (aromatic)
 * C-**H** bonds (aldehyde)
 * O-**H** bonds 
 * N-**H** bonds
   
that can be present in pyrrole, hydroxyl or aldehyde groups in PAHs. It was developed for studying the combustion of PAHs in radical-rich 
atmospheres. The initial steps are usually hydrogen abstractions from the edge of the PAHs by hydrogen atoms or hydroxyl radicals that create 
aryl or oxo radicals. These steps usually require considerable amounts of computation time throughout reactive forcefield molecular dynamics 
simulations. Applying **HyATraX** can skip these initial abstractions by generating reactive intermediates that can be used as input 
structures for MD simulations that already went through common reaction channels. 
## How to use

**HyATrax** requires an input of the Cartesian coordinates (xyz format) for the molecule of interest. The **main
routine** is invoked by
 * `python3 main_H.py [coord.xyz]` for abstractions by H atoms or
 * `python3 main_OH.py [coord.xyz]` for abstractions by hydroxyl radicals.
   
Then, the program generates for each R-H bond (R=C,N,O) a respective xyz-file with the index of the hydrogen
atom with index [index]:
 * coord_hidx_[index].xyz
   
The user can then refine the transition state guess in three steps using the TURBOMOLE program package.
 * default: DFT-D3 TPSS/TZVP
   
If the user decides to change the basis set and functional, the reference values for the bond distances must be
adjusted accordingly.
 * calulate Hessian using **ridft** and **aoforce** (submit_level1.sh)
 * transition state optimization using by **jobex** (submit_level2.sh)
 * verification by calcualting Hessian of optimized transition state using **ridft** and **aoforce** (submit_level3.sh)
   
Use transition state refinement scripts (bash) are invoked by:
 * `./submit_level[1,2,3].sh`

Remember to make the submission scripts executable with:
 * `chmod u+x submit_level[1,2,3].sh`

In the first refinement step **HyATraX** generates for each detected R-H bond a directory
named `hidx_[index].xyz` in which the quantum chemical refinement is performed. 

For each abstracted hydrogen atom the reaction product is also saves as `idx_[index].xyz` in `Products/`. 
The reaction products can be optimized within the `Products/` directory in two levels:
 * geometry optimization using **jobex** (products_1.sh)
 * single point energy using **ridft** and **aoforce** (products_2.sh)
   
The product refinement scripts (bash) can be invoked by:
 * `./products_[1,2].sh`
   
Remember to make the submission scripts executable with:
 * `chmod u+x products_[1,2].sh`
   
After successful refinement of the transition states and reaction products, the directories can be  cleaned using the `clean.sh` script.
Clean the directory by invoking:
 * `./clean.sh`
   
Lastly, the single point energies of the transition states and reaction products can be extracted and used for the calculation of
activation and reaction energies of H-abstraction reactions.

## Requirements

**HyATraX** requires the following Python modules/packages:
 * numpy
 * sys
 * os
 * matplotlib
 * rdkit