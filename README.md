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

## How does it work?

**HyATraX** generates approximate transition state geometries that are further optimized with trust-region image minimization. It creates a mol object of the given molecule, detects C-H, N-H, and O-H bonds, and subdivides the C-H bonds in aromatic and aldehyde functionalities. Then, the former R-H bond is elongated to a fixed transition state length before a radical (hydrogen atoms or hydroxyl radicals) is added to the hydrogen atom that is about to be abstracted. The transition state bond lengths and angles were determined empirically for the DFT-D3 TPSS/TZVP method and can be changed in `main.py`. Before the refinement process, the user can visualize all R-H vectors along with the added radicals.  Moreover, various functional groups, abstraction radicals, and modifications of bond lengths and angles for the transition states can be easily added to the source code.

## How to use

**HyATrax** requires an input of the Cartesian coordinates (xyz format) for the molecule of interest. The **main
routine** is invoked by
 * `python3 main_H.py [coord.xyz]` for abstractions by H atoms or
 * `python3 main_OH.py [coord.xyz]` for abstractions by hydroxyl radicals.
   
Then, the program generates for each R-H bond (R=C,N,O) a respective xyz-file with the index of the hydrogen
atom with index [index]:
 * coord_hidx_[index].xyz
   
The user can then refine the transition state guess in three steps using the [TURBOMOLE](https://www.turbomole.org/) program package.
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

## Adjustments 
 * If your molecule is charged adjust this parameter in rdDetermineBonds.DetermineBonds(mol,charge=[your charge])
 * Change RAM, number of Cores/Nodes and [TURBOMOLE](https://www.turbomole.org/) paths in the submission scripts appropriately to your local machine or cluster ([user-guide](https://arc-user-guide.readthedocs.io/en/latest/arc-script-basics.html)) 
 * Check the paths for the input files in the submission scripts as well (default user is "/home/USERNAME")
   
## Requirements

**HyATraX** requires the following Python modules/packages:
 * [numpy](https://numpy.org/)
 * [os](https://github.com/python/cpython/blob/3.12/Lib/os.py)
 * [matplotlib](https://matplotlib.org/)
 * [rdkit](https://www.rdkit.org/docs/index.html)
 * sys

## Acknowledgments
The **HyATraX** code uses the [RDKit](https://www.rdkit.org/docs/index.html) package, an open-source cheminformatics software, 
to detect and analyze R-H bonds in PAHs.  
