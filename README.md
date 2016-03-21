MechAFM
=======

Mechanical AFM (implementation based on the Probe Particle Model [1] and inspired by [2]).

Installation
============
Clone the repository to the desired folder and open it.

```
git clone https://github.com/SINGROUP/MechAFM
cd MechAFM
```

Run

```
make omp
```
to create the openMP version of the executable or

```
make mpi
```
to create the hybrid openMP + openMPI version of the executable.

It's recommended to use the openMP version in a single machine enviroments and only use MPI to scale to multiple machines.

To run the openMP version type

```
bin/mechafm-omp INPUT-FILE [OUTPUT-FOLDER]
```
where INPUT-FILE is the input file you want to use (for example: examples/input.scan) and the optional OUTPUT-FOLDER defines the folder where the simulation output will be written (current directory by default).

##Windows

The openMP version of MechAFM has been tested to be working on Windows with [MinGW](http://www.mingw.org/). Instructions for installing MinGW can be found [here](http://www.mingw.org/wiki/Getting_Started). Otherwise just follow the instructions above to compile.


Input files
===========

The simulation requires three input files to run: the input file defines the input options of the simulation, the parameter file defines the atom and bond properties of the system and the xyz file defines the atomic structure of the system to be scanned.

##Input file

###Mandatory options
    xyzfile: The xyz file to be used.
  
    paramfile: The parameter file to be used.
  
    tipatom: The type of the tip atom as it is found on the parameter file.
  
    dummyatom: The type of the dummy atom as it is found on the parameter file.
  
    minterm: Term used to check whether the system has converged to a minimum.
             (options: e (energy), f (force) or ef (energy and force))
  
###Other options
  
    units: The units used by the simulation. Affects some constants. (default: kcal/mol)
  
    coulomb: Defines whether coulomb interactions are on or off. (default: off)
    
    use_external_potential: Defines whether the tip atom interacts with an external electrostatic
                            potential. Works only for rigid systems with coulomb off. (default: off)
  
    e_potential_file: The file from which the electrostatic (Hartree) potential is read if 
                      use_external_potential is on. Only cube files supported at the moment.
  
    area: Defines the size of the simulation area in x and y. (default: 10.0 10.0)

    center: Defines the position of the molecules center of mass in x and y. 
            (default: area / 2)
    
    zlow: Defines the lowest z point for the dummy atom. (default: 6.0)
    
    zhigh: Defines the starting z for the dummy atom. (default: 10.0)

    dx: Defines the x distance between the simulation points. (default: 0.1)
    
    dy: Defines the y distance between the simulation points. (default: 0.1)
    
    dz: Defines the z distance between the simulation points. (default: 0.1)
    
    etol: Defines the energy value used to check for convergence. (default: 0.01)
    
    ftol: Defines the force value used to check for convergence. (default: 0.01)
    
    dt: Defines the time step of the simulation. (default: 0.001)
    
    maxsteps: Defines the maximum number of minimisation steps used for single tip position.
              (default 5000)
    
    bufsize: Defines the size of the output buffer. Ie. how many xyz points can be stored before they're 
             written to a file. (default: 1000)
    
    gzip: Defines whether the output files are gzipped or not. (default: on)
    
    flexible: Defines whether the whole system is allowed to move or just the tip. (default: off)
    
    rigidgrid: Defines whether the tip forces are precomputed on a grid or not. 
               Can only be used on rigid systems! (default: off)
    
    minimiser_type: Defines the type of a minimiser to be used. (Options: SD (Steepest Descent), FIRE)
                    (default: FIRE)
    
    integrator_type: Defines the type of a integrator to be used. Doesn't do anything with 
                     Steepest Descent minimisation. (Options: euler, midpoint, rk4 (Runge-Kutta 4))
                     (default: midpoint)
                     

##Parameter file

###General parameters

General parameters are used by all simulations.

####Atom
```
atom name epsilon sigma mass charge
```
Atom defines the parameters for a atom with 'name'. Each atom type in the xyz file must have parameters set.

####Pair overwrite (optional)
```
pair_ovwrt atom1 atom2 morse de a re
pair_ovwrt atom1 atom2 lj epsilon sigma
```
Pair overwrite gives specific parameters for the interactions between atoms of type 'atom1' and 'atom2'. Morse option will use morse potential between the atoms and lj option will use Lennard-Jones potential. Without overwrites Lennard-Jones potentials will be used and the parameters will be read from the corresponding atom parameters.

####Tip harmonic
```
harm k distance
```
Tip harmonic defines the harmonic constant for the tip harmonic constraint in the xy-plane. Distance defines the rest xy-distance of the tip atom.

###Flexible parameters

Flexible parameters are only used by simulations with flexibility enabled.

####Topological bonds
```
topobond atom1 atom2 length
```
Topological bonds define the types of bonds present in the system. Atom1 and atom2 give the types of atoms the bonds connect and length gives the expected length of the bond.

####Harmonic bond
```
bond k
```
Bond defines the harmonic constant for harmonic bond interactions.

####Harmonic angle
```
angle k
```
Angle defines the harmonic constant for harmonic angle interactions.


####Harmonic dihedral
```
dihedral k
```
Dihedral defines the harmonic constant for harmonic dihedral interactions.


####Substrate support
```
substrate epsilon sigma lambda r_cut (lateral_k)
```
Substrate defines the force constants for the substrate support. Epsilon, sigma, lambda and r_cut define parameters for the  10-4 wall potential [3] and the lateral_k defines the harmonic constant used for the optional xy-plane harmonic potential (see XYZ-file options).


##XYZ-file
The simulation reads xyz files in the following format:
```
atom x y z charge fixed
```
where atom gives the type of the atom, x y z the position of the atom. The molecule will be placed in the simulation such that the xy center of the molecule will be at the position defined by the center input variable and the z value is defined by placing the lowest non-hydrogen atom at z = 0. The Charge defines the charge of the atom and fixed whether the atom is fixed in place (1), free to move (0) or if it's bound by harmonic potential in xy-plane (2) (only used by flexible simulations).

Output format
=============

The simulation will produce output in the following format:
```
z_index x_index y_index pos.x pos.y pos.z f.x f.y f.z r.x r.y r.z r.len angle energy n_steps
```
where pos is the position of the dummy atom, f is the force on the tip caused by the surface, r is the difference between tip and dummy positions, angle is the angle defined by the r-vector and z-axis, energy is the energy of the tip caused by the surface and n_steps the amount of minimisation steps required for this point.

References
==========

1. P. Hapala, G. Kichin, C. Wagner, F.S. Tauts, R. Temirov, and P. Jelinek, *Mechanism of high-resolution STM/AFM imaging with functionalized tips*, Phys. Rev. B, 90:085421 (2014)
2. S. Hamalainen, N. van der Heijden, J. van der Lit, S. den Harto, P. Liljeroth, and I. Swart, *Intermolecular Contrast in Atomic Force Microscopy Images without Intermolecular Bonds*, Phys. Rev. Lett, 113:186102 (2014).
3. P. Spijker, H.M.M. ten Eikelder, A.J. Markvoort, S.V. Nedea, and P.A.J. Hilbers, *Implicit Particle Wall Boundary Condition in Molecular Dynamics*, J. Mech. Eng. Sci., 222:855 (2008).
