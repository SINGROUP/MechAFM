MechAFM
=======

Mechanical AFM (implementation based on Hapala et al., Phys. Rev. B, 90:085421, 2014)

Installation
------------
Clone the repository to the desired folder and open it.

```
git clone https://github.com/SINGROUP/MechAFM
cd MechAFM
```

Run 

```
make serial
```
to create the openMP version of the executable or 

```
make mpi
```
to create the hybrid openMP + openMPI version of the executable.

It's recommended to use the openMP version in a single machine enviroments and only use MPI to scale to multiple machines.

To run the openMP version type

```
bin/mechafm-serial INPUT-FILE [OUTPUT-FOLDER]
```
where INPUT-FILE is the input file you want to use (for example: examples/input.scan) and the optional OUTPUT-FOLDER defines the folder where the simulation output will be written (current directory by default).


Input files
-----------

The simulation requires three input files to run: the input file defines the input options of the simulation, the parameter file defines the atom and bond properties of the system and the xyz file defines the atomic structure of the system to be scanned.

###Input file

####Mandatory options
    xyzfile: The xyz file to be used.
  
    paramfile: The parameter file to be used.
  
    tipatom: The type of the tip atom as it is found on the parameter file.
  
    dummyatom: The type of the dummy atom as it is found on the parameter file.
  
    minterm: Term used to check whether the system has converged to a minimum.
             (options: e (energy), f (force) or ef (energy and force))

Either *planeatom* or *zplane* must be defined.

    planeatom: The type of a atom that defines the z-axis zero.
  
    zplane: The value for the z-axis zero.
  
####Other options
  
    units: The units used by the simulation. Affects some constants. (default: kcal/mol)
  
    coulomb: Defines whether coulomb interactions are on or off. (default: off)
  
    dx: Defines the x distance between the simulation points. (default: 0.1)
    
    dy: Defines the y distance between the simulation points. (default: 0.1)
    
    dz: Defines the z distance between the simulation points. (default: 0.1)
    
    zlow: Defines the lowest z point for the dummy atom. (default: 6.0)
    
    zhigh: Defines the starting z for the dummy atom. (default: 10.0)
    
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
                     Steepest Descent minimisation. (Options: euler, rk4 (Runge-Kutta4))
                     (default: rk4)
                     
-
###Parameter file

####General parameters

#####Box
```
box x y z
```

#####Atom
```
atom name epsilon sigma mass charge
```

#####Pair overwrite
```
pair_ovwrt atom1 atom2 morse de a re
pair_ovwrt atom1 atom2 lj epsilon sigma
```
#####Tip harmonic
```
harm name k distance
```

####Flexible parameters


#####Topological bonds
```
topobond atom1 atom2 length
```

#####Harmonic bond
```
bond k
```

#####Harmonic angle
```
angle k
```

#####Harmonic dihedral
```
dihedral k
```

#####Substrate support
```
substrate k
```
-
###XYZ-file

```
atom x y z charge fixed
```
