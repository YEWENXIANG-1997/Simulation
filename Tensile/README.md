# tensile test for methane hydrates

## About <a name = "about"></a>

Tensile test for methane hydrates using esysparticle.

## Usage <a name = "usage"></a>

This project consists of the followings:
 - tensile.py : a main program for the tensile test
 - tensileRunnable.py : a runnable class for moving particles
 - conditions.csv : a input file
 - plot.py: a program for plots.

To run this program, execute the following command.
```
mpirun -np 2 esysparticle tensile.py
```
After a computation, results will be output in "out_data" directory.
Its contents are :
 - nbonds.dat, the number of bonds.
 - wall_Force.dat, the force acting on the walls.
 - wall_Position.dat, the position of the wall the walls.

The bond will be broken when the tensile stress exceeds 0.2MPa, 
but the compuation will continue 
until the strain reaches the value set in conditions.csv(endstrain).
Because results after bonds broken are just meaningless, 
please hault the computation by "Ctrl + c".
To check if the bond is broken or not, check "nbonds.dat".
If its value is 0, the bond is broken.



### conditions.csv
The contens in conditions.csv are explained below.

|variable|explanation|
|---|---|
|dT[s]|delta t[s]|
|incT[s]|the interval for output of data[s]|
|incT_snap[s]|the interval for output of snapshot|
|speed[m/s]|the speed of the particle[m/s]|
|H[m]|the height of the cylinder|
|R[m]|the radius of the cylinder|
|r_p[m]|the radius of particles|
|endstrain|max of the strain.|
|density[kg/m^3]|density of particles|
|maxTensile[Pa]|max tensile stress that the bond will break|
|Sh|the saturation of methane hydrates. Set as 1.0.|
|bondModulus|Young's modulus of methane hydrates. Set as 300e+6.|
|bondPoissonRatio|Poisson ratio of methane hydrates. Set as 0.31403.|
|bondCohesion|Cohesion of methane hydrates. Set as 2.81e+6.|
|bondTanAngle|Angle of MohrCoulomb criterion. Set as 4.0.|
|beta1|a param for BrittleBeamIG. Set as 1.0.|
|beta2|a param for BrittleBeamIG. Set as 1.0.|
|k0|a cohesion coefficient for MHbond. Set as 6.25.|
|k1|a Young's modules coefficient for MHbond. Set as 12.5.|
|normalK_wall[N/m]|Normal kn of wall. Set as 10.0e8.|
|normalK_sand[N/m]|Normal kn of particles. Set as 3.5e9.|
|kn/ks|a ratio of kn to ks of particles. Set as 10.0.|
|dynamicMu|Dynamic frictional coefficient. Set as 0.001.|
|staticMu|Static frictional coefficient. Set as 0.4.|
|viscosity_damp_nr|Viscousity. Set as 1.00|


### plot
plot.py is for plotting some figures.
Just exexute the following command.
```
python plot.py
```
Then, some figures are output in "/fig" directory. 
This program uses matplotlib. 
You can edit the code as you 
want.
