# Packing
## Table of Contents

- [About](#about)
- [Usage](#usage)
- [Log](#log)

## About <a name = "about"></a>

This program is for packing hydrate-particles in cylindrical shape.

## Usage <a name = "usage"></a>
The basical usage is as follows:
1. [make a cylindrical mesh ](#1)
2. [set the condition in condition.csv](#2)
3. [run packing.py](#3)
4. [output a geo file](#4)

### make a cylindrical mesh <a name = "1"></a>
Since this program packs particles into a cylindrical mesh, a cylindrical mesh must be created in advance.

In this repository, a cylindrical mesh with radius=2.5mm and height=2.7mm simulating Santamarina's paper is already created and located in ./mesh. But, just in case, how to creating mesh is explained below.


In creating a mesh, [gmsh](http://gmsh.info) is used.
Please install gmsh to your machine(not the server machine.  It requires a display.)


After you installed gmsh, please open ./mesh/cylinder.geo by a texteditor.
Radius and height of a cylinder are defined in the top of .geo file like below. Please change them as you want.
```cylinder.geo
//radius [m]
r = 0.0025;

//height [m]
h = 0.0027;
```

Then, please open cylinder.geo by gmsh.
After opening geo file, select mesh->2D, then cylinder is splitted as a mesh.

<img width="747" alt="mesh" src="https://user-images.githubusercontent.com/50572759/94533443-15351c80-027a-11eb-9055-ebd4e635fc58.png">

Next, export the mesh. Select File->Export, and fill the name as "cylinder" and chose "INRIA Medit" in Save As. Click Save button, then "cylinder" will be output.
<img width="953" alt="save" src="https://user-images.githubusercontent.com/50572759/94533573-444b8e00-027a-11eb-91c5-fa916e9183d9.png">

"cylinder" has to be converted so as esysparticle can read.
medit2timesh.py can do this.
Please just run the below command:
```
python medit2timesh.py
```
Then, cylinder.lsm, which is a mesh file for esysparticle, will be output.




### set the condition in condition.csv <a name = "2"></a>


|variable|explanation|
|---|---|
|endTs|the time step stops the simulation|
|incTs_snap|the interval time step for output of snapshot|
|dT[s]|delta t[s]|
|R[m]|the radius of the cylinder[m]|
|H[m]|the height of the cylinder[m]|
|r_h[m]|the radius of hydrates[m]|
|normalK_wall[N/m]|Normal k of wall|
|normalK_sand[N/m]|Normal k of particles|
|shearK_sand[N/m]|Shear k of particles|
|dynamicMu|Dynamic frictional coefficient|
|staticMu|Static frictional coefficient|
|porosity|porosity|
|viscosity_damp|Viscousity|
|initialfactor|the factor multipled to the particle size at the start of the calculation|

### run packing.py <a name = "3"></a>
To run the packing program, please run the below command.
```
mpirun -np #proc esysparticle packing.py
```

In this program, the radius expansion technique is used to prevent particle overlap. First, the number of particles satisfying the target porosity is calculated and randomly placed in a cylinder. The calculation is started with each particle size multiplied by a certain magnification(initialfactor) to reduce it, and the size of the particles is returned to their original size over time. After the particles have settled, the calculation is finished.
The standard output shows "scaling" when the particle is expanding, and "easing" when it is doing nothing and waiting for the particle to settle down.

There are no clear criteria as to when to stop calculating after "easing" is indicated.
In general, you can visualize the snapshot with paraview and if the particles are not moving, you can stop it.



### output a geo file <a name = "4"></a>
To use the particles packed by this program in the tensile program, you need to output a geo file by
```
dump2geo -i snapshot_packing -o "outputname" -rot -t "target number" 1 1
```
The usage of dump2geo is the same to dump2vtk.

However, there is a bug in dump2geo that it lacks the information "Dimension 3D" in the header and outputs only "Dimension", so you need to add "3D" manually after outputting the geo file. (I will fix this bug in the near future.)
<img width="1235" alt="geo" src="https://user-images.githubusercontent.com/50572759/94570577-f484bb00-02a9-11eb-80fc-a9599449dfbb.png">

Note that the center of the cylinder output from this program goes through the z-axis.

<img width="569" alt="particles" src="https://user-images.githubusercontent.com/50572759/94575968-c904cf00-02af-11eb-9591-e9f9c5707cb4.png">


After you created geo file, please use it in your tensile program.

## Log <a name = "log"></a>
2020.09.29
- First commit.
____