import os
import shutil
import numpy as np
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *
from packingRunnable import packingRunnable
import packing_lib as plib

wallNames = []
wallNormals = []
interactionNames = []
vPlates = []
contactAreas = []

out_dat_dir = "out_data_packing"

if os.path.exists("./"+out_dat_dir):
    print(out_dat_dir, "exist. Cleared contents.")
    shutil.rmtree(out_dat_dir)
    os.mkdir(out_dat_dir)
else:
    print(out_dat_dir, "doesn't exist. made")
    os.mkdir(out_dat_dir)

con_d = np.genfromtxt(
    "conditions.csv", delimiter=" ", dtype="U", autostrip=True)
for c in con_d:
    c = c.split(",")
    print(c)
    if c[0] == "endTs":
        endTs = int(c[1])
    if c[0] == "incTs_snap":
        incT_snap = int(c[1])
    if c[0] == "dT[s]":
        dT = float(c[1])
    if c[0] == "H[m]":
        H = float(c[1])
    if c[0] == "R[m]":
        R = float(c[1])
    if c[0] == "r_h[m]":
        r_h = float(c[1])
    if c[0] == "density[kg/m^3]":
        density = float(c[1])
    if c[0] == "normalK_wall[N/m]":
        normalK_wall = float(c[1])
    if c[0] == "normalK_sand[N/m]":
        normalK_sand = float(c[1])
    if c[0] == "shearK_sand[N/m]":
        shearK_sand = float(c[1])
    if c[0] == "dynamicMu":
        dynamicMu = float(c[1])
    if c[0] == "staticMu":
        staticMu = float(c[1])
    if c[0] == "porosity":
        porosity = float(c[1])
    if c[0] == "viscosity_damp":
        viscosity_damp = float(c[1])
    if c[0] == "initialfactor":
        initialfactor = float(c[1])

endT_pack = endTs
d_h = r_h*1e+3*2

pobj = plib.Pack()
pobj.setConstDiam(d=d_h)
pobj.setInitFac(initialfactor)
pobj.setCylinder(R=R*1e+3, H=H*1e+3)
pobj.setTargetPorosity(poro=porosity)
pobj.makeParticle_constD()
pobj.makeGeofile()

ep = 0
speed = 0.0
# top
wn_top = "top_wall"
in_top = "tw_repell"
vp_top = Vec3(0.0, 0.0, -speed)
nm_top = Vec3(0.0, 0.0, -1.0)
ca_top = np.pi*R**2
pos_top = Vec3(0.0000, 0.0000, H+ep)
# pos_top = Vec3(0.0000, 0.0000, H)
wallNames.append(wn_top)
wallNormals.append(nm_top)
interactionNames.append(in_top)
vPlates.append(vp_top)
contactAreas.append(ca_top)

# bottom
wn_bottom = "bottom_wall"
in_bottom = "btw_repell"
vp_bottom = Vec3(0.0, 0.0, speed)
nm_bottom = Vec3(0.0, 0.0, 1.0)
ca_bottom = np.pi*R**2
pos_bottom = Vec3(0.0000, 0.0000, 0.0000-ep)
# pos_bottom = Vec3(0.0000, 0.0000, 0.0000)
wallNames.append(wn_bottom)
wallNormals.append(nm_bottom)
interactionNames.append(in_bottom)
vPlates.append(vp_bottom)
contactAreas.append(ca_bottom)


# instantiate a simulation object and
# initialise the neighbour search algorithm:
sim = LsmMpi(numWorkerProcesses=4, mpiDimList=[1, 1, 4])
sim.initNeighbourSearch(
    particleType="RotSphere",
    gridSpacing=2.5*r_h,
    verletDist=0.2*r_h
)

# sim.setVerbosity(True)


# specify the number of timesteps and the timestep increment:
sim.setNumTimeSteps(endT_pack)
sim.setTimeStepSize(dT)

# ep = 0.0
ep = H*0.01
# specify the spatial domain for the simulation:
domain = BoundingBox(Vec3(-R-ep, -R-ep, 0-ep), Vec3(R+ep, R+ep, H+ep))
# domain = BoundingBox(Vec3(-H, -H, -H), Vec3(H, H, H))
sim.setSpatialDomain(domain)

sim.readGeometry("ini.geo")
p_list = sim.getParticleList()
Np = len(p_list)
for i, p in enumerate(p_list):
    pid = p.getId()
    sim.setParticleTag(pid, 1)
    p.setTag(1)
sim.setParticleDensity(
    tag=1,
    mask=-1,
    Density=2500e+9
)


# add a back wall to the model:
sim.createWall(
    name=wn_bottom,
    posn=pos_bottom,
    normal=nm_bottom
)

# add a front wall to the model:
sim.createWall(
    name=wn_top,
    posn=pos_top,
    normal=nm_top
)


# specify that particles undergo elastic repulsion
# from the front wall:
sim.createInteractionGroup(
    NRotElasticWallPrms(
        name=in_bottom,
        wallName=wn_bottom,
        normalK=normalK_wall*1
    )
)


# specify that particles undergo elastic repulsion
# from the front wall:
sim.createInteractionGroup(
    NRotElasticWallPrms(
        name=in_top,
        wallName=wn_top,
        normalK=normalK_wall*1
    )
)


sim.readMesh(
    fileName="mesh/cylinder.lsm",
    meshName="cylinder"
)


sim.createInteractionGroup(
    RotFrictionPrms(
        name="friction",
        normalK=normalK_sand,
        dynamicMu=dynamicMu,
        staticMu=staticMu,
        shearK=shearK_sand
    )
)

sim.createInteractionGroup(
    LinDampingPrms(
        name="damping",
        viscosity=viscosity_damp,
        maxIterations=10
    )
)


sim.createInteractionGroup(
    NRotElasticTriMeshPrms(
        name="cylinder_repell",
        meshName="cylinder",
        normalK=normalK_wall*1
    )
)


packing = packingRunnable(
    LsmMpi=sim,
    porosity=porosity,
    R=R,
    H=H
)
sim.addPreTimeStepRunnable(packing)


# add a CheckPointer to store simulation data:
sim.createCheckPointer(
    RestartCheckPointPrms(
        fileNamePrefix="snapshot_packing",
        beginTimeStep=0,
        endTimeStep=endTs,
        timeStepIncr=incT_snap
    )
)


# execute the simulation:
sim.run()
