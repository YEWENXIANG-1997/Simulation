
import os
import shutil
import numpy as np
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *
from tensileRunnable import tensileClass

def makedirectory(dirname):
    if os.path.exists(dirname):
        print(dirname, "exist. Cleared contents.")
        shutil.rmtree(dirname)
        os.mkdir(dirname)
    else:
        print(dirname, "doesn't exist. made")
        os.mkdir(dirname)


bBond = False

out_dat_dir = "out_data"
disp_dir = "displacement"
force_dir = "nforce"
position_dir = "position"
sigmad_dir = "sigmaD"

dirname = "./"+out_dat_dir
makedirectory(dirname)
# dirname = "./"+out_dat_dir + "/" + disp_dir
# makedirectory(dirname)
dirname = "./"+out_dat_dir + "/" + force_dir
makedirectory(dirname)
# dirname = "./"+out_dat_dir + "/" + position_dir
# makedirectory(dirname)
# dirname = "./"+out_dat_dir + "/" + sigmad_dir
# makedirectory(dirname)


con_d = np.genfromtxt(
    "conditions.csv", delimiter=" ", dtype="U", autostrip=True)
for c in con_d:
    c = c.split(",")
    print(c)
    if c[0] == "dT[s]":
        dT = float(c[1])
    if c[0] == "incT[s]":
        incT_data = float(c[1])
        incTs_data = int(incT_data/dT)
    if c[0] == "incT_snap[s]":
        incT_snap = float(c[1])
        incTs_snap = int(incT_snap/dT)
    if c[0] == "speed[m/s]":
        speed = float(c[1])
    if c[0] == "H[m]":
        H = float(c[1])
        iniD = H
    if c[0] == "R[m]":
        R = float(c[1])
    if c[0] == "endstrain":
        t = iniD*float(c[1]) / (2.*speed)
        endTs = int(t / dT)
        # endTs = 100000
    if c[0] == "r_p[m]":
        r_p = float(c[1])
    if c[0] == "density[kg/m^3]":
        density = float(c[1])
    if c[0] == "maxTensile[Pa]":
        maxTensile = float(c[1])
    if c[0] == "Sh":
        Sh = float(c[1])
    if c[0] == "k0":
        k0 = float(c[1])
    if c[0] == "k1":
        k1 = float(c[1])
    if c[0] == "bondModulus":
        bondModulus = float(c[1])
    if c[0] == "bondPoissonRatio":
        bondPoissonRatio = float(c[1])
    if c[0] == "bondCohesion":
        bondCohesion = float(c[1])
    if c[0] == "bondTanAngle":
        bondTanAngle = np.tan(np.radians(float(c[1])))
    if c[0] == "beta1":
        beta1 = float(c[1])
    if c[0] == "beta2":
        beta2 = float(c[1])
    if c[0] == "normalK_wall[N/m]":
        normalK_wall = float(c[1])
    if c[0] == "normalK_sand[N/m]":
        normalK_sand = float(c[1])
    if c[0] == "kn/ks":
        alpha = float(c[1])
        shearK_sand = normalK_sand / alpha
    if c[0] == "dynamicMu":
        dynamicMu = float(c[1])
    if c[0] == "staticMu":
        staticMu = float(c[1])
    if c[0] == "viscosity_damp_nr":
        viscosity_damp_nr = float(c[1])

ep = 1e-2*H

# instantiate a simulation object and
# initialise the neighbour search algorithm:
sim = LsmMpi(numWorkerProcesses=1, mpiDimList=[1, 1, 1])
sim.initNeighbourSearch(
    particleType="RotSphere",
    gridSpacing=2.5*r_p,
    verletDist=0.2*r_p
)


#sim.setVerbosity(True)

# specify the number of timesteps and the timestep increment:
sim.setNumTimeSteps(endTs)
sim.setTimeStepSize(dT)
#
ep = 0.01*H
# ep = 0.0
# specify the spatial domain for the simulation:
domain = BoundingBox(Vec3(-R-ep, -R-ep, -0-ep), Vec3(R+ep, R+ep, H+ep))
sim.setSpatialDomain(domain)



sim.readGeometry("R3H3.geo")

plist = []

for p in sim.getParticleList():
    x, y, z = p.getCenter().toList()
    pid = p.getId()
    m = p.getMass()
    rad = p.getRad()
    particle = SimpleSphere(id=pid, posn=Vec3(
        x, y, z), radius=rad, mass=m)
    plist.append(particle)


Np = len(plist)
ep = H*0.1
print("Number of Particle:{}".format(Np))

print("setting Tag")
for i, p in enumerate(plist):
    x, y, z = p.getCentre().toList()
    pid = p.getId()
    if z < 0 + 1.00*r_p or z > H - 1.00*r_p:
        p.setTag(2)
        sim.setParticleTag(pid, 2)
    else:
        p.setTag(1)
        sim.setParticleTag(pid, 1)



sim.setParticleDensity(tag=1, mask=-1, Density=density)
sim.setParticleDensity(tag=2, mask=-1, Density=density)

sim.setParticleNonRotational(tag=1)
sim.setParticleNonRotational(tag=2)





sim.createInteractionGroup(
    RotFrictionPrms(
        name="friction",
        normalK=normalK_sand,
        dynamicMu=dynamicMu,
        staticMu=staticMu,
        shearK=shearK_sand,
        scaling=True,
        rigid=False,
        meanR_scaling=True
    )
)


sim.createConnections(
    ConnectionFinder(
        maxDist=r_p*0.1,
        bondTag=1,
        pList=plist
    )
)
sim.createInteractionGroup(
    BrittleBeamPrms(
        name="ppbond",
        youngsModulus=bondModulus*k1*Sh,
        poissonsRatio=bondPoissonRatio,
        cohesion=bondCohesion*k0*Sh,
        tanAngle=bondTanAngle,
        tag=1,
        meanR_scaling=True,
        truncated=maxTensile,
        beta1=1.0,
        beta2=1.0
    )
)

sim.createExclusion(
    interactionName1="ppbond",
    interactionName2="friction"
)

sim.createInteractionGroup(
    LinDampingPrms(
        name="damping",
        viscosity=viscosity_damp_nr,
        maxIterations=100
    )
)

sim.createFluidForce(
    FluidForcePrms(
        name="lbm"
    )
)


# add a back wall to the model:
sim.createWall(
    name="top",
    posn=Vec3(0,0,H),
    normal=Vec3(0,0,-1)
)

# add a back wall to the model:
sim.createWall(
    name="bottom",
    posn=Vec3(0,0,0),
    normal=Vec3(0,0,+1)
)

sim.createInteractionGroup(
    NRotBondedWallPrms(
        name="ptop",
        wallName="top",
        normalK=normalK_wall,
        particleTag=2
    )
)
sim.createInteractionGroup(
    NRotBondedWallPrms(
        name="pbottom",
        wallName="bottom",
        normalK=normalK_wall,
        particleTag=2
    )
)



# # add a wall loader to move the back wall:
tensile = tensileClass(
    LsmMpi=sim,
    speed=speed,
    endTs=endTs,
    iniD=iniD
)
sim.addPreTimeStepRunnable(tensile)


# add a CheckPointer to store simulation data:
sim.createCheckPointer(
    CheckPointPrms(
        fileNamePrefix="snapshot",
        beginTimeStep=0,
        endTimeStep=endTs,
        timeStepIncr=incTs_snap
    )
)

# add a CheckPointer to store simulation data:

sim.createFieldSaver (
    InteractionScalarFieldSaverPrms(
        interactionName="ppbond",
        fieldName="count",
        fileName=out_dat_dir + "/nbonds.dat",
        fileFormat="SUM",
        beginTimeStep=0,
        endTimeStep=endTs,
        timeStepIncr=incTs_data
   )
)

# sim.createFieldSaver (
#     InteractionVectorFieldSaverPrms(
#         interactionName="ppbond",
#         fieldName="normal_force",
#         fileName=out_dat_dir + "/" + force_dir + "/nforce",
#         fileFormat="RAW2",
#         beginTimeStep=0,
#         endTimeStep=endTs,
#         timeStepIncr=incTs_data
#    )
# )



# sim.createFieldSaver(
#     ParticleVectorFieldSaverPrms(
#         fieldName="displacement",
#         fileName=out_dat_dir + "/" + disp_dir + "/displacement",
#         fileFormat="RAW2",
#         beginTimeStep=0,
#         endTimeStep=endTs,
#         timeStepIncr=incTs_data
#     )
# )

# sim.createFieldSaver(
#     ParticleVectorFieldSaverPrms(
#         fieldName="force",
#         fileName=out_dat_dir + "/" + force_dir + "/force",
#         fileFormat="RAW2",
#         beginTimeStep=0,
#         endTimeStep=endTs,
#         timeStepIncr=incTs_data
#     )
# )

# sim.createFieldSaver(
#     ParticleVectorFieldSaverPrms(
#         fieldName="position",
#         fileName=out_dat_dir + "/" + position_dir + "/position",
#         fileFormat="RAW2",
#         beginTimeStep=0,
#         endTimeStep=endTs,
#         timeStepIncr=incTs_data
#     )
# )

# sim.createFieldSaver(
#     ParticleScalarFieldSaverPrms(
#         # fieldName="sigma_xx_2d",
#         fieldName="sigma_d",
#         fileName=out_dat_dir + "/" + sigmad_dir + "/sigma_d",
#         fileFormat="RAW_WITH_POS_ID",
#         beginTimeStep=0,
#         endTimeStep=endTs,
#         timeStepIncr=incTs_data
#     )
# )


# create a FieldSaver to wall forces:
force_saver = WallVectorFieldSaverPrms(
    wallName=["top", "bottom"],
    fieldName="Force",
    fileName=out_dat_dir+"/wall_Force.dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=endTs,
    timeStepIncr=incTs_data
)
sim.createFieldSaver(force_saver)


# create a FieldSaver to wall positions:
posn_saver = WallVectorFieldSaverPrms(
    wallName=["top", "bottom"],
    fieldName="Position",
    fileName=out_dat_dir+"/wall_Position.dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=endTs,
    timeStepIncr=incTs_data
)
sim.createFieldSaver(posn_saver)




# execute the simulation:
sim.run()
