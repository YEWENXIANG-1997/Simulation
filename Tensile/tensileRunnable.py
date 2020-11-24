# import the division module for compatibility between Python 2 and Python 3
from __future__ import division
# import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *

import os
import sys
import time
import numpy as np

class tensileClass(Runnable):
    def __init__(self,
                 LsmMpi=None,
                 speed=1e-5,
                 endTs=1000,
                 iniD=1e-3
                 ):
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.speed = speed
        self.endTs = endTs
        self.iniD = iniD
        self.dt = self.sim.getTimeStepSize()
        self.dz = self.dt*self.speed
        self.Lini = self.iniD
        self.p_bottom = self.sim.getWallPosition("bottom").toList()[2]
        self.p_top = self.sim.getWallPosition("top").toList()[2]
        self.strain = 0.0

    def run(self):
        self.ts = self.sim.getTimeStep()
        self.p_bottom = self.sim.getWallPosition("bottom").toList()[2]
        self.p_top = self.sim.getWallPosition("top").toList()[2]
        Lnow = self.p_top - self.p_bottom
        self.strain = (Lnow - self.Lini) / (self.Lini) * 100
        self.sim.moveWallBy("top", Vec3(0, 0, +self.dz))
        self.sim.moveWallBy("bottom", Vec3(0, 0, -self.dz))
        if self.strain > 100.0:
            self.sim.exit()

        if self.ts % 100 == 0:
            print("________________________")
            print("ts            :{:>8}".format(self.ts))
            print("end ts        :{:>8}".format(self.endTs))
            print("dz         [m]:{:>8}".format(self.dz))
            print("strain rate[%]:{:>8}".format(self.strain))
        

