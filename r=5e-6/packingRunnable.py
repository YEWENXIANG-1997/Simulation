# https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
from __future__ import division
from esys.lsm import *
from esys.lsm.util import *
import time
import random
import math
import numpy as np


class packingRunnable(Runnable):
    def __init__(self,
                 LsmMpi=None,
                 porosity=0.4,
                 R=1e-3,
                 H=4e-3
                 ):
        """
        Subroutine to initialise the Runnable and store parameter values.
        """
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.dt = self.sim.getTimeStepSize()
        self.tporo = porosity
        self.poro = 1.0
        self.ts = 0
        self.R = R
        self.H = H
        self.V_all = np.pi * self.R**2 * self.H
        self.fac = 1.0
        self.endTs = 10000000
        self.addTs = 10000
        self.bswitch = True

    def run(self):
        self.ts = self.sim.getTimeStep()
        # self.sim.getParticleList()[0].setTag(self.ts)
        # print(self.sim.getParticleList()[0].getTag())
        if self.ts%100 == 0:
            print("________________________")
            print("ts              :{:>8}".format(self.ts))
            print("num of particles:{:>8}".format(len(self.sim.getParticleList())))
            print("target porosity :{:>8}".format(self.tporo))
            print("current porosity:{:>8}".format(self.poro))
            if self.poro > self.tporo:
                print("current status  :{:>8}".format("scaling"))
            else:
                print("current status  :{:>8}".format("easing"))

        if self.poro > self.tporo:
            self.scaling()
        else:
            self.plist = self.sim.getParticleList()

            if self.bswitch is True:
                self.setEndTs()
                self.bswitch = False

        if self.ts == self.endTs:
            self.sim.exit()

    def scaling(self):
        self.plist = self.sim.getParticleList()
        # print("nump", len(self.plist))
        # self.fac += 0.01
        beta = 0.3
        gamma = 1
        self.fac = 1.0 + beta / self.ts**gamma
        self.sim.setParticleRadiusFactor(self.fac)
        V = 0.0
        for p in self.plist:
            r = p.getRad()
            #print(self.fac)
            #print(r)
            V += 4.*np.pi*r**3/3.0
        self.poro = (self.V_all - V)/self.V_all

    def setEndTs(self):
        self.entTs = self.ts + self.addTs
