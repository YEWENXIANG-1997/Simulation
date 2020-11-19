from scipy.optimize import curve_fit
from scipy.special import erf
from scipy import stats
import numpy as np


class Pack(object):
    def __init__(self):
        self.R = 4.0
        self.H = 18.0
        self.V_all = np.pi*self.R**2*self.H
        self.tporo = 0.4
        self.rvs_size = 1000
        self.mode = 0.21
        self.median = 0.22
        self.d0 = 0.05
        self.d100 = 0.3
        self.mu = np.log(self.median)
        self.sigma = (np.log(self.median) - np.log(self.mode))**0.5
        self.xl = []
        self.yl = []
        self.zl = []
        self.rl = []
        self.inifac = 0.1

    def setCdfParm(self, path):
        data = np.genfromtxt(path, delimiter=",")
        self.mode = data[0]
        self.median = data[1]
        self.mu = np.log(self.median)
        self.sigma = (np.log(self.median) - np.log(self.mode))**0.5
        print("loaded ", path)
        print("mode:", self.mode)
        print("median:", self.median)
        print("mu:", self.mu)
        print("sigam:", self.sigma)

    def setMaxMinDiam(self, d0, d100):
        self.d0 = d0
        self.d100 = d100

    def setConstDiam(self, d):
        self.d_const = d

    def setInitFac(self, inifac):
        self.inifac = inifac

    def setCylinder(self, R, H):
        self.R = R
        self.H = H
        self.V_all = np.pi*self.R**2*self.H

    def setTargetPorosity(self, poro):
        self.tporo = poro

    def setMedianMode(self, median=0.22, mode=0.21):
        self.mode = mode
        self.median = median
        self.mu = np.log(self.median)
        self.sigma = (np.log(self.median) - np.log(self.mode))**0.5

    def makeParticle_constD(self):
        V = 0.0
        count = 0
        for i in range(1000000):
            r = 0.5*self.d_const
            v = 4.*np.pi*r**3/3.0
            V += v
            self.rl.append(r*1e-3)
            rp = self.R*np.random.rand()*1e-3
            theta = 2.0*np.random.rand()*np.pi
            self.xl.append(rp*np.cos(theta))
            self.yl.append(rp*np.sin(theta))
            self.zl.append(self.H*np.random.rand()*1e-3)
            count += 1
            poro = (self.V_all-V)/self.V_all
            print(poro*100)
            if poro < self.tporo:
                print(count)
                break

    def makeParticle(self):
        V = 0.0
        count = 0
        for i in range(1000000):
            dl = stats.lognorm(
                s=self.sigma, loc=0, scale=np.exp(self.mu)).rvs(size=self.rvs_size)
            for d in dl:
                if d < self.d0 or d > self.d100:
                    continue
                    # break
                r = 0.5*d
                v = 4.*np.pi*r**3/3.0
                V += v
                self.rl.append(r*1e-3)
                rp = self.R*np.random.rand()*1e-3
                theta = 2.0*np.random.rand()*np.pi
                self.xl.append(rp*np.cos(theta))
                self.yl.append(rp*np.sin(theta))
                self.zl.append(self.H*np.random.rand()*1e-3)
                count += 1
            poro = (self.V_all-V)/self.V_all
            print(poro*100)
            if poro < self.tporo:
                print(count*self.rvs_size)
                break

    def getPXList(self):
        return self.xl

    def getPYList(self):
        return self.yl

    def getPZList(self):
        return self.zl

    def getPRList(self):
        return self.rl

    def makeGeofile(self):
        lines = []
        ep = self.H*1e-3*0.01
        r = self.R*1e-3
        h = self.H*1e-3
        s = "LSMGeometry 1.2"
        lines.append(s)
        s = "BoundingBox {} {} {} {} {} {}".format(
            -r-ep, -r-ep, -ep, r+ep, r+ep, h+ep)
        lines.append(s)
        s = "PeriodicBoundaries 0 0 0"
        lines.append(s)
        s = "Dimension 3D"
        lines.append(s)
        s = "BeginParticles"
        lines.append(s)
        s = "Simple"
        lines.append(s)
        s = str(len(self.rl))
        lines.append(s)
        for i in range(0, len(self.rl)):
            s = "{} {} {} {} {} -1".format(
                self.xl[i], self.yl[i], self.zl[i], self.rl[i]*self.inifac, i)
            lines.append(s)
        s = "EndParticles"
        lines.append(s)
        s = "BeginConnect"
        lines.append(s)
        s = "0"
        lines.append(s)
        s = "EndConnect"
        lines.append(s)
        with open("ini.geo", mode='w') as f:
            f.write('\n'.join(lines))
