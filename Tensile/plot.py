import sys, os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
from sys import platform

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

if platform == "darwin":
    matplotlib.use('Qt5Agg')
# import matplotlib.backends.backend_qt5agg
#matplotlib.use('Qt5Agg')
# rc('text', usetex=True)

plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['axes.linewidth'] = 0.7
# plt.rcParams['axes.grid'] = True
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['lines.markersize'] = 4.5
plt.rcParams['lines.markeredgewidth'] = 0.5
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.3
plt.rcParams["legend.fontsize"] = 7
# plt.rcParams["legend.markerscale"] = 1
# plt.rcParams["legend.fancybox"] = False
# plt.rcParams["legend.framealpha"] = 1
# plt.rcParams["legend.edgecolor"] = 'black'


class Data(object):
    """docstring for Data"""

    def __init__(self, path):
        super(Data, self).__init__()
        self.path = path

    def loadCondition(self):
        con_d = np.genfromtxt(
            self.path+"conditions.csv", delimiter=" ", dtype="U", autostrip=True)
        for c in con_d:
            c = c.split(",")
            print(c)
            if c[0] == "dT[s]":
                self.dT = float(c[1])
            if c[0] == "incT[s]":
                self.incT_data = float(c[1])
                self.incTs_data = int(self.incT_data/self.dT)
            if c[0] == "incT_snap[s]":
                self.incT_snap = float(c[1])
                self.incTs_snap = int(self.incT_snap/self.dT)
            if c[0] == "speed[m/s]":
                self.speed = float(c[1])
            if c[0] == "R[m]":
                self.R = float(c[1])
            if c[0] == "initialDistance[m]":
                self.iniD = float(c[1])
            if c[0] == "endstrain":
                t = self.iniD*float(c[1]) / (2.*self.speed)
                self.endTs = int(t / self.dT)
            if c[0] == "r_p[m]":
                self.r_p = float(c[1])
            if c[0] == "density[kg/m^3]":
                self.density = float(c[1])
            if c[0] == "maxTensile[Pa]":
                self.maxTensile = float(c[1])
            if c[0] == "Sh_mh":
                self.Sh = float(c[1])
            if c[0] == "AreaFacter":
                self.Afac = float(c[1])
            if c[0] == "bondModulus":
                self.bondModulus = float(c[1])
            if c[0] == "bondPoissonRatio":
                self.bondPoissonRatio = float(c[1])
            if c[0] == "bondCohesion":
                self.bondCohesion = float(c[1])
            if c[0] == "bondTanAngle":
                self.bondTanAngle = np.tan(np.radians(float(c[1])))
            if c[0] == "beta1":
                self.beta1 = float(c[1])
            if c[0] == "beta2":
                self.beta2 = float(c[1])
            if c[0] == "frictionCoeff":
                self.frictionCoeff = float(c[1])
            if c[0] == "normalK_sand[N/m]":
                self.normalK_sand = float(c[1])
            if c[0] == "kn/ks":
                self.alpha = float(c[1])
                self.shearK_sand = self.normalK_sand / self.alpha
            if c[0] == "dynamicMu":
                self.dynamicMu = float(c[1])
            if c[0] == "staticMu":
                self.staticMu = float(c[1])
            if c[0] == "viscosity_damp_nr":
                self.viscosity_damp_nr = float(c[1])
            if c[0] == "viscosity_damp_r":
                self.viscosity_damp_r = float(c[1])

    def loadData(self):
        dir_data = self.path + "out_data/"
        fn_con = self.path + "conditions.csv"

        # fn_ek = dir_data + "ekin.dat"
        # fn_ep = dir_data + "epot.dat"
        fn_nb = dir_data + "nbonds.dat"
        fn_of_comp = dir_data + "wall_Force.dat"
        fn_op_comp = dir_data + "wall_Position.dat"

        con_d = np.genfromtxt(
            fn_con, delimiter=" ", dtype="U", autostrip=True)

        # self.ekin = np.genfromtxt(fn_ek, delimiter=" ")
        # self.epit = np.genfromtxt(fn_ep, delimiter=" ")
        self.nbonds = np.genfromtxt(fn_nb, delimiter=" ")
        self.outF_comp = np.genfromtxt(fn_of_comp, delimiter=" ")
        self.outPos_comp = np.genfromtxt(fn_op_comp, delimiter=" ")


        con_l = []
        for c in con_d:
            c = c.split(",")
            print(c)
            if c[0] == "dT[s]":
                self.dT = float(c[1])
            if c[0] == "speed[m/s]":
                self.speed = float(c[1])*2
            if c[0] == "incT[s]":
                incT_data = float(c[1])
                self.tsInc= int(incT_data/self.dT)
            if c[0] == "H[m]":
                self.H = float(c[1])
            if c[0] == "R[m]":
                self.R = float(c[1])
            if c[0] == "endstrain":
                t = self.H*float(c[1]) / self.speed
                self.eTs = int(t / self.dT)
        self.sTs = 0
        self.csTs = 0
        self.cp = 0.0
        self.calPos()
        self.calF()
        self.calStrain()
        self.setT()
        self.cutT()


    def calF(self):
        self.area_bottom = self.R**2*np.pi

        self.F_top_comp = self.outF_comp[:, 2]
        self.F_bottom_comp = self.outF_comp[:, 6]
        self.P_top_comp = np.abs(self.F_top_comp/self.area_bottom)
        self.P_bottom_comp = np.abs(self.F_bottom_comp/self.area_bottom)

    def calPos(self):
        self.Pos_top_comp = self.outPos_comp[:, 2]
        self.Pos_bottom_comp = self.outPos_comp[:, 6]

    def calStrain(self):
        self.s1 = self.cp
        self.s3 = (self.P_top_comp+self.P_bottom_comp)*0.5
        self.stress = (self.s3 - self.s1)

        L0 = self.Pos_top_comp[0] - self.Pos_bottom_comp[0]
        dzt = abs(self.Pos_top_comp - self.Pos_top_comp[0])
        dzb = abs(self.Pos_bottom_comp - self.Pos_bottom_comp[0])
        self.strain = (dzt+dzb)/L0
        self.devStress_max = np.max(self.stress)

    def setT(self):
        self.nT = int((self.eTs)/self.tsInc + 1)
        self.time = np.zeros([self.nT])
        for it in range(0, self.nT):
            self.time[it] = (self.tsInc * it) * self.dT

    def cutT(self):
        N = np.shape(self.F_top_comp)[0]
        self.time = self.time[0:N]


    def loadData_for2particles(self):
        data = self._loadRAW_ig("./out_data/nforce", mode = "normal_force")
        self.position0 = []
        self.position1 = []
        self.nforce = []
        self.stress = []
        self.displacement = []
        self.strain = []
        for i in range(0, len(data)):
            self.position0.append(data[i][0])
            self.position1.append(data[i][1])
            self.nforce.append(-data[i][2])     #tensile +
            self.stress.append(-data[i][2]/(np.pi*self.r_p**2))
            print(-data[i][2]/(np.pi*self.r_p**2))
            self.displacement.append((data[i][0]-data[i][1])-(data[0][0]-data[0][1]))
            self.strain.append(self.displacement[i]/self.iniD)
        self.nbonds = np.genfromtxt("./out_data/nbonds.dat", delimiter=" ")
        self.time = [i * self.dT for i in range(0, len(data))] 

    def loadData_old(self):
        l_disp = self._loadRAW("./out_data/displacement", mode = "disp")
        l_force = self._loadRAW("./out_data/nforce", mode = "normal_force")
        # l_force = self._loadRAW("./out_data/force", mode = "force")
        l_position = self._loadRAW("./out_data/position", mode = "position")
        self.time = [i * self.dT for i in range(0, len(l_disp))] 

        l_disp = [l / self.iniD for l in l_disp]
        self.displacement = l_disp
        l_force = [l/ (np.pi*self.r_p**2) for l in l_force]
        self.stress = l_force
        print(len(l_position))
        print(len(l_position[0]))
        self.position0 = []
        self.position1 = []
        for i in range(0, len(l_position)):
            self.position0.append(l_position[i][0])
            self.position1.append(l_position[i][1])
        for i in range(0, len(self.position0)):
            self.displacement[i] = self.position0[i] - self.position1[1]
        # for i in range(0, len(self.stress)):
        #     print(i, self.displacement[i])
        #     print(i, self.stress[i])
        # print(l_disp)
        # print(l_force)

    def _loadRAW_ig(self, directory, mode):
        lfiles =sorted(os.listdir(directory)) # this doesn't sort correctly
        nfiles = len(lfiles)

        # print(lfiles)
        ld = []
        for i in range(0,len(lfiles)):
            if mode == "normal_force":
                fn = "nforce." + str(i) + ".dat"
                print("loading ", fn)
                data = np.genfromtxt(directory+"/"+fn, delimiter=" ")
                if len(data)==0:
                    p0 = np.nan
                    p1 = np.nan
                    d0 = np.nan
                    # print(lfiles[i])
                    # print(i, [p0, p1, d0])
                    ld.append([p0, p1, d0])
                else:
                    p0 = data[2]
                    p1 = data[6]
                    d0 = data[-1]
                    # print(lfiles[i])
                    # print(i, [p0, p1, d0])
                    ld.append([p0, p1, d0])
        return ld

    def _loadRAW(self, directory, mode):
        lfiles =sorted(os.listdir(directory))
        nfiles = len(lfiles)
        # print(lfiles)
        ld = []
        # self.time = [i * self.dT for i in range(0, nfiles)] 
        # print(self.time)
        for fn in lfiles:
            if mode == "normal_force":
                data = np.genfromtxt(directory+"/"+fn, delimiter=" ")
                d0 = data[-1]
                ld.append(d0)

            else:
                data = np.genfromtxt(directory+"/"+fn, delimiter=" ")
                d0 = data[0,-1]
                d1 = data[1,-1]
                if mode == "disp":
                    ld.append(d1-d0)
                elif mode == "force":
                    ld.append(d1)
                elif mode == "position":
                    ld.append([d0,d1])
        # print(ld)
        return ld
        
    def m2mm(self, d):
        dd = []
        for i in range(0, len(d)):
            dd.append(d[i] * 10**3)
        return dd

    def pa2MPa(self, d):
        dd = []
        for i in range(0, len(d)):
            dd.append(d[i] * 10**-6)
        return dd

    def plotNumBonds(self):
        fig = plt.figure(figsize=(3.14, 3.14), dpi=300,
                         facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.nbonds, '-')
        ax.set_ylabel(r"$\mathrm{\#\ of\ bonds}$")
        ax.set_xlabel(r"$\mathrm{time\ [s]}$")
        # ax.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("fig/numBonds_time.png", dpi=300)


    def plotDisp(self):
        fig = plt.figure(figsize=(3.14, 3.14), dpi=300,
                         facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.m2mm(self.displacement), '-')
        ax.set_ylabel(r"$\mathrm{displacement\ [mm]}$")
        ax.set_xlabel(r"$\mathrm{time\ [s]}$")
        # ax.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("fig/disp_time.png", dpi=300)
        
    def plotStress(self):
        fig = plt.figure(figsize=(3.14, 3.14), dpi=300,
                         facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.pa2MPa(self.stress), '-')
        ax.set_ylabel(r"$\mathrm{stress\ [MPa]}$")
        ax.set_xlabel(r"$\mathrm{time\ [s]}$")
        # ax.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("fig/stress_time.png", dpi=300)

    def plotNormalForce(self):
        fig = plt.figure(figsize=(3.14, 3.14), dpi=300,
                         facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.nforce, '-')
        # ax.plot(self.time, self.position1, '-')
        ax.set_ylabel(r"$\mathrm{normal\ force\ [N]}$")
        ax.set_xlabel(r"$\mathrm{time\ [s]}$")
        # ax.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("fig/nforce_time.png", dpi=300)

    def plotPosition(self):
        fig = plt.figure(figsize=(3.14, 3.14), dpi=300,
                         facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.m2mm(self.position0), '-')
        # ax.plot(self.time, self.position1, '-')
        ax.set_ylabel(r"$\mathrm{position\ [mm]}$")
        ax.set_xlabel(r"$\mathrm{time\ [s]}$")
        # ax.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("fig/position_time.png", dpi=300)


    def plotStrain(self):
        fig = plt.figure(figsize=(3.14, 3.14), dpi=300,
                         facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.strain*100, self.pa2MPa(self.stress), '-')
        # ax.plot(self.strain, self.stress, '-')
        ax.set_ylabel(r"$\mathrm{stress\ [MPa]}$")
        ax.set_xlabel(r"$\mathrm{strain\ [\%]}$")
        # ax.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("fig/stress_strain.png", dpi=300)


def main():
    dl = []

    d = Data("./")
    dl.append(d)

    for d in dl:
        # d.loadCondition()
        d.loadData()
        # d.plotPosition()
        d.plotNumBonds()
        # d.plotDisp()
        # d.plotStress()
        # d.plotNormalForce()
        d.plotStrain()

if __name__ == "__main__":
    main()
