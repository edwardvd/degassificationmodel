# -*- coding: utf-8 -*-
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as mpl
import numpy
"""
Created on Tue Jan 23 12:12:11 2018

@author: NL14168
"""

"""
General setup
"""

mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 14


"""
Define reactor

    L: reactor height (m)
    J: number of cells
    j: cell index number, with 0 = top and J is bottom of the reactor
    T: total time to run (s)
    N: number of timesteps
"""

"""
Initialize model¶
"""


"""
Setup model type
"""
# logistic function to make a smooth step function


def logistic(x, L, k, x0):
    return numpy.power(numpy.exp(-k*(x-x0))+1, -1) * L


class ADM():
    def __init__(self):
        self.Dax = 0.
        self.a = 0.
        self.aklaref = 0.
        self.akla = 0.

        self.Tc = 13.5
        self.Patm = 1.

        # Constants
        self.g = 9.81  # gravity
        self.rho = 1000  # density

        '''
        Constants for Henry's law
        '''
        self.frn = 0.79         # [mol/mol] gas fraction N2 gas in air
        self.H = 1600.          # [L.atm/mol] Henry coefficient N2 gas
        self.Hc = 1300.         # [K] temperature compensation for Henry coefficient N2 gas
        self.M = 28.            # [g/mol] molair weight of N2 gas

        self.DN = 1.88E-5
        self.DO = 2.10E-5
        self.Tref = 298.15

    def f_vec(self):
        assert False

    def dim(self, L, J, T, dt):
        N = T/dt + 1
        print N
        self.L, self.J, self.T, self.N = L, J, T, N

    def redim(self):
        self.H2 = self.H/self.M           # [L.atm/g] Henry coefficient N2 gas on weight base
        self.Tk = 273.15 + self.Tc

        self.dx = float(self.L)/float(self.J-1)
        x_grid = numpy.array([j*self.dx for j in range(self.J)])

        self.x_grid = x_grid
        # initialize timestep
        self.dt = float(self.T)/float(self.N-1)
        self.t_grid = numpy.array([n*self.dt for n in range(self.N)])

        # pressure
        self.P = self.x_grid * self.g * self.rho/1E5 + self.Patm  # water pressure
        self.Pp = self.P * self.frn  # partial pressure

        # concentrations
        print "REDIM Tk %f" % self.Tk
        self.cs = self.P/self.fH(self.Tk)*1000.
        self.ca = self.Pp/self.fH(self.Tk)*1000.
        self.cb = None

    def setup(self):
        self.redim()
        sigma = float(self.Dax*self.dt)/float(2.*self.dx*self.dx)
        rho = float(self.a*self.dt)/float(4.*self.dx)
        J = self.J

        self.A = numpy.diagflat([-sigma+rho for i in range(J-1)], -1) +\
            numpy.diagflat([1.+sigma+rho]+[1.+2.*sigma for i in range(J-2)]+[1.+sigma-rho]) +\
            numpy.diagflat([-(sigma+rho) for i in range(J-1)], 1)

        self.B = numpy.diagflat([sigma-rho for i in range(J-1)], -1) +\
            numpy.diagflat([1.-sigma-rho]+[1.-2.*sigma for i in range(J-2)]+[1.-sigma+rho]) +\
            numpy.diagflat([sigma+rho for i in range(J-1)], 1)

    def run(self):
        self.setup()
        C_record = []
        t_record = []
        self.C_total = []
        self.C_total.append(self.C)
        for ti in range(1, self.N):
            if not self.cb is None:
                self.C[-1] = self.cb  # set boundary
            C_new = numpy.linalg.solve(self.A, self.B.dot(self.C) + self.f_vec())
            self.C = C_new
            t_record.append(self.t_grid[ti])
            C_record.append(self.C[-1])
            self.C_total.append(self.C)

    # function for temperature correction Henry coefficient
    def fH(self, Tk):
        return self.H2 * numpy.exp(-self.Hc * (1/Tk - 1/(298.15)))

    # function for temperature correction of Diffusion coefficient
    def fD(self, D, Tk):
        return D * Tk / self.Tref * self.mu(self.Tref)/self.mu(Tk)

    def mu(self, Tk):
        return 2.414E-5 * 10.**(247.8/(Tk-140.))

    # denitrification
    def Xrate(self, Tc, X):
        return logistic(Tc, 0.5, 0.5, 11)*X/3600.

    def plot(self, tarray):
        C_total = numpy.array(self.C_total)
        t_grid = self.t_grid
        x_grid = self.x_grid
        X_grid = self.X_grid
        for tval in tarray:

            it = numpy.argwhere(t_grid == tval)[0][0]

            # mpl.figure()

            fig, ax1 = mpl.subplots()
            lns = ax1.plot(x_grid, self.cs, 'k--')
            lns += ax1.plot(x_grid, C_total[it, :], 'k')
            lns += ax1.plot(x_grid, self.ca, 'k-.')

            ax2 = ax1.twinx()
            lns += ax2.plot(x_grid, X_grid, color='tab:gray')
            ax1.set_xlabel("reactor depth [m]")
            ax1.set_ylabel("concentration [mg L$^{-1}$]")
            ax2.set_ylabel("biomass concentration [g L$^{-1}$]")
            mpl.title("%d min" % (tval/60))

            mpl.legend(lns, ("saturation", "actual", "equilibrium", "biomass"), loc=0)
            #mpl.legend(("saturation", "actual", "equilibrium", "biomass"))

            minorLocator = MultipleLocator(0.1)
            mpl.gca().xaxis.set_minor_locator(minorLocator)
            mpl.gca().xaxis.set_tick_params(which='minor', direction='in')

            minorLocator = MultipleLocator(0.5)
            mpl.gca().yaxis.set_minor_locator(minorLocator)
            mpl.gca().yaxis.set_tick_params(which='minor', direction='in')
            mpl.savefig("feeding_withstrip%04d.eps" % tval)
        # mpl.show()

    def plot2(self, tarray):
        C_total = numpy.array(self.C_total)
        t_grid = self.t_grid
        x_grid = self.x_grid
        X_grid = self.X_grid
        for tval in tarray:

            it = numpy.argwhere(t_grid == tval)[0][0]

            # mpl.figure()

            fig, ax1 = mpl.subplots()
            lns = ax1.plot(self.cs, x_grid, 'k--')
            lns += ax1.plot(C_total[it, :], x_grid, 'k')
            lns += ax1.plot(self.ca, x_grid, 'k-.')
            ax1.set_ylabel("reactor depth [m]")
            ax1.set_xlabel("concentration [mg L$^{-1}$]")
            ax1.invert_yaxis()

            if False:
                ax2 = ax1.twinx()
                lns += ax2.plot(x_grid, X_grid, color='tab:gray')
                ax2.set_ylabel("biomass concentration [g L$^{-1}$]")

            mpl.title("%d min" % (tval/60))

            mpl.legend(lns, ("saturation", "actual", "equilibrium", "biomass"), loc=0)
            #mpl.legend(("saturation", "actual", "equilibrium", "biomass"))

            minorLocator = MultipleLocator(0.1)
            mpl.gca().xaxis.set_minor_locator(minorLocator)
            mpl.gca().xaxis.set_tick_params(which='minor', direction='in')

            minorLocator = MultipleLocator(0.5)
            mpl.gca().yaxis.set_minor_locator(minorLocator)
            mpl.gca().yaxis.set_tick_params(which='minor', direction='in')
            mpl.savefig("feeding_withstrip%04d.eps" % tval)

    def plot3(self, tarray):
        C_total = numpy.array(self.C_total)
        t_grid = self.t_grid
        x_grid = self.x_grid
        X_grid = self.X_grid
        mpl.rcParams['figure.figsize'] = 12, 10
        label = ("A", "B", "C", "D", "E", "F")
        ilabel = 0
        fig, ax = mpl.subplots(2, 2, sharex=True, sharey=True)

        fig.subplots_adjust(hspace=0.1, wspace=0.05,
                            left=0.1, right=0.95,
                            top=0.95, bottom=0.1)

        for tval in tarray:
            ax1 = ax[ilabel/2, ilabel % 2]
            it = numpy.argwhere(t_grid == tval)[0][0]

            # mpl.figure()

            lns = []
            #lns = ax1.plot(x_grid, self.cs, 'k--')
            ymin = numpy.min(self.ca)
            ymax = numpy.max(self.cs)
            ymin = 16
            ymax = 36
            # ax1.fill_between(x_grid, ymax, self.cs,facecolor='#e5e5e5')
            # ax1.fill_between(x_grid, ymin, self.ca,facecolor='#e7f0f9')
            ax1.fill_between(x_grid, ymax, self.cs, facecolor='#CC0000', hatch='')
            # ax1.fill_between(x_grid, ymax, self.cs,facecolor='None', hatch='x', edgecolor = "#CC0000")
            ax1.fill_between(x_grid, ymin, self.ca, facecolor='#66CC00', hatch='')
            # ax1.fill_between(x_grid, ymin, self.ca, facecolor="none", hatch="X",
            #                 edgecolor="b", linewidth=0.0)
            lns += ax1.plot(x_grid, C_total[it, :], 'k')
            #lns += ax1.plot(x_grid,self.ca, 'k-.')

            # t =ax1.text(0, ymax, "\n  (%s)\n  over saturated" % label[ilabel],
            #         fontsize=16,verticalalignment='top')
            ax1.text(0, ymax, "\n over saturated",
                     fontsize=16, verticalalignment='top')
            t = ax1.text(0.2, ymin+0.6, "%s" % label[ilabel],
                         fontsize=16, verticalalignment='bottom')
            t.set_bbox(dict(facecolor='white', alpha=0.75, edgecolor='black'))

            ax1.text(7, ymin, "max deficit  \n",
                     fontsize=16, verticalalignment='bottom', horizontalalignment='right')
            mpl.grid(axis='None')

            # mpl.ylim((numpy.min(self.ca),numpy.max(self.cs)))
            mpl.ylim((ymin, ymax))
            ax1.yaxis.set_major_locator(MultipleLocator(2))

            if (ilabel/2 == 1):
                ax1.set_xlabel("reactor depth [m]")
            if (ilabel % 2 == 0):
                ax1.set_ylabel("concentration [mg L$^{-1}$]")

            if False:
                ax2 = ax1.twinx()
                lns += ax2.plot(x_grid, X_grid, color='tab:gray')

                if (ilabel % 2 == 1):
                    ax2.set_ylabel("biomass concentration [g L$^{-1}$]")
                #mpl.title("%d min" % (tval/60))
                ax1.legend(lns, ("N$_2$ concentration", "MLSS concentration"), loc=1)
                ax2.set_ylim(0, 25)
            else:
                ax1.legend(lns, ("N$_2$ concentration",), loc=1)
            ax1.set_xlim((0, 7))
            ax1.set_ylim()
            if False:

                #mpl.legend(("saturation", "actual", "equilibrium", "biomass"))

                minorLocator = MultipleLocator(0.1)
                mpl.gca().xaxis.set_minor_locator(minorLocator)
                mpl.gca().xaxis.set_tick_params(which='minor', direction='in')

                minorLocator = MultipleLocator(0.5)
                mpl.gca().yaxis.set_minor_locator(minorLocator)
                mpl.gca().yaxis.set_tick_params(which='minor', direction='in')

            mpl.savefig("Feeding%04d.png" % tval)
            ilabel += 1


class ADM_Strip(ADM):
    def f_vec(self):
        # gas transfer: normal henry + additional 10 times factor for oversaturation
        return (-self.akla * (self.C-self.ca) + -self.akla*10*numpy.max(self.C-self.cs, 0)) * self.dt

    def redim(self):
        ADM.redim(self)
        self.C = self.cs
        self.X_grid = self.x_grid*0. + 10.
        self.akla = self.aklaref * numpy.sqrt(self.fD(self.DN, self.Tk)/self.fD(self.DO, 293.15))
        print "akla: %f" % (self.akla * 3600.)
        print numpy.sqrt(self.fD(self.DO, self.Tk)/self.fD(self.DO, 293.15))


class ADM_Feed(ADM):
    def f_vec(self):
        return (self.Xrate(self.Tc, self.X_grid)) * self.dt

    def redim(self):
        ADM.redim(self)
        self.X_grid = logistic(self.x_grid, 10.*7/4, 10, 3)
        #self.C = self.ca
        self.cb = self.frn/self.fH(self.Tk)*1000.


print "Running stripping"
model1 = ADM_Strip()
model1.Tc = 13.
model1.dim(L=7., J=71, T=5*60, dt=5)
model1.Dax = 0.01
model1.a = 0.
# TODO temperature correction
model1.aklaref = 5./3600.
model1.run()
# model1.plot3((0,300,900,3600,))

print "Running feeding"
model = ADM_Feed()
model.Tc = 13.
model.dim(L=7., J=71, T=7200, dt=5)
model.Dax = 0.0001
model.a = 4./3600
model.C = model1.C
model.run()
model.plot3((0, 900, 1800, 3000))

mpl.show()


print "Done!"
