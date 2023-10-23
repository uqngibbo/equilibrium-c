"""
Application script: Compute a normal shock with ionization

v2: With pressure error method.

@author: Nick Gibbons
"""

from numpy import zeros, array, linspace, meshgrid, sqrt, floor, ceil, log10
from collections import namedtuple
from scipy.optimize import newton 
import matplotlib.pyplot as plt
from matplotlib import colors, cm, rcParams
from atmo import US76Atmo
from normalshock import *
import pyeq

plt.rcParams.update({'font.size': 12})
plt.rcParams['svg.fonttype'] = 'none'

spnames = ['N2', 'N2+', 'NO', 'NO+', 'O2', 'O2+', 'N', 'N+', 'O', 'O+', 'e-']
X0 = array([0.77, 0.0, 0.0, 0.0, 0.23, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
X0/=(X0.sum())
ceq = pyeq.EqCalculator(spnames)

atmo = US76Atmo()

hhs = linspace(10e3, 60e3, 50)
vvs = linspace(2e3, 8e3, 50)

vs,hs = meshgrid(vvs, hhs)

ends= zeros(hs.shape)
ifs = zeros(hs.shape)
Ms  = zeros(hs.shape)
rhos= zeros(hs.shape)
ps  = zeros(hs.shape)
Ts  = zeros(hs.shape)
pvs = zeros(hs.shape)
Xs  = zeros((hs.shape[0], hs.shape[1], len(spnames)))

ni,nj = hs.shape
lowest = 1e99

ii = 0
for i in range(ni):
    for j in range(nj):
        h = hs[i,j]
        V = vs[i,j]
        rho0 = atmo.rho(h)
        T0   = atmo.T(h)
        p0 = rho0*287.0*T0
        M = V/sqrt(1.4*287.0*T0)

        preshock = GasState.from_pTv(p=p0, T=T0, v=V, X0=X0, ceq=ceq)
        postshock = normal_shock(preshock)

        ionisation_fraction = postshock.X[spnames.index('e-')]*2

        Nav = 6.02214076e+23
        Y = ceq.XtoY(postshock.X)
        eidx = spnames.index('e-')
        Me = ceq.M[eidx]
        rhoe = postshock.rho*Y[eidx]
        ne = rhoe*Nav/Me

        #print("ii={} v={} h={} p={} T={} -> ne: {:e}".format(ii, V, h, postshock.p, postshock.T, ne))

        if (ne!=0.0): lowest = min(lowest, log10(ne))
        ends[i,j] = log10(max(ne,1e10))
        ifs[i,j] = log10(max(ionisation_fraction,1e-10))
        Ms[i,j] = M
        rhos[i,j] = postshock.rho
        ps[i,j] = postshock.p
        Ts[i,j] = postshock.T
        pvs[i,j] = postshock.v
        Xs[i,j,:] = postshock.X
        ii += 1

thing = ends
print("Lowest ne was :", lowest)
print("min: ", thing.min(), "max: ", thing.max())
minlevel = floor(thing.min())
maxlevel = floor(thing.max())+1
nlevels = int(round((maxlevel-minlevel), 0))
print("minlevel: ", minlevel, " maxlevel: ", maxlevel, "nlevels: ", nlevels)
levels = linspace(minlevel, maxlevel, nlevels+1)

fig = plt.figure(figsize=((8.5,5)))
axes = fig.subplots(1,1)
cs = axes.contourf(vs,hs/1e3,thing,cmap=cm.autumn, levels=levels)

axes.set_ylabel('Altitude (km)')
axes.set_xlabel('Velocity (m/s)')

nlevels = int(round((maxlevel-minlevel), 0))
clevels = linspace(minlevel, maxlevel, nlevels+1)
cs2 = axes.contour(vs,hs/1e3,thing,levels=clevels, colors='black')
axes.clabel(cs2, inline=1, fontsize=10, fmt=" %d ")
#axes.clabel(cs2, inline=1, inline_spacing=10, fontsize=12, fmt=lambda s : "$10^{:d}$".format(int(s)))

#cs2= [
#    axes.contour(vs,hs/1e3,Ms,levels=list(range(8,26,4)), colors='black', linestyles='dotted'),
#]
#axes.clabel(cs2[0], inline=1, fontsize=10, inline_spacing=10, fmt=" %d ")


#axes.set_title("Postshock Ionisation Map: Contours of Mach Number")
axes.set_title("Postshock Ionisation Map: Electron Number Density")

vmin = thing.min()
vmax = thing.max()

norm = colors.Normalize(vmin=vmin, vmax=vmax)
cs.set_norm(norm)
#fig.colorbar(cs[1], ax=axes, shrink=0.9, orientation='horizontal')
#cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.spring),
#             ax=axes, orientation='horizontal', label='Ionisation Fraction',
#             format=FuncFormatter(lambda s, pos : "$10^{:d}$".format(int(s))))
cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.autumn),
             ax=axes, orientation='horizontal', label='log10(Electron Number Density)')
             #format=FuncFormatter(lambda s, pos : "$10^{:d}$".format(int(s))))
cbar.set_ticks(levels)
plt.tight_layout()
plt.show()

