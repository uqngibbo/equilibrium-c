"""
Handbook of Geophysics Atmosphere Model from 1985

Created on Tue Jul 23 14:23:29 2013

Originally for MECH4800 Space Engineering Assignment One
 - Updated to python3 in March 2021

@author: Nick N. Gibbons
"""

from numpy import array, arange, ceil, floor, int32, linspace, sqrt, exp
from scipy.interpolate import interp1d
from pylab import plot,show

class WMOJanuaryAtmo(object):
    """
    Dynamic Container for interpolation of atmospheric data.
    Uses Reference Data at 0 Latitude in January average conditions.

    References:
    -Champion K.S.W. et. al. "Handbook of Geophysics: Chapter 14 Standard and 
    Reference Atmospheres" 
    Online at http://www.cnofs.org/Handbook_of_Geophysics_1985/Chptr14.pdf
    Accessed on 23/7
    """
    def __init__(self,):
        
        h = [2000.*i for i in range(46)] # m
        
        T = [299.15,288.78,278.41,268.06,256.40,240.49,224.58,208.69,197.34,
             195.03,201.37,207.71,214.05,220.38,226.71,232.70,236.65,240.60,
             244.55,249.66,254.98,260.30,265.61,270.92 ,272.15 ,272.15 ,269.12,
             265.01 ,260.71 ,255.61 ,250.52 ,245.44 ,240.35 ,235.27 ,230.19 ,
             225.12 ,217.25, 206.72 ,196.20 ,195.65 ,195.65 ,195.65 ,195.65 ,
             194.59 ,191.29 ,187.99]    # K
        
        p = [1.0100e3,8.0106e2,6.3009e2,4.9119e2,3.7918e2,2.8842e2,2.1535e2,
             1.5741e2,1.1258e2,7.9414e1,5.6424e1,4.0525e1,2.9403e1,2.1538e1,
             1.5920e1,1.1867e1,8.9018 ,6.7103 ,5.0825,3.8695,2.9634,2.2824,
             1.7675,1.3759,1.0749,8.4007e-1,6.5594e-1,5.1036e-1,3.9560e-1,
             3.0525e-1,2.3434e-1,1.7897e-1,1.3593e-1,1.0265e-1,7.7069e-2,
             5.7500e-2,4.2580e-2,3.1098e-2,2.2347e-2,1.5906e-2,1.1324e-2,
             8.0642e-3,5.7436e-3,4.0905e-3,2.9008e-3,2.0453e-3] # mBar
             
        rho = [1.1761,9.6636e-1,7.8839e-1,6.3834e-1,5.1518e-1,4.1779e-1,
               3.3404e-1,2.6276e-1,1.9873e-1,1.4185e-1,0.9761e-1,6.7966e-2,
               4.7853e-2,3.4046e-2,2.4462e-2,1.7766e-2,1.3103e-2,0.9715e-2,
               7.2401e-3,5.3993e-3,4.0487e-3,3.0546e-3,2.3182e-3,1.7692e-3,
               1.3760e-3,1.0753e-3,8.4906e-4 ,6.7090e-4,5.2860e-4,4.1600e-4,
               3.2586e-4,2.5402e-4,1.9702e-4,1.5200e-4,1.1663e-4,8.8979e-5,
               6.8278e-5,5.2406e-5,3.9678e-5,2.8323e-5,2.0164e-5,1.4358e-5,
               1.0226e-5,7.3229e-6,5.2827e-6,3.7900e-6] # kg/m3
               
        self.T = interp1d(h,T)
        p_kPa = [i/10. for i in p]    #convert to kPa
        self.p = interp1d(h,p_kPa)
        self.rho = interp1d(h,rho)
        return

class US76Atmo(object):
    """
    US Standard Atmosphere, 1976 version. 

    References:
      - https://en.wikipedia.org/wiki/U.S._Standard_Atmosphere (Accessed 21/06/24)
      - "U.S. Standard Atmosphere, 1976", NASA-TM-X-74335, 1976

    """
    headings = ['b', 'Geopotential Height (m)',  'Geopotential Height (ft)', 'Static Pressure (Pa)',
                'Static Pressure (inHg)', 'Standard Temperature (K)',
                'Temperature Lapse Rate (K/m)',   'Temperature Lapse Rate (K/ft)']

    data = """
        0   0   0   101325  29.92126  288.15  -0.0065   -0.001981
        1   11000  36089  22632.1   6.683245  216.65  0.0   0.0
        2   20000  65617  5474.89   1.616734  216.65  0.001   0.0003048
        3   32000  104987   868.019   0.2563258   228.65  0.0028  0.0008534
        4   47000  154199   110.906   0.0327506   270.65  0.0   0.0
        5   51000  167323   66.9389   0.01976704  270.65  -0.0028   -0.0008534
        6   71000  232940   3.95642   0.00116833  214.65  -0.002  -0.0006096
        6   86000  282152   0.30233   0.00008928  184.65   0.000  0.0"""

    # TODO: The last row of the table is made up except for the altitudes
    # Also: we might want to use the exact pressures computed by the previous
    # section instead of the tabulated ones.
    R = 8.31432/28.9644e-3 # USAA1976 uses a nominal dry air weight
    r0 = 6356766 # Nominal earth radius from page 8 of the Standard.
    g0 = 9.80665 # nominal standard gravitational acceleration

    def __init__(self):
        data = list(filter(bool, __class__.data.splitlines()))
        data = [list(map(float,line.strip().split())) for line in data]

        self.sections = []
        for d in data:
            self.sections.append(dict(zip(__class__.headings, d)))

        self.geopotentials = []
        for section in self.sections:
            section['k'] = self.compute_k(section)
            self.geopotentials.append(section['Geopotential Height (m)'])

        self.geopotentials = array(self.geopotentials)
        idxs = arange(self.geopotentials.size)
        self.interpidx_from_h = interp1d(self.geopotentials, idxs)

    def check_limits_h(self, h):
        if h>86000.0:
            raise Exception("Error requested altitude is above model validity limits.")
        if h<0.0:
            raise Exception("Error requested altitude is below model validity limits.")

    def loweridx_from_h(self, h):
        self.check_limits_h(h)
        return floor(self.interpidx_from_h(h)).astype(int32)

    def upperidx_from_h(self, h):
        self.check_limits_h(h)
        return ceil(self.interpidx_from_h(h)).astype(int32)

    def test(self):
        for section0, section1 in zip(self.sections[:-1], self.sections[1:]):
            print("Section: ", section0['b'])
            print("Table: p0 {} T0 {}: ".format(section0['Static Pressure (Pa)'], section0['Standard Temperature (K)']))
            print("Computed p0 {} T0: {}".format(self.p_from_section(section0, section0['Geopotential Height (m)']), self.T_from_section(section0, section0['Geopotential Height (m)'])))
            print("Computed rho0: {}".format(self.rho(section0['Geopotential Height (m)'])))

            print("Table: p1 {} T1 {}: ".format(section1['Static Pressure (Pa)'], section1['Standard Temperature (K)']))
            print("Computed p1 {} T1: {}".format(self.p_from_section(section0, section1['Geopotential Height (m)']), self.T_from_section(section0, section1['Geopotential Height (m)'])))
            print(" ")

    def compute_k(self, section):
        """
        The separable differential equation for pressure from Temperature and altitude has a 
        constant k, that requires the pressure at some altitude to solve.

        See: https://www.wolframalpha.com/input/?i=dy%2Fdx+%3D+-a*y%2F%28b%2Bc*x%29
        Annoyingly, the L==0 solution is very different to the L!=0 one.
        """
        p0 = section['Static Pressure (Pa)']
        T0 = section['Standard Temperature (K)']
        L = section['Temperature Lapse Rate (K/m)']
        h0 = section['Geopotential Height (m)']
        R = __class__.R
        g = __class__.g0
        if L==0.0:
            k = p0/exp(-g*h0/R/T0)
        else:
            k = p0/(R*T0)**(-g/R/L)
        return k

    def altitude_from_geopotential(self, h):
        z = self.r0*h/(self.r0 - h)
        return z

    def geopotential_from_altitude(self, z):
        h = self.r0*z/(self.r0 + z)
        return z
        
    def p(self, z):
        h = self.geopotential_from_altitude(z)
        l = self.loweridx_from_h(h)
        section = self.sections[l]
        return self.p_from_section(section, h)

    def T(self, z):
        h = self.geopotential_from_altitude(z)
        l = self.loweridx_from_h(h)
        section = self.sections[l]
        return self.T_from_section(section, h)

    def rho(self, z):
        h = self.geopotential_from_altitude(z)
        l = self.loweridx_from_h(h)
        section = self.sections[l]
        T = self.T_from_section(section, h)
        p = self.p_from_section(section, h)
        rho = p/__class__.R/T
        return rho

    def p_from_section(self, section, h):
        p0 = section['Static Pressure (Pa)']
        T0 = section['Standard Temperature (K)']
        L = section['Temperature Lapse Rate (K/m)']
        h0 = section['Geopotential Height (m)']
        k = section['k']
        R = __class__.R
        g = __class__.g0
        if L==0.0:
            p = k*exp(-g*h/R/T0)
        else:
            p = k*(R*T0 + R*L*(h-h0))**(-g/R/L)
        return p

    def T_from_section(self, section, h):
        T0 = section['Standard Temperature (K)']
        L = section['Temperature Lapse Rate (K/m)']
        h0 = section['Geopotential Height (m)']
        T = L*(h-h0) + T0
        return T

    
if __name__=="__main__":
    atmo = US76Atmo()
    atmo.test()
    rho = []
    zkm = linspace(0,80,10000)
    for zkmi in zkm:
        z = zkmi*1000.0
        rho.append(atmo.rho(z))

    rhoswitch = []
    zswitch = atmo.altitude_from_geopotential(atmo.geopotentials[:-1])
    for z in zswitch:
        rhoswitch.append(atmo.rho(z))
    plot(zkm, rho, 'b-')
    plot(zswitch/1000, rhoswitch, 'ko')
    show()

        

