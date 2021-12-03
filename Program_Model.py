import sys
import math
from scipy.interpolate import interp1d

class Model():
    def __init__(self, __mdot):
        if (__mdot > 1000):
            self.__mdot = __mdot / 3600.0
        else:
            self.__mdot = __mdot 

        self.__thi = 433   # Temp hot in, K
        self.__tci = 289   # Temp cold in, K
        self.__nt = 11     # Number of tubes
        self.__np = 4      # Number of passes
        self.__di = 0.0229  # Diameter inside, m
        self.__do = 0.0254  # Diameter outside, m
        self.__l = 17.0716 # Tube Length, m
        self.__rf = 0.002  # Fouling Factor, m^2 * K / W  
        self.__hh = 400.0      # HT coeff for hot side, W / m^2 
        self.area_outside = (math.pi * self.__do) * self.__np * self.__nt * self.__l

        self.__temps = [275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0,
                         315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 
                         355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 
                         395.0, 400.0, 410.0, 420.0, 430.0]

        self.__cps = [4211.0, 4198.0, 4189.0, 4184.0, 4181.0, 4179.0, 4178.0, 
                        4178.0, 4179.0, 4180.0, 4182.0, 4184.0, 4186.0, 4188.0, 
                        4191.0, 4195.0, 4199.0, 4203.0, 4209.0, 4214.0, 4217.0, 
                        4220.0, 4226.0, 4232.0, 4239.0, 4256.0, 4278.0, 4302.0, 
                        4331.0]
        self.__mus = [0.001652, 0.001422, 0.001225, 0.001080, 0.000959, 0.000855, 
                        0.000769, 0.000695, 0.000631, 0.000577, 0.000528, 0.000489, 
                        0.000453, 0.000420, 0.000389, 0.000365, 0.000343, 0.000324, 
                        0.000306, 0.000289, 0.000279, 0.000274, 0.000260, 0.000248, 
                        0.000237, 0.000217, 0.000200, 0.000185, 0.000173]

        self.__prs = [12.22, 10.26, 8.810, 7.560, 6.620, 5.830, 5.200, 4.620, 
                        4.160, 3.770, 3.420, 3.150, 2.880, 2.660, 2.450, 2.290, 
                        2.140, 2.020, 1.910, 1.800, 1.760, 1.700, 1.610, 1.530, 
                        1.470, 1.340, 1.240, 1.160, 1.090]

        self.__ks = [0.574, 0.582, 0.590, 0.598, 0.606, 0.613, 0.620, 0.628, 
                        0.634, 0.640, 0.645, 0.650, 0.656, 0.660, 0.668, 0.668, 
                        0.671, 0.674, 0.677, 0.679, 0.680, 0.681, 0.683, 0.685, 
                        0.686, 0.688, 0.688, 0.688, 0.685]

        self.__vfs = [1.000E-3, 1.000E-3, 1.000E-3, 1.001E-3, 1.002E-3, 1.003E-3, 1.005E-3, 
                        1.007E-3, 1.009E-3, 1.011E-3, 1.013E-3, 1.016E-3, 1.018E-3, 1.021E-3, 
                        1.024E-3, 1.027E-3, 1.030E-3, 1.034E-3, 1.038E-3, 1.041E-3, 1.044E-3, 
                        1.045E-3, 1.049E-3, 1.053E-3, 1.058E-3, 1.067E-3, 1.077E-3, 1.088E-3, 
                        1.099E-3]    

        self.__tcoReady = False
        self.__thoReady = False

        #Initial temperature guesses
        self.__tho_est = self.__thi - 65   # Temperature hot out, K
        self.__tco_est = self.__tci + 70  # Temperature cold out, K

        while True:
            self.__thavg = (self.__thi + self.__tho_est) / 2.0
            self.__tcavg = (self.__tci + self.__tco_est) / 2.0

            self.__interMu(self.__tcavg)
            self.__interCp(self.__tcavg)
            self.__interK(self.__tcavg)
            self.__interPr(self.__tcavg)
            self.__interVf(self.__tcavg)

            self.__funcRe_d(self.__mdot, self.__di, self.__mu, self.__nt)
            self.__funcNu(self.__Re_d, self.__pr)
            self.__funcHi(self.__Nu, self.__k, self.__di)
            self.__funcU(self.__hi, self.__do, self.__di, self.__rf, self.__hh)
            self.__funcCwater(self.__mdot, self.__cp)
            self.__funcCoil()
            self.__findCmin(self.__coil, self.__cw)
            self.__findCr(self.__coil, self.__cw)
            self.__funcNTU(self.__U, self.area_outside, self.__cmin)
            self.__funcEff(self.__NTU, self.__cr)
            self.__funcQmax(self.__cmin, self.__thi, self.__tci)
            self.__funcQ(self.__eff, self.__qmax)
            self.__funcThotOut(self.__thi, self.__q, self.__coil)
            self.__funcTcoldOut(self.__tci, self.__q, self.__cw)
            self.__funcCheckTout()
            

            if (self.__funcDiff(self.__tco_est,self.__tco) < 0.1):
                self.__tcoReady = True
            else:
                self.__tco_est += 5.0  
            if (self.__funcDiff(self.__tho_est,self.__tho) < 0.1):
                self.__thoReady = True
            else:
                self.__tho_est -= 5.0
            
            if ((self.__tcoReady and self.__thoReady) == True):
                break   
        
        self.__findF()
        self.__funcUm()
        self.__funcDeltaP()
        self.__funcThermalR()

    # Calculation Functions

    def __funcRe_d(self, m_dot, D_i, __mu, Nt):
        self.__Re_d = (4.0 * m_dot) / (math.pi * D_i * __mu * Nt)
        return self.__Re_d

    def __funcNu(self, __Re_d, Pr):
        if (__Re_d < 2300):
            self.__Nu = 4.36
        else:
            self.__Nu = 0.023 * (__Re_d ** 0.8) * (Pr ** 0.4)
        return self.__Nu

    def __funcHi(self, __Nu, __k, Di):
        self.__hi = (__Nu * __k) / Di
        return self.__hi

    def __funcU(self, __hi, Do, Di, Rf, ho):
        self.__U = ((Do / (__hi * Di)) + ((Do * Rf) / Di) + 1.0/ho) ** -1
        return self.__U

    def __funcCwater(self, m_dot, __cp):
        self.__cw = m_dot * __cp
        return self.__cw

    def __funcCoil(self):
        self.__coil = self.__cw * ((self.__tco_est - self.__tci) / (self.__thi - self.__tho_est))
        return self.__coil

    def __funcNTU(self, __U, A, __cmin):
        self.__NTU = (__U * A) / __cmin
        return self.__NTU

    def __funcEff(self, ntu, cr):  
        z = (1.0 + cr**2.0)**(0.5)
        y = math.exp(-ntu*z)
        u = (1.0 + y)/(1.0 - y)
        t = u * z
        w = (1.0 + cr + t)
        self.__eff = 2.0 * (w ** (-1))
        return self.__eff

    def __funcQmax(self,__cmin,__thi,__tci):
        self.__qmax = __cmin*(__thi-__tci)
        return self.__qmax

    def __funcQ(self,__eff,__qmax):
        self.__q = __eff * __qmax
        return self.__q

    def __funcThotOut(self, __thi, __q, ch):
        self.__tho = __thi - (__q/ch)
        return self.__tho

    def __funcTcoldOut(self, __tci, __q, cc):
        self.__tco = (__q/cc) + __tci
        return self.__tco
    
    def __funcUm(self):
        n = (self.__mdot * self.__vf * 4.0)
        d = (math.pi * (self.__di ** 2.0) * self.__nt)
        self.__um =  n / d
        return self.__um

    def __funcCheckTout(self):
        if (self.__tco > self.__tho):
            output = (self.__tci + self.__thi) / 2.0
            self.__tco = output
            self.__tho = output
            return True
        return False

    def __funcDeltaP(self):
        self.__dp = (self.__f * self.__l * self.__np * self.__um **2) / (2.0 * self.__di * self.__vf)
        return self.__dp

    def __funcThermalR(self):
        self.__rtherm = self.__q / (self.__tco - self.__tci)
        return self.__rtherm

    def __funcDiff(self,guess,result):
        return abs((result-guess) / guess)

    # Finding functions
    def __findCmin(self,cHot,cCold):
        self.__cmin = cHot if (cHot < cCold) else cCold
        self.cmax = cHot if (cHot > cCold) else cCold
        return self.__cmin

    def __findCr(self,cHot,cCold):
        if (cHot > cCold):
            self.__cr = cCold/cHot
        elif (cHot < cCold):
            self.__cr = cHot/cCold
        else:
            self.__cr = 1   
        return self.__cr

    def __findF(self):
        if (self.__Re_d < 2300):
            self.__f = 64 / self.__Re_d
        else:
            self.__f = (0.79 * math.log(self.__Re_d) - 1.64) ** -2.0
        return self.__f

    # Interpolation Functions
    def __interCp(self,T):
        cpInter = interp1d(self.__temps,self.__cps)
        self.__cp = cpInter(T)
        return self.__cp

    def __interMu(self,T):
        muInter = interp1d(self.__temps,self.__mus)
        self.__mu = muInter(T)
        return self.__mu

    def __interK(self,T):
        kInter = interp1d(self.__temps,self.__ks)
        self.__k = kInter(T)
        return self.__k

    def __interPr(self,T):
        prInter = interp1d(self.__temps,self.__prs)
        self.__pr = prInter(T)
        return self.__pr

    def __interVf(self,T):
        vfInter = interp1d(self.__temps,self.__vfs)
        self.__vf = vfInter(T)
        return self.__vf

    # Variable Access Functions
    def giveTcout(self):
        return round(self.__tco,6)
    def giveThout(self):
        return round(self.__tho,6)
    def giveCc(self):
        return round(self.__cw,6)
    def giveCh(self):
        return round(self.__coil,6)
    def giveRe_d(self):
        return round(self.__Re_d,6)
    def giveDp(self):
        return round(self.__dp,6)
    def giveNu(self):
        return round(self.__Nu,6)
    def giveH(self):
        return round(self.__hi,6)
    def giveThermR(self):
        return round(self.__rtherm,6)
    def giveNTU(self):
        return round(self.__NTU,6)
    def giveEff(self):
        return round(self.__eff,6)
    def giveQ(self):
        return round(self.__q,6)
    def giveMdot(self):
        return round(self.__mdot * 3600,2)
        
    def giveAll(self):
        array = [self.__mdot*3600, self.__tco, self.__tho, self.__cw, self.__coil, self.__Re_d, self.__dp, self.__Nu, self.__hi, self.__rtherm, self.__U, self.__NTU, self.__eff, self.__q]
        return array