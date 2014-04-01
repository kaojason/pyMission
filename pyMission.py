import numpy
import mission
import system

class flightAnalysis(object):
    def __init__(self, x, h, M, kw):
        
        # Dictionary setup:
        # kw = {
        #     'alpha': alpha,     (angle of attack)       [sizeN] {rad}
        #     'tau': tau,         (throttle setting)      [sizeN] {%}
        #     'eta': eta,         (tail rotation angle)   [sizeN] {rad}
        #     'Wf': Wf,           (weight of fuel)        [sizeN] {Newtons}
        #     'S': S,             (surface area of wing)  [scal]  {m^2}
        #     'Wac': Wac,         (weight without fuel)   [scal]  {Newtons}
        #     'cThrust': cThrust, (total thrust of eng)   [scal]  {Newtons}
        #     'SFC': SFC,         (TSFC at sea level)     [scal]  {kg/(s N)}
        #     'chord': chord,     (M.A.C)                 [scal]  {m}
        #     'inertia': inertia, (pitching MoI)          [scal]  {kg m^2}
        #     'AR': AR,           (AR of wing)            [scal]  {m}
        #     'e': e,             (Oswald's efficiency)   [scal]
        #     'CL': CL,           (lift coeff)            [sizeN]
        #     'CD': CD,           (drag coeff)            [sizeN]
        #     'CM': CM,           (pitching mmt coeff)    [sizeN]
        #     'Thrust' : Thrust,  (Thrust)                [sizeN] {Newtons}
        #     'numInt': numInt    (# integration/segment) [scal]
        # }

        self.x = x
        self.h = h
        self.M = M
        self.alpha = kw['alpha']
        self.tau = kw['tau']
        self.eta = kw['eta']
        self.Wf = kw['Wf']
        self.S = kw['S']
        self.Wac = kw['Wac']
        self.cThrustSL = kw['cThrust']
        self.SFCSL = kw['SFC']
        self.chord = kw['chord']
        self.AR = kw['AR']
        self.e = kw['e']
        self.inertia = kw['inertia']
        self.CL = kw['CL']
        self.CD = kw['CD']
        self.CM = kw['CM']
        self.Thrust = kw['Thrust']
        self.numInt = kw['numInt']
        self.numSeg = len(numpy.atleast_1d(x))

        self.g = 9.81

        self.dv = numpy.zeros(self.numSeg)
        self.gamma = numpy.zeros(self.numSeg)
        self.dGamma = numpy.zeros(self.numSeg)
        self.d2Gamma = numpy.zeros(self.numSeg)
        self.ddhdh = numpy.zeros(self.numSeg)
        self.ddvdv = numpy.zeros(self.numSeg)
        self.dgammadh = numpy.zeros(self.numSeg)
        self.dd2hdh = numpy.zeros(self.numSeg)
        self.ddgammadh = numpy.zeros(self.numSeg)
        self.dd3hdh = numpy.zeros(self.numSeg)
        self.dd2gammadh = numpy.zeros(self.numSeg)

    def getDerivatives(self):
        # compute the derivatives of v, h to get dvdx, gamma,
        # dGammadx, d2Gammadx2

        [self.dv,self.gamma,self.dGamma,
         self.d2Gamma] = mission.getderivatives(self.numSeg,
                                                self.x,self.h,self.v)

    def getDDerivatives(self):
        # compute the derivatives of the derivatives

        [self.ddhdh,self.ddvdv,self.dgammadh,self.dd2hdh,
         self.ddgammadh,self.dd3hdh,
         self.dd2gammadh] = mission.getdderivatives(self.numSeg,self.x,self.h,
                                                    self.v)

    def getSFC(self):
        # compute the SFC (TSFC) using a simple altitude relation

        SFC = (self.SFCSL + (6.39e-13)*self.h)
        return SFC

    def getDSFC(self):
        # compute the derivatives of TSFC

        self.dSFCdSFCSL = numpy.ones(self.numSeg)
        self.dSFCdH = numpy.ones(self.numSeg)*(6.39e-13)
        return

    def getV(self):
        # compute the flight velocity

        V = numpy.sqrt(1.4*287*self.T)*self.M
        return V
    
    def getDV(self):
        # compute the derivatives of flight velocity

        self.dVdTemp = self.Mach*numpy.power(1.4*287*self.T,-0.5)*1.4*287
        self.dVdM = numpy.sqrt(1.4*287*self.T)
        return

    def getTemp(self):
        # compute the atmospheric temperature

        Temp = (288.16 - 6.5e-3 * self.h)
        return Temp

    def getDTemp(self):
        # compute the derivatives of the atmosphreic temperature

        self.dTempdH = numpy.ones(self.numSeg)*(-6.5e-3)
        return

    def getRho(self):
        # compute the atmospheric density

        rho = (1.225*numpy.power(self.T/288.16,-(self.g/(-6.5e-3*287)+1)))
        return rho

    def getDRho(self):
        # compute the derivatives of the atmospheric density

        self.dRhodTemp = 1.225*numpy.power(self.T/288.16,-(self.g/(-6.5e-3*287)+2))*(1/288.16)
        return

    def getAlpha(self):
        # compute alpha from a given CL

        alpha = numpy.zeros(self.numSeg)
        alpha = mission.getalpha(self.numSeg, self.eta, self.CL)
        return alpha

    def getCD(self):
        # compute the drag coefficient

        CD = numpy.zeros(self.numSeg)
        CD = mission.getcd(self.numSeg, self.AR, self.e, self.CL)
        return CD

    def getEta(self):
        # compute eta from a given pitching moment coefficient

        eta = numpy.zeros(self.numSeg)
        eta = mission.geteta(self.numSeg, self.alpha, self.eta, self.CM)
        return eta

    def getTau(self):
        # compute the throttle setting from a given Thrust

        tau = numpy.zeros(self.numSeg)
        tau = mission.gettau(self.numSeg, self.cThrustSL, 
                             self.x, self.h, self.Thrust)
        return tau

    def getWf(self):
        # compute the fuel weight

        Wf = numpy.zeros(self.numSeg)
        R1 = numpy.linspace(0.0,1.0,num=self.numInt)
        R2 = numpy.linspace(1.0,0.0,num=self.numInt)

        Wf = mission.getwf(self.numInt,self.numSeg,self.x,self.v,
                           self.gamma,self.Thrust,self.SFC,
                           self.g,R1,R2)
        return Wf

    def getCLRes(self):
        # compute the residuals of the lift coefficient

        CLRes = numpy.zeros(self.numSeg)
        CLRes = mission.getclres(self.numSeg, self.numInt, self.Wac,
                                 self.S, self.g, self.x, self.v,
                                 self.rho, self.CL, self.Wf,
                                 self.gamma, self.Thrust, self.alpha,
                                 self.dGamma)
        return CLRes

    def getThrustRes(self):
        # compute the residuals of thrust

        ThrustRes = numpy.zeros(self.numSeg)
        ThrustRes = mission.getthrustres(self.numSeg, self.numInt, self.S,
                                         self.Wac, self.g, self.x, self.v,
                                         self.rho, self.gamma, self.dv, self.CD,
                                         self.Wf, self.alpha, self.Thrust)
        return ThrustRes

    def getCMRes(self):
        # compute the residuals of the pitching moment coefficient

        CMRes = numpy.zeros(self.numSeg)
        CMRes = mission.getcmres(self.numSeg, self.numInt, self.S,
                                 self.chord, self.inertia, self.x,
                                 self.v, self.rho, self.gamma,
                                 self.dGamma, self.d2Gamma, self.dv,
                                 self.CM)
        return CMRes

    def getDAlpha(self):
        # compute the derivatives of alpha

        self.dAlphadEta = numpy.zeros(self.numSeg)
        self.dAlphadCL = numpy.zeros(self.numSeg)
        
        [self.dAlphadEta, 
         self.dAlphadCL] = mission.getdalpha(self.numSeg, self.alpha, self.eta,
                                             self.CL)
        
        return

    def getDCD(self):
        # compute the derivatives of the drag coefficient

        self.dCDdAR = numpy.zeros(self.numSeg)
        self.dCDdE = numpy.zeros(self.numSeg)
        self.dCDdCL = numpy.zeros(self.numSeg)

        [self.dCDdAR, self.dCDdE, self.dCDdCL] \
            = mission.getdcd(self.numSeg, self.AR, self.e, self.CL, 
                             self.CD)

        return

    def getDEta(self):
        # compute the derivatives of eta

        self.dEtadAlpha = numpy.zeros(self.numSeg)
        self.dEtadCM = numpy.zeros(self.numSeg)

        [self.dEtadAlpha,
         self.dEtadCM] = mission.getdeta(self.numSeg, self.alpha, self.eta, 
                                         self.CM)

        return

    def getDTau(self):
        # compute the derivatives of the throttle setting

        self.dTaudCThrustSL = numpy.zeros(self.numSeg)
        self.dTaudH = numpy.zeros(self.numSeg)
        self.dTaudThrust = numpy.zeros(self.numSeg)

        [self.dTaudCThrustSL, self.dTaudH,
         self.dTaudThrust] = mission.getdtau(self.numSeg, self.cThrustSL, 
                                             self.x, self.h, self.tau, 
                                             self.Thrust)

        return

    def getDWf(self):
        # compute the derivatives of the fuel weight

        self.dWfdSFC1 = numpy.zeros(self.numSeg)
        self.dWfdSFC2 = numpy.zeros(self.numSeg)
        self.dWfdThrust1 = numpy.zeros(self.numSeg)
        self.dWfdThrust2 = numpy.zeros(self.numSeg)
        self.dWfdV1 = numpy.zeros(self.numSeg)
        self.dWfdV2 = numpy.zeros(self.numSeg)
        self.dWfdGamma1 = numpy.zeros(self.numSeg)
        self.dWfdGamma2 = numpy.zeros(self.numSeg)

        R1 = numpy.linspace(0.0,1.0,num=self.numInt)
        R2 = numpy.linspace(1.0,0.0,num=self.numInt)

        [self.dWfdSFC1, self.dWfdSFC2, self.dWfdThrust1,
         self.dWfdThrust2, self.dWfdV1, self.dWfdV2, 
         self.dWfdGamma1, 
         self.dWfdGamma2] = mission.getdwf(self.numSeg, self.numInt,
                                           self.g, self.x, self.v, self.gamma,
                                           self.Thrust, self.SFC,
                                           R1, R2)

        return

    def getDCL(self):
        # compute the derivatives of the residuals of the lift coefficient

        self.dCLdWac = numpy.zeros(self.numSeg)
        self.dCLdS = numpy.zeros(self.numSeg)
        self.dCLdV1 = numpy.zeros(self.numSeg)
        self.dCLdV2 = numpy.zeros(self.numSeg)
        self.dCLdV3 = numpy.zeros(self.numSeg)
        self.dCLdRho1 = numpy.zeros(self.numSeg)
        self.dCLdRho2 = numpy.zeros(self.numSeg)
        self.dCLdRho3 = numpy.zeros(self.numSeg)
        self.dCLdCL1 = numpy.zeros(self.numSeg)
        self.dCLdCL2 = numpy.zeros(self.numSeg)
        self.dCLdCL3 = numpy.zeros(self.numSeg)
        self.dCLdWf1 = numpy.zeros(self.numSeg)
        self.dCLdWf2 = numpy.zeros(self.numSeg)
        self.dCLdWf3 = numpy.zeros(self.numSeg)
        self.dCLdGamma1 = numpy.zeros(self.numSeg)
        self.dCLdGamma2 = numpy.zeros(self.numSeg)
        self.dCLdGamma3 = numpy.zeros(self.numSeg)
        self.dCLdThrust1 = numpy.zeros(self.numSeg)
        self.dCLdThrust2 = numpy.zeros(self.numSeg)
        self.dCLdThrust3 = numpy.zeros(self.numSeg)
        self.dCLdAlpha1 = numpy.zeros(self.numSeg)
        self.dCLdAlpha2 = numpy.zeros(self.numSeg)
        self.dCLdAlpha3 = numpy.zeros(self.numSeg)
        self.dCLddGamma1 = numpy.zeros(self.numSeg)
        self.dCLddGamma2 = numpy.zeros(self.numSeg)
        self.dCLddGamma3 = numpy.zeros(self.numSeg)

        [self.dCLdWac, self.dCLdS, self.dCLdV1, self.dCLdV2, self.dCLdV3,
         self.dCLdRho1, self.dCLdRho2, self.dCLdRho3, self.dCLdCL1,
         self.dCLdCL2, self.dCLdCL3, self.dCLdWf1, self.dCLdWf2, self.dCLdWf3,
         self.dCLdGamma1, self.dCLdGamma2, self.dCLdGamma3, self.dCLdThrust1,
         self.dCLdThrust2, self.dCLdThrust3, self.dCLdAlpha1, self.dCLdAlpha2,
         self.dCLdAlpha3, self.dCLddGamma1, self.dCLddGamma2,
         self.dCLddGamma3] = mission.getdcl(self.numSeg, self.numInt,
                                            self.Wac, self.S, self.g,
                                            self.x, self.v, self.rho,
                                            self.CL, self.Wf, self.gamma,
                                            self.Thrust, self.alpha,
                                            self.dGamma)
        return

    def getDThrust(self):
        # compute the derivatives of the residuals of thrust

        self.dTdS = numpy.zeros(self.numSeg)
        self.dTdWac = numpy.zeros(self.numSeg)
        self.dTdV1 = numpy.zeros(self.numSeg)
        self.dTdV2 = numpy.zeros(self.numSeg)
        self.dTdV3 = numpy.zeros(self.numSeg)
        self.dTdRho1 = numpy.zeros(self.numSeg)
        self.dTdRho2 = numpy.zeros(self.numSeg)
        self.dTdRho3 = numpy.zeros(self.numSeg)
        self.dTdGamma1 = numpy.zeros(self.numSeg)
        self.dTdGamma2 = numpy.zeros(self.numSeg)
        self.dTdGamma3 = numpy.zeros(self.numSeg)
        self.dTddV1 = numpy.zeros(self.numSeg)
        self.dTddV2 = numpy.zeros(self.numSeg)
        self.dTddV3 = numpy.zeros(self.numSeg)
        self.dTdCD1 = numpy.zeros(self.numSeg)
        self.dTdCD2 = numpy.zeros(self.numSeg)
        self.dTdCD3 = numpy.zeros(self.numSeg)
        self.dTdWf1 = numpy.zeros(self.numSeg)
        self.dTdWf2 = numpy.zeros(self.numSeg)
        self.dTdWf3 = numpy.zeros(self.numSeg)
        self.dTdAlpha1 = numpy.zeros(self.numSeg)
        self.dTdAlpha2 = numpy.zeros(self.numSeg)
        self.dTdAlpha3 = numpy.zeros(self.numSeg)
        self.dTdThrust1 = numpy.zeros(self.numSeg)
        self.dTdThrust2 = numpy.zeros(self.numSeg)
        self.dTdThrust3 = numpy.zeros(self.numSeg)

        [self.dTdS, self.dTdWac, self.dTdV1, self.dTdV2, self.dTdV3,
         self.dTdRho1, self.dTdRho2, self.dTdRho3, self.dTdGamma1,
         self.dTdGamma2, self.dTdGamma3, self.dTddV1, self.dTddV2,
         self.dTddV3, self.dTdCD1, self.dTdCD2, self.dTdCD3, self.dTdWf1,
         self.dTdWf2, self.dTdWf3, self.dTdAlpha1, self.dTdAlpha2,
         self.dTdAlpha3, self.dTdThrust1, self.dTdThrust2, 
         self.dTdThrust3] = mission.getdthrust(self.numSeg, self.numInt, self.S,
                                               self.Wac, self.g, self.x, self.v,
                                               self.rho, self.gamma, self.dv,
                                               self.CD, self.Wf, self.alpha,
                                               self.Thrust)
        return

    def getDCM(self):
        # compute the derivatives of the residuals of the pitching moment coefficient

        self.dCMdS = numpy.zeros(self.numSeg)
        self.dCMdC = numpy.zeros(self.numSeg)
        self.dCMdI = numpy.zeros(self.numSeg)
        self.dCMdV1 = numpy.zeros(self.numSeg)
        self.dCMdV2 = numpy.zeros(self.numSeg)
        self.dCMdV3 = numpy.zeros(self.numSeg)
        self.dCMdRho1 = numpy.zeros(self.numSeg)
        self.dCMdRho2 = numpy.zeros(self.numSeg)
        self.dCMdRho3 = numpy.zeros(self.numSeg)
        self.dCMdGamma1 = numpy.zeros(self.numSeg)
        self.dCMdGamma2 = numpy.zeros(self.numSeg)
        self.dCMdGamma3 = numpy.zeros(self.numSeg)
        self.dCMddGamma1 = numpy.zeros(self.numSeg)
        self.dCMddGamma2 = numpy.zeros(self.numSeg)
        self.dCMddGamma3 = numpy.zeros(self.numSeg)
        self.dCMdd2Gamma1 = numpy.zeros(self.numSeg)
        self.dCMdd2Gamma2 = numpy.zeros(self.numSeg)
        self.dCMdd2Gamma3 = numpy.zeros(self.numSeg)
        self.dCMddV1 = numpy.zeros(self.numSeg)
        self.dCMddV2 = numpy.zeros(self.numSeg)
        self.dCMddV3 = numpy.zeros(self.numSeg)
        self.dCMdCM1 = numpy.zeros(self.numSeg)
        self.dCMdCM2 = numpy.zeros(self.numSeg)
        self.dCMdCM3 = numpy.zeros(self.numSeg)

        [self.dCMdS, self.dCMdC, self.dCMdI, self.dCMdV1,
         self.dCMdV2, self.dCMdV3, self.dCMdRho1, self.dCMdRho2, self.dCMdRho3,
         self.dCMdGamma1, self.dCMdGamma2, self.dCMdGamma3, self.dCMddGamma1,
         self.dCMddGamma2, self.dCMddGamma3, self.dCMdd2Gamma1,
         self.dCMdd2Gamma2, self.dCMdd2Gamma3, self.dCMddV1, self.dCMddV2,
         self.dCMddV3, self.dCMdCM1, self.dCMdCM2,
         self.dCMdCM3] = mission.getdcm(self.numSeg, self.numInt, self.S,
                                        self.chord, self.inertia, self.x,
                                        self.v, self.rho, self.gamma,
                                        self.dGamma, self.d2Gamma,
                                        self.dv, self.CM)
        return

class Var_SFC(ExplicitSystem):
    
    def _declare_global(self):
        return 'SFC', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*8.951e-6,
                                     lb=numpy.zeros(len(segs)))
        self._declare_local_argument('SFCSL')
        self._declare_local_argument('h', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        self.uVec()[:] = fa.getSFC()
        fa.SFC = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDSFC()

        if mode == 'fwd':
            if self._ID('SFCSL') in arguments:
                self.dgVec()[:] += fa.dSFCdSFCSL * self.dpVec('SFCSL')[:]
            if self._ID('h') in arguments:
                self.dgVec()[:] += fa.dSFCdH * self.dpVec('h')[:]
        if mode == 'rev':
            if self._ID('SFCSL') in arguments:
                self.dpVec('SFCSL')[:] = fa.dSFCdSFCSL * self.dgVec()[:]
            if self._ID('h') in arguments:
                self.dpVec('h')[:] = fa.dSFCdH * self.dgVec()[:]

class Var_v(ExplicitSystem):

    def _declare_global(self):
        return 'v', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*300.0,
                                     lb=numpy.ones(len(segs)))
        self._declare_local_argument('T', indices=segs)
        self._declare_local_argument('M', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        self.uVec()[:] = fa.getV()
        fa.v = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDV()

        if mode == 'fwd':
            if self._ID('T') in arguments:
                self.dgVec()[:] += fa.dVdTemp * self.dpVec('T')[:]
            if self._ID('M') in arguments:
                self.dgVec()[:] += fa.dVdM * self.dpVec('M')[:]
        if mode == 'rev':
            if self._ID('T') in arguments:
                self.dpVec('T')[:] += fa.dVdTemp * self.dgVec()[:]
            if self._ID('M') in arguments:
                self.dpVec('M')[:] += fa.dVdM * self.dgVec()[:]

class Var_temp(ExplicitSystem):

    def _declare_global(self):
        return 'temp', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*200.0,
                                     lb=numpy.ones(len(segs)))
        self._declare_local_argument('h', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        self.uVec()[:] = fa.getTemp()
        fa.T = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDTemp()

        if mode == 'fwd':
            if self._ID('h') in arguments:
                self.dgVec()[:] += fa.dTempdH * self.dpVec('h')[:]
        if mode == 'rev':
            if self._ID('h') in arguments:
                self.dpVec('h')[:] += fa.dTempdH * self.dgVec()[:]

class Var_rho(ExplicitSystem):
    
    def _declare_global(self):
        return 'rho', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs)),
                                     lb=numpy.ones(len(segs))*0.01)
        self._declare_local_argument('T', indices=segs)
        self._declare_local_argument('g')

    def _apply_G(self):
        fa = self.kwargs['fa']
        self.uVec()[:] = fa.getRho()
        fa.rho = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDRho()

        if mode == 'fwd':
            if self._ID('T') in arguments:
                self.dgVec()[:] += fa.dRhodTemp * self.dpVec('T')[:]
        if mode == 'rev':
            if self._ID('T') in arguments:
                self.dfVec('T')[:] += fa.dRhodTemp * self.dgVec()[:]

class Var_gamma(ExplicitSystem):

    def _declare_global(self):
        return 'gamma', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.zeros(len(segs)))
        self._declare_local_argument('x', indices=segs)
        self._declare_local_argument('h', indices=segs)
        self._declare_local_argument('v', indices=segs)

    def _apply_G(self):

class Var_alpha(ExplicitSystem):

    def _declare_global(self):
        return 'alpha', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*pi/180.0)
        self._declare_local_argument('eta', indices=segs)
        self._declare_local_argument('CL', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        fa.eta = self.pVec('eta')[:]
        fa.CL = self.pVec('CL')[:]
        self.vVec()[:] = fa.getAlpha()
        fa.alpha = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDAlpha()

        if mode == 'fwd':
            if self._ID('eta') in arguments:
                self.dgVec()[:] += fa.dAlphadEta * self.dpVec('eta')[:]
            if self._ID('CL') in arguments:
                self.dgVec()[:] += fa.dAlphadCL * self.dpVec('CL')[:]
        if mode == 'rev':
            if self._ID('eta') in arguments:
                self.dpVec('eta')[:] += fa.dAlphadEta * self.dgVec()[:]
            if self._ID('CL') in arguments:
                self.dpVec('CL')[:] += fa.dAlphadCL * self.dgVec()[:]

class Var_CD(ExplicitSystem):
    
    def _declare_global(self):
        return 'CD', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*1e-5,
                                     lb=numpy.zeros(len(segs)))
        self._declare_local_argument('AR')
        self._declare_local_argument('e')
        self._declare_local_argument('CL', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        fa.CL = self.pVec('CL')[:]
        self.vVec()[:] = fa.getCD()
        fa.CD = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDCD()

        if mode == 'fwd':
            if self._ID('AR') in arguments:
                self.dgVec()[:] += fa.dCDdAR * self.dpVec('AR')[:]
            if self._ID('e') in arguments:
                self.dgVec()[:] += fa.dCDdE * self.dpVec('e')[:]
            if self._ID('CL') in arguments:
                self.dgVec()[:] += fa.dCDdCL * self.dpVec('CL')[:]
        if mode == 'rev':
            if self._ID('AR') in arguments:
                self.dpVec('AR')[:] += fa.dCDdAR * self.dgVec()[:]
            if self._ID('e') in arguments:
                self.dpVec('e')[:] += fa.dCDdE * self.dgVec()[:]
            if self._ID('CL') in arguments:
                self.dpVec('CL')[:] += fa.dCDdCL * self.dgVec()[:]

class Var_eta(ExplicitSystem):
    
    def _declare_global(self):
        return 'eta', 1
    
    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*(-pi/180.0))
        self._declare_local_argument('alpha', indices=segs)
        self._declare_local_argument('CM', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        fa.alpha = self.pVec('alpha')[:]
        fa.CM = self.pVec('CM')[:]
        self.uVec()[:] = fa.getEta()
        fa.eta = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDEta()

        if mode == 'fwd':
            if self._ID('alpha') in arguments:
                self.dgVec()[:] += fa.dEtadAlpha * self.dpVec('alpha')[:]
            if self._ID('CM') in arguments:
                self.dgVec()[:] += fa.dEtadCM * self.dpVec('CM')[:]
        if mode == 'rev':
            if self._ID('alpha') in arguments:
                self.dpVec('alpha')[:] += fa.dEtadAlpha * self.dgVec()[:]
            if self._ID('CM') in arguments:
                self.dpVec('CM')[:] += fa.dEtadCM * self.dgVec()[:]

class Var_tau(ExplicitSystem):
    
    def _declare_global(self):
        return 'tau', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*0.5,
                                     lb=numpy.zeros(len(segs)),ub=numpy.ones(len(segs)))
        self._declare_local_argument('cThrustSL')
        self._declare_local_argument('x', indices=segs)
        self._declare_local_argument('h', indices=segs)
        self._declare_local_argument('Thrust', indices=segs)

    def _apply_G(self):
        fa = self.kwargs['fa']
        fa.Thrust = self.pVec('Thrust')[:]
        self.uVec()[:] = fa.getTau()
        fa.tau = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDTau()

        if mode == 'fwd':
            if self._ID('cThrustSL') in arguments:
                self.dgVec()[:] += fa.dTaudCThrustSL * self.dpVec('cThrustSL')[:]
            if self._ID('h') in arguments:
                self.dgVec()[:] += fa.dTaudH * self.dpVec('h')[:]
            if self._ID('Thrust') in arguments:
                self.dgVec()[:] += fa.dTaudThrust * self.dpVec('Thrust')[:]
        if mode == 'rev':
            if self._ID('cThrustSL') in arguments:
                self.dpVec('cThrustSL')[:] += fa.dTaudCThrustSL * self.dgVec()[:]
            if self._ID('h') in arguments:
                self.dpVec('h')[:] += fa.dTaudH * self.dgVec()[:]
            if self._ID('Thrust') in arguments:
                self.dpVec('Thrust')[:] += fa.dTaudThrust * self.dgVec()[:]

class Var_Wf(ExplicitSystem):
    
    def _declare_global(self):
        return 'Wf', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.linspace(20000.0,0.0,len(segs)),
                                     lb=numpy.zeros(len(segs)))
        self._declare_local_argument('x', indices=segs)
        self._declare_local_argument('v', indices=segs)
        self._declare_local_argument('gamma', indices=segs)
        self._declare_local_argument('Thrust', indices=segs)
        self._declare_local_argument('SFC', indices=segs)
        self._declare_local_argument('g')

    def _apply_G(self):
        fa = self.kwargs['fa']
        fa.Thrust = self.pVec('Thrust')
        self.uVec()[:] = fa.getWf()
        fa.Wf = self.uVec()[:]

    def _apply_dGdp(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDWf()

        if mode == 'fwd':
            if self._ID('v') in arguments:
                self.dgVec()[:] += (fa.dWfdV1+fa.dWfdV2) * self.dpVec('v')[:]
            if self._ID('gamma') in arguments:
                self.dgVec()[:] += (fa.dWfdGamma1+fa.dWfdGamma2) * self.dpVec('gamma')[:]
            if self._ID('Thrust') in arguments:
                self.dgVec()[:] += (fa.dWfdThrust1+fa.dWfdThrust2) * self.dpVec('Thrust')[:]
            if self._ID('SFC') in arguments:
                self.dgVec()[:] += (fa.dWfdSFC1+fa.dWfdSFC2) * self.dpVec('SFC')[:]
        if mode == 'rev':
            if self._ID('v') in arguments:
                self.dpVec('v')[:] += (fa.dWfdV1+fa.dWfdV2) * self.dgVec()[:]
            if self._ID('gamma') in arguments:
                self.dpVec('gamma')[:] += (fa.dWfdGamma1+fa.dWfdGamma2) * self.dgVec()[:]
            if self._ID('Thrust') in arguments:
                self.dpVec('Thrust')[:] += (fa.dWfdThrust1+fa.dWfdThrust2) * self.dgVec()[:]
            if self._ID('SFC') in arguments:
                self.dpVec('SFC')[:] += (fa.dWfdSFC1+fa.dWfdSFC2) * self.dgVec()[:]

class Var_CL(ImplicitSystem):

    def _declare_global(self):
        return 'CL', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs)))
        self._declare_local_argument('Wac')
        self._declare_local_argument('S')
        self._declare_local_argument('g')
        self._declare_local_argument('x', indices=segs)
        self._declare_local_argument('v', indices=segs)
        self._declare_local_argument('rho', indices=segs)
        self._declare_local_argument('Wf', indices=segs)
        self._declare_local_argument('gamma', indices=segs)
        self._declare_local_argument('Thrust', indices=segs)
        self._declare_local_argument('alpha', indices=segs)
        self._declare_local_argument('dGamma', indices=segs)

    def _apply_F(self):
        fa = self.kwargs['fa']
        fa.CL = self.uVec()[:]
        fa.Wf = self.pVec('Wf')[:]
        fa.Thrust = self.pVec('Thrust')[:]
        fa.alpha = self.pVec('alpha')[:]
        self.fVec()[:] = fa.getCLRes()

    def _apply_dFdpu(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDCL()

        if mode == 'fwd':
            if self._ID('Wac') in arguments:
                self.dfVec()[:] += fa.dCLdWac * self.dpVec('Wac')[:]
            if self._ID('S') in arguments:
                self.dfVec()[:] += fa.dCLdS * self.dpVec('S')[:]
            if self._ID('v') in arguments:
                self.dfVec()[:] += (fa.dCLdV1+fa.dCLdV2+fa.dCLdV3) * self.dpVec('v')[:]
            if self._ID('rho') in arguments:
                self.dfVec()[:] += (fa.dCLdRho1+fa.dCLdRho2+fa.dCLdRho3) * self.dpVec('rho')[:]
            if self._ID('Wf') in arguments:
                self.dfVec()[:] += (fa.dCLdWf1+fa.dCLdWf2+fa.dCLdWf3) * self.dpVec('Wf')[:]
            if self._ID('gamma') in arguments:
                self.dfVec()[:] += (fa.dCLdGamma1+fa.dCLdGamma2+fa.dCLdGamma3) * self.dpVec('gamma')[:]
            if self._ID('Thrust') in arguments:
                self.dfVec()[:] += (fa.dCLdThrust1+fa.dCLdThrust2+fa.dCLdThrust3) * self.dpVec('Thrust')[:]
            if self._ID('alpha') in arguments:
                self.dfVec()[:] += (fa.dCLdAlpha1+fa.dCLdAlpha2+fa.dCLdAlpha3) * self.dpVec('alpha')[:]
            if self._ID('dGamma') in arguments:
                self.dfVec()[:] += (fa.dCLddGamma1+fa.dCLddGamma2+fa.dCLddGamma3) * self.dpVec('dGamma')[:]
            if self._ID('CL') in arguments:
                self.dfVec()[:] += (fa.dCLdCL1+fa.dCLdCL2+fa.dCLdCL3) * self.duVec('CL')[:]

        elif mode == 'rev':
            if self._ID('Wac') in arguments:
                self.dpVec('Wac')[:] += fa.dCLdWac * self.dfVec()[:]
            if self._ID('S') in arguments:
                self.dpVec('S')[:] += fa.dCLdS * self.dfVec()[:]
            if self._ID('v') in arguments:
                self.dpVec('v')[:] += (fa.dCLdV1+fa.dCLdV2+fa.dCLdV3) * self.dfVec()[:]
            if self._ID('rho') in arguments:
                self.dpVec('rho')[:] += (fa.dCLdRho1+fa.dCLdRho2+fa.dCLdRho3) * self.dfVec()[:]
            if self._ID('Wf') in arguments:
                self.dpVec('Wf')[:] += (fa.dCLdWf1+fa.dCLdWf2+fa.dCLdWf3) * self.dfVec()[:]
            if self._ID('gamma') in arguments:
                self.dpVec('gamma')[:] += (fa.dCLdGamma1+fa.dCLdGamma2+fa.dCLdGamma3) * self.dfVec()[:]
            if self._ID('Thrust') in arguments:
                self.dpVec('Thrust')[:] += (fa.dCLdThrust1+fa.dCLdThrust2+fa.dCLdThrust3) * self.dfVec()[:]
            if self._ID('alpha') in arguments:
                self.dpVec('alpha')[:] += (fa.dCLdAlpha1+fa.dCLdAlpha2+fa.dCLdAlpha3) * self.dfVec()[:]
            if self._ID('dGamma') in arguments:
                self.dpVec('dGamma')[:] += (fa.dCLddGamma1+fa.dCLddGamma2+fa.dCLddGamma3) * self.dfVec()[:]
            if self._ID('CL') in arguments:
                self.duVec('CL')[:] += (fa.dCLdCL1+fa.dCLdCL2+fa.dCLdCL3) * self.dfVec()[:]

class Var_Thrust(ImplicitSystem):
    
    def _declare_global(self):
        return 'Thrust', 1
    
    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.ones(len(segs))*400.0e3,
                                     lb=numpy.zeros(len(segs)))
        self._declare_local_argument('S')
        self._declare_local_argument('Wac')
        self._declare_local_argument('g')
        self._declare_local_argument('x', indices=segs)
        self._declare_local_argument('v', indices=segs)
        self._declare_local_argument('rho', indices=segs)
        self._declare_local_argument('gamma', indices=segs)
        self._declare_local_argument('dv', indices=segs)
        self._declare_local_argument('CD', indices=segs)
        self._declare_local_argument('Wf', indices=segs)
        self._declare_local_argument('alpha', indices=segs)

    def _apply_F(self):
        fa = self.kwargs['fa']
        fa.CD = self.pVec('CD')[:]
        fa.Wf = self.pVec('Wf')[:]
        fa.alpha = self.pVec('alpha')[:]
        fa.Thrust = self.uVec()[:]
        self.fVec()[:] = fa.getThrustRes()

    def _apply_dFdpu(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDThrust()

        if mode == 'fwd':
            if self._ID('S') in arguments:
                self.dfVec()[:] += fa.dTdS * self.dpVec('S')[:]
            if self._ID('Wac') in arguments:
                self.dfVec()[:] += fa.dTdWac * self.dpVec('Wac')[:]
            if self._ID('v') in arguments:
                self.dfVec()[:] += (fa.dTdV1+fa.dTdV2+fa.dTdV3) * self.dpVec('v')[:]
            if self._ID('rho') in arguments:
                self.dfVec()[:] += (fa.dTdRho1+fa.dTdRho2+fa.dTdRho3) * self.dpVec('rho')[:]
            if self._ID('gamma') in arguments:
                self.dfVec()[:] += (fa.dTdGamma1+fa.dTdGamma2+fa.dTdGamma3) * self.dpVec('gamma')[:]
            if self._ID('dv') in arguments:
                self.dfVec()[:] += (fa.dTddV1+fa.dTddV2+fa.dTddV3) * self.dpVec('dv')[:]
            if self._ID('CD') in arguments:
                self.dfVec()[:] += (fa.dTdCD1+fa.dTdCD2+fa.dTdCD3) * self.dpVec('CD')[:]
            if self._ID('WF') in arguments:
                self.dfVec()[:] += (fa.dTdWf1+fa.dTdWf2+fa.dTdWf3) * self.dpVec('Wf')[:]
            if self._ID('alpha') in arguments:
                self.dfVec()[:] += (fa.dTdAlpha1+fa.dTdAlpha2+fa.dTdAlpha3) * self.dpVec('alpha')[:]
            if self._ID('Thrust') in arguments:
                self.dfVec()[:] += (fa.dTdThrust1+fa.dTdThrust2+fa.dTdThrust3) * self.duVec('Thrust')[:]

        elif mode == 'rev':
            if self._ID('S') in arguments:
                self.dpVec('S')[:] += fa.dTdS * self.dfVec()[:]
            if self._ID('Wac') in arguments:
                self.dpVec('Wac')[:] += fa.dTdWac * self.dfVec()[:]
            if self._ID('v') in arguments:
                self.dpVec('v')[:] += (fa.dTdV1+fa.dTdV2+fa.dTdV3) * self.dfVec()[:]
            if self._ID('rho') in arguments:
                self.dpVec('rho')[:] += (fa.dTdRho1+fa.dTdRho2+fa.dTdRho3) * self.dfVec()[:]
            if self._ID('gamma') in arguments:
                self.dpVec('gamma')[:] += (fa.dTdGamma1+fa.dTdGamma2+fa.dTdGamma3) * self.dfVec()[:]
            if self._ID('dv') in arguments:
                self.dpVec('dv')[:] += (fa.dTddV1+fa.dTddV2+fa.dTddV3) * self.dfVec()[:]
            if self._ID('CD') in arguments:
                self.dpVec('CD')[:] += (fa.dTdCD1+fa.dTdCD2+fa.dTdCD3) * self.dfVec()[:]
            if self._ID('Wf') in arguments:
                self.dpVec('Wf')[:] += (fa.dTdWf1+fa.dTdWf2+fa.dTdWf3) * self.dfVec()[:]
            if self._ID('alpha') in arguments:
                self.dpVec('alpha')[:] += (fa.dTdAlpha1+fa.dTdAlpha2+fa.dTdAlpha3) * self.dfVec()[:]
            if self._ID('Thrust') in arguments:
                self.duVec('Thrust')[:] += (fa.dTdThrust1+fa.dTdThrust2+fa.dTdThrust3) * self.dfVec()[:]

class Var_CM(ImplicitSystem):
    
    def _declare_global(self):
        return 'CM', 1

    def _declare_local(self):
        segs = range(self.kwargs['numSeg'])
        self._declare_local_variable(len(segs),val=numpy.zeros(len(segs)))
        self._declare_local_argument('S')
        self._declare_local_argument('chord')
        self._declare_local_argument('inertia')
        self._declare_local_argument('x', indices=segs)
        self._declare_local_argument('v', indices=segs)
        self._declare_local_argument('rho', indices=segs)
        self._declare_local_argument('gamma', indices=segs)
        self._declare_local_argument('dGamma', indices=segs)
        self._declare_local_argument('d2Gamma', indices=segs)
        self._declare_local_argument('dv', indices=segs)

    def _apply_F(self):
        fa = self.kwargs['fa']
        fa.CM = self.vVec()[:]
        self.cVec()[:] = fa.getCMRes()

    def _apply_dFdpu(self, mode, arguments):
        fa = self.kwargs['fa']
        fa.getDCM()

        if mode == 'fwd':
            if self._ID('S') in arguments:
                self.dfVec()[:] += fa.dCMdS * self.dpVec('S')[:]
            if self._ID('chord') in arguments:
                self.dfVec()[:] += fa.dCMdC * self.dpVec('chord')[:]
            if self._ID('inertia') in arguments:
                self.dfVec()[:] += fa.dCMdI * self.dpVec('inertia')[:]
            if self._ID('v') in arguments:
                self.dfVec()[:] += (fa.dCMdV1+fa.dCMdV2+fa.dCMdV3) * self.dpVec('v')[:]
            if self._ID('rho') in arguments:
                self.dfVec()[:] += (fa.dCMdRho1+fa.dCMdRho2+fa.dCMdRho3) * self.dpVec('rho')[:]
            if self._ID('gamma') in arguments:
                self.dfVec()[:] += (fa.dCMdGamma1+fa.dCMdGamma2+fa.dCMdGamma3) * self.dpVec('gamma')[:]
            if self._ID('dGamma') in arguments:
                self.dfVec()[:] += (fa.dCMddGamma1+fa.dCMddGamma2+fa.dCMddGamma3) * self.dpVec('dGamma')[:]
            if self._ID('d2Gamma') in arguments:
                self.dfVec()[:] += (fa.dCMdd2Gamma1+fa.dCMdd2Gamma2+fa.dCMdd2Gamma3) * self.dpVec('d2Gamma')[:]
            if self._ID('dv') in arguments:
                self.dfVec()[:] += (fa.dCMddV1+fa.dCMddV2+fa.dCMddV3) * self.dpVec('dv')[:]
            if self._ID('CM') in arguments:
                self.dfVec()[:] += (fa.dCMdCM1+fa.dCMdCM2+fa.dCMdCM3) * self.duVec('CM')[:]

        elif mode == 'rev':
            if self._ID('S') in arguments:
                self.dpVec('S')[:] += fa.dCMdS * self.dfVec()[:]
            if self._ID('chord') in arguments:
                self.dpVec('chord')[:] += fa.dCMdC * self.dfVec()[:]
            if self._ID('inertia') in arguments:
                self.dpVec('inertia')[:] += fa.dCMdI * self.dfVec()[:]
            if self._ID('v') in arguments:
                self.dpVec('v')[:] += (fa.dCMdV1+fa.dCMdV2+fa.dCMdV3) * self.dfVec()[:]
            if self._ID('rho') in arguments:
                self.dpVec('rho')[:] += (fa.dCMdRho1+fa.dCMdRho2+fa.dCMdRho3) * self.dfVec()[:]
            if self._ID('gamma') in arguments:
                self.dpVec('gamma')[:] += (fa.dCMdGamma1+fa.dCMdGamma2+fa.dCMdGamma3) * self.dfVec()[:]
            if self._ID('dGamma') in arguments:
                self.dpVec('dGamma')[:] += (fa.dCMddGamma1+fa.dCMddGamma2+fa.dCMddGamma3) * self.dfVec()[:]
            if self._ID('d2Gamma') in arguments:
                self.dpVec('d2Gamma')[:] += (fa.dCMdd2Gamma1+fa.dCMdd2Gamma2+fa.dCMdd2Gamma3) * self.dfVec()[:]
            if self._ID('dv') in arguments:
                self.dpVec('dv')[:] += (fa.dCMddV1+fa.dCMddV2+fa.dCMddV3) * self.dfVec()[:]
            if self._ID('CM') in arguments:
                self.duVec('CM')[:] += (fa.dCMdCM1+fa.dCMdCM2+fa.dCMdCM3) * self.dfVec()[:]

x = numpy.linspace(0.0,1000.0,100)*1000.0
h = numpy.zeros(100)
h[0:4] = numpy.linspace(0.0,2000.0,4)
h[4:9] = numpy.ones(5)*2000.0
h[9:17] = numpy.linspace(2000.0,6500.0,8)
h[17:22] = numpy.ones(5)*6500.0
h[22:30] = numpy.linspace(6500.0,11000.0,8)
h[30:75] = numpy.ones(45)*11000.0
h[75:100] = numpy.linspace(11000.0,0.0,25)
M = numpy.zeros(100)
M[0:30] = numpy.linspace(0.2,0.84,30)
M[30:75] = numpy.ones(45)*0.84
M[75:100] = numpy.linspace(0.84,0.2,25)
S = 427.8
Wac = 185953*9.81
cThrust = 1020000.0
SFCSL = 8.9513-6
chord = 8.15
inertia = 4.1e7
AR = 8.68
e = 0.8
g = 9.81

params = {'alpha':numpy.ones(100)*3.0*numpy.pi/180.0,
          'tau':numpy.ones(100)*0.3,
          'eta':numpy.zeros(100),
          'Wf':numpy.linspace(100000,0.0,100)*9.81,
          'S':427.8,
          'Wac':185953*9.81,
          'cThrust':1020000,
          'SFC':8.951e-6,
          'chord':8.15,
          'inertia':4.1e7,
          'AR':8.68,
          'e':0.8,
          'CL':numpy.ones(100)*0.3,
          'CD':numpy.ones(100)*0.0015,
          'CM':numpy.zeros(100),
          'Thrust':numpy.ones(100)*0.3*1020000,
          'numInt':100
          }

flightProblem = flightAnalysis(x,h,M,params)


Var_x = IndependentSystem('x',0,value=x)
Var_h = IndependentSystem('h',0,value=h)
Var_M = IndependentSystem('M',0,value=M)
Var_S = IndependentSystem('S',0,value=S,size=1)
Var_Wac = IndependentSystem('Wac',0,value=Wac,size=1)
Var_cThrustSL = IndependentSystem('cThrustSL',0,value=cThrustSL,size=1)
Var_SFCSL = IndependentSystem('SFCSL',0,value=SFCSL,size=1)
Var_chord = IndependentSystem('chord',0,value=chord,size=1)
Var_inertia = IndependentSystem('inertia',0,value=inertia,size=1)
Var_AR = IndependentSystem('AR',0,value=AR,size=1)
Var_e = IndependentSystem('e',0,value=e,size=1)
Var_g = IndependentSystem('g',0,value=g,size=1)

main = SerialSystem('main',[
        Var_x,
        Var_h,
        Var_M,
        Var_S,
        Var_Wac,
        Var_cThrustSL,
        Var_SFCSL,
        Var_chord,
        Var_inertia,
        Var_AR,
        Var_e,
        Var_g,
        Var_temp(),
        Var_rho(),
        Var_v(),
        Var_SFC(),
        Var_gamma(),
        Var_dGamma(),
        Var_d2Gamma(),
        Var_dv(),
        SerialSystem('aeroFlight',[
                SerialSystem('aero',[
                        Var_alpha(),
                        Var_CD(),
                        Var_Eta(),
                        Var_Tau(),
                        ]),
                SerialSystem('flight',[
                        Var_CL(),
                        Var_Thrust(),
                        Var_CM(),
                        ])
                ])
        ])
