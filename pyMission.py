import numpy
import mission

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

        [self.dv,self.gamma,self.dGamma,
         self.d2Gamma] = mission.getderivatives(self.numSeg,
                                                self.x,self.h,self.v)

        [self.ddhdh,self.ddvdv,self.dgammadh,self.dd2hdh,
         self.ddgammadh,self.dd3hdh,
         self.dd2gammadh] = mission.getdderivatives(self.numSeg,self.x,self.h,
                                                    self.v)

    def getSFCRes(self):
        # compute the residuals of SFC (TSFC) using a simple altitude relation

        SFCRes = self.SFC - (self.SFCSL + (6.39e-13)*self.h)
        return SFCRes

    def getDSFC(self):
        # compute the derivatives of the residual of TSFC

        self.dSFCdSFC = numpy.ones(self.numSeg)
        self.dSFCdSFCSL = numpy.ones(self.numSeg)*-1
        self.dSFCdh = numpy.ones(self.numSeg)*(-6.39e-13)
        return

    def getVRes(self):
        # compute the residuals of the flight velocity

        VRes = self.v - numpy.sqrt(1.4*287*self.T)*self.M
        return VRes
    
    def getDV(self):
        # compute the derivatives of the residual of flight velocity

        self.dVdV = numpy.ones(self.numSeg)
        self.dVdTemp = self.Mach*numpy.power(1.4*287*self.T,-0.5)*1.4*287
        self.dVdM = numpy.sqrt(1.4*287*self.T)
        return

    def getTempRes(self):
        # compute the residuals of the atmospheric temperature

        TempRes = self.T - (288.16 - 6.5e-3 * self.h)
        return TempRes

    def getDTemp(self):
        # compute the derivatives of the residual of the atmosphreic temperature

        self.dTempdTemp = numpy.ones(self.numSeg)
        self.dTempdh = numpy.ones(self.numSeg)*6.5e-3
        return

    def getRhoRes(self):
        # compute the residuals of the atmospheric density

        RhoRes = self.rho - (1.225*numpy.power(self.T/288.16,-(self.g/(-6.5e-3*287)+1)))
        return RhoRes

    def getDRho(self):
        # compute the derivatives of the residual of the atmospheric density

        self.dRhodRho = numpy.ones(self.numSeg)
        self.dRhodTemp = -1.225*numpy.power(self.T/288.16,-(self.g/(-6.5e-3*287)+2))*(1/288.16)
        return

    def getAlphaRes(self):
        # compute the residuals of alpha from a given CL

        alphaRes = numpy.zeros(self.numSeg)
        alphaRes = mission.getalphares(self.numSeg, self.alpha, self.eta, self.CL)
        return alphaRes

    def getCDRes(self):
        # compute the residuals of the drag coefficient

        CDRes = numpy.zeros(self.numSeg)
        CDRes = mission.getcdres(self.numSeg, self.AR, self.e, self.CL, 
                                   self.CD)
        return CDRes

    def getEtaRes(self):
        # compute the residuals of eta from a given pitching moment coefficient

        etaRes = numpy.zeros(self.numSeg)
        etaRes = mission.getetares(self.numSeg, self.alpha, self.eta, self.CM)
        return etaRes

    def getTauRes(self):
        # compute the residuals of the throttle setting from a given Thrust

        tauRes = numpy.zeros(self.numSeg)
        tauRes = mission.gettaures(self.numSeg, self.cThrustSL, 
                                         self.x, self.h, self.tau, 
                                         self.Thrust)
        return tauRes

    def getWfRes(self):
        # compute the residuals of the fuel weight

        WfRes = numpy.zeros(self.numSeg)
        R1 = numpy.linspace(0.0,1.0,num=self.numInt)
        R2 = numpy.linspace(1.0,0.0,num=self.numInt)

        WfRes = mission.getwfres(self.numInt,self.numSeg,self.x,self.v,
                                 self.gamma,self.tau,self.Thrust,self.SFC,
                                 self.g,R1,R2,self.Wf)
        return WfRes

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
        # compute the residuals of the thrust

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
        # compute the derivatives of the residual of alpha

        self.dAlphadAlpha = numpy.zeros(self.numSeg)
        self.dAlphadEta = numpy.zeros(self.numSeg)
        self.dAlphadCL = numpy.zeros(self.numSeg)
        
        [self.dAlphadAlpha, self.dAlphadEta, 
         self.dAlphadCL] = mission.getdalpha(self.numSeg, self.alpha, self.eta, 
                                             self.CL)
        
        return

    def getDCD(self):
        # compute the derivatives of the residual of the drag coefficient

        self.dCDdAR = numpy.zeros(self.numSeg)
        self.dCDdE = numpy.zeros(self.numSeg)
        self.dCDdCL = numpy.zeros(self.numSeg)
        self.dCDdCD = numpy.zeros(self.numSeg)

        [self.dCDdAR, self.dCDdE, self.dCDdCL,
         self.dCDdCD] = mission.getdcd(self.numSeg, self.AR, self.e, self.CL, 
                                       self.CD)

        return

    def getDEta(self):
        # compute the derivatives of the residual of eta

        self.dEtadAlpha = numpy.zeros(self.numSeg)
        self.dEtadEta = numpy.zeros(self.numSeg)
        self.dEtadCM = numpy.zeros(self.numSeg)

        [self.dEtadAlpha, self.dEtadEta,
         self.dEtadCM] = mission.getdeta(self.numSeg, self.alpha, self.eta, 
                                         self.CM)

        return

    def getDTau(self):
        # compute the derivatives of the residual of the throttle setting

        self.dTaudCThrustSL = numpy.zeros(self.numSeg)
        self.dTaudH = numpy.zeros(self.numSeg)
        self.dTaudTau = numpy.zeros(self.numSeg)
        self.dTaudThrust = numpy.zeros(self.numSeg)

        [self.dTaudCThrustSL, self.dTaudH, 
         self.dTaudTau, 
         self.dTaudThrust] = mission.getdtau(self.numSeg, self.cThrustSL, 
                                             self.x, self.h, self.tau, 
                                             self.Thrust)

        return

    def getDWf(self):
        # compute the derivatives of the residual of the fuel weight

        self.dWfdCT1 = numpy.zeros(self.numSeg)
        self.dWfdCT2 = numpy.zeros(self.numSeg)
        self.dWfdThrust1 = numpy.zeros(self.numSeg)
        self.dWfdThrust2 = numpy.zeros(self.numSeg)
        self.dWfdX1 = numpy.zeros(self.numSeg)
        self.dWfdX2 = numpy.zeros(self.numSeg)
        self.dWfdV1 = numpy.zeros(self.numSeg)
        self.dWfdV2 = numpy.zeros(self.numSeg)
        self.dWfdGamma1 = numpy.zeros(self.numSeg)
        self.dWfdGamma2 = numpy.zeros(self.numSeg)
        self.dWfdWf1 = numpy.zeros(self.numSeg)
        self.dWfdWf2 = numpy.zeros(self.numSeg)

        R1 = numpy.linspace(0.0,1.0,num=self.numInt)
        R2 = numpy.linspace(1.0,0.0,num=self.numInt)

        [self.dWfdCT1, self.dWfdCT2, self.dWfdThrust1,
         self.dWfdThrust2, self.dWfdV1, self.dWfdV2, 
         self.dWfdGamma1, self.dWfdGamma2, self.dWfdWf1,
         self.dWfdWf2] = mission.getdwf(self.numSeg, self.numInt,
                                        self.g, self.x, self.v, self.gamma,
                                        self.tau, self.Thrust, self.SFC,
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
        # compute the derivatives of the residual of the thrust

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
        # compute the derivatives of the residual of the pitching moment coefficient

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

params = {'alpha':numpy.ones(100)*3.0*numpy.pi/180.0,
          'tau':numpy.ones(100)*0.3,
          'eta':numpy.zeros(100),
          'Wf':numpy.linspace(100000,0.0,100)*9.81,
          'S':525.0,
          'Wac':185953*9.81,
          'cThrust':1020000,
          'SFC':8.951e-6,
          'chord':8.15,
          'inertia':4.1e7,
          'AR':7.9,
          'e':0.8,
          'CL':numpy.ones(100)*0.3,
          'CD':numpy.ones(100)*0.0015,
          'CM':numpy.zeros(100),
          'Thrust':numpy.ones(100)*0.3*1020000,
          'numInt':100
          }

flightProblem = flightAnalysis(x,h,M,params)
xin = numpy.zeros([3*flightProblem.numSeg])
xin[0::3] = numpy.ones(100)*3*numpy.pi/180.0
xin[1::3] = numpy.ones(100)*0.3
xin[2::3] = numpy.zeros(100)

flightProblem.CM = numpy.ones(100)*0.1
flightProblem.getDEta()
print flightProblem.dEdGamma1
print flightProblem.dEdGamma2
print flightProblem.dEdGamma3
