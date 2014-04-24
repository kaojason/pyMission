from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
import numpy
import copy
from framework import *
import mission
import pylab


class Sys_h(ImplicitSystem):
    def _declare(self):
        h_IC = self.kwargs['h_IC']
        self.numPts = self.kwargs['numElem']+1
        self.numInt = self.kwargs['numInt']
        lower = 0.0
        iPts = range(self.numPts)

        self._declare_variable('h', size=self.numPts, 
                               val=h_IC, lower=lower)
        self._declare_argument('h_ends', indices=[0,1])
        #self._declare_argument(['v_ends',-1], indices=[0,1])
        self._declare_argument(['x',0], indices=[self.copy,
                                                 self.copy+1])
        #self._declare_argument(['h_dot',-1], indices=[0])
        self._declare_argument('Wf', indices=iPts)
        self._declare_argument('CT', indices=iPts)
        self._declare_argument('alpha', indices=iPts)
        self._declare_argument('CD', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        self._declare_argument('v', indices=iPts)
        self._declare_argument(['S', 0], indices=[0])
        self._declare_argument(['Wac', 0], indices=[0])

    def apply_F(self):
        self.climb = self.kwargs['climb']
        if self.climb == -1:
            pass
        else:
            #self.opt = self.kwargs['opt']
            p = self.vec['p']
            u = self.vec['u']
            f = self.vec['f']
            '''
            h_ends = p(['h_ends',-1])
            h = numpy.linspace(h_ends[0],h_ends[1],self.numPts)
            u(['h',-1])[:] = h
            '''
            h_ends = p('h_ends') * 1e3
            #v_ends = p(['v_ends',-1])
            #h_dot = p(['h_dot',-1])
            x_ends = p(['x',0]) * 1e6
            Wf = p('Wf') * 1e6
            CT = p('CT') * 1e-1
            alpha = p('alpha') * 1e-1
            CD = p('CD') * 1e-1
            rho = p('rho')
            v = p('v') * 1e2
            S = p(['S',0]) * 1e2
            h = u('h') * 1e3
            Wac = p(['Wac',0]) * 1e6

            h_new = mission.get_h(self.numPts-1, self.numInt, S, Wac, x_ends, 
                                  h_ends, Wf, CT, alpha, CD, rho, v)
            
            for i in xrange(self.numPts):
                if h_new[i] < 0:
                    h_new[i] = 0
                    #raise Exception('negative altitude!')
            
            #print "current h: ", h
            f('h')[:] = (h - h_new) / 1e3
    
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    
    def apply_dFdpu(self, arguments):
        self._lin_init()
        
        if self.climb == -1:
            pass

        else:
            p = self.vec['p']
            dp = self.vec['dp']
            df = self.vec['df']
            du = self.vec['du']

            numElem = self.numPts-1
            numInt = self.numInt
            S = p(['S',0]) * 1e2
            x_ends = p(['x',0]) * 1e6
            h_ends = p('h_ends') * 1e3
            Wf = p('Wf') * 1e6
            CT = p('CT') * 1e-1
            alpha = p('alpha') * 1e-1
            CD = p('CD') * 1e-1
            rho = p('rho')
            v = p('v') * 1e2
            Wac = p(['Wac',0]) * 1e6

            dS = dp(['S',0])
            dx_ends = dp(['x',0])
            dh_ends = dp('h_ends')
            dWf = dp('Wf')
            dCT = dp('CT')
            dalpha = dp('alpha')
            dCD = dp('CD')
            drho = dp('rho')
            dv = dp('v')
            dWac = dp(['Wac',0])

            [dhdS, dhdWac, dhdx0, dhdx1, dhdWf1, dhdWf2, dhdRho1, 
             dhdRho2, dhdV1, dhdV2,
             dhdCT1, dhdCT2, dhdA1, dhdA2, dhdCD1, 
             dhdCD2] = mission.get_h_d(numElem, numInt, S, Wac,  x_ends, 
                                       h_ends, Wf, CT, alpha, CD, rho, v)
             
            if self.mode == 'fwd':
                df('h')[:] = 0.0
                if self.get_id('S') in arguments:
                    df('h')[:] += dhdS * dS * 1e2/1e3
                if self.get_id('Wac') in arguments:
                    df('h')[:] += dhdWac * dWac * 1e6/1e3
                if self.get_id('x') in arguments:
                    df('h')[:] += (dhdx0*dx_ends[0] + 
                                   dhdx1*dx_ends[1]) * 1e6/1e3
                if self.get_id('Wf') in arguments:
                    df('h')[:] -= (dhdWf1+dhdWf2)*dWf * 1e6/1e3
                if self.get_id('rho') in arguments:
                    df('h')[:] -= (dhdRho1+dhdRho2)*drho * 1/1e3
                if self.get_id('v') in arguments:
                    df('h')[:] -= (dhdV1+dhdV2)*dv * 1e2/1e3
                if self.get_id('CT') in arguments:
                    df('h')[:] -= (dhdCT1+dhdCT2)*dCT * 1e-1/1e3
                if self.get_id('alpha') in arguments:
                    df('h')[:] -= (dhdA1+dhdA2)*dalpha * 1e-1/1e3
                if self.get_id('CD') in arguments:
                    df('h')[:] -= (dhdCD1+dhdCD2)*dCD * 1e-1/1e3
                if self.get_id('h_ends') in arguments:
                    df('h')[:] -= numpy.ones(self.numPts)* \
                        dh_ends[0]
                if self.get_id('h') in arguments:
                    df('h')[:] += du('h')

            if self.mode == 'rev':
                dx_ends[:] = 0.0
                dWf[:] = 0.0
                drho[:] = 0.0
                dv[:] = 0.0
                dalpha[:] = 0.0
                dpCT[:] = 0.0
                dpCD[:] = 0.0
                dS[:] = 0.0
                dWac[:] = 0.0
                du('h')[:] = 0.0
                if self.get_id('x') in arguments:
                    dx_ends[0] = dhdx0*df('h')
                    dx_ends[1] = dhdx1*df('h')
                if self.get_id('Wf') in arguments:
                    dWf[:] = (dhdWf1 + dhdWf2)*df('h')
                if self.get_id('rho') in arguments:
                    drho[:] = (dhdRho1 + dhdRho2)*df('h')
                if self.get_id('v') in arguments:
                    dv[:] = (dhdV1 + dhdV2)*df('h')
                if self.get_id('CT') in arguments:
                    dCT[:] = (dhdCT1 + dhdCT2)*df('h')
                if self.get_id('alpha') in arguments:
                    dalpha[:] = (dhdA1 + dhdA2)*df('h')
                if self.get_id('CD') in arguments:
                    dCD[:] = (dhdCD1 + dhdCD2)*df('h')

class Sys_SFC(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        SFC_IC = self.kwargs['SFC_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('SFC', size=numPts, val=SFC_IC) 
        self._declare_argument('h', indices=iPts)
        self._declare_argument(['SFCSL',0], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p('h') * 1e3
        SFCSL = p(['SFCSL',0]) * 1e-6

        SFC = mission.get_sfc(self.numElem, SFCSL, h)
        u('SFC')[:] = SFC / 1e-6
    
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    
    def apply_dGdp(self, arguments):
        mode = self.mode
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        numPts = self.numElem+1
        dSFC_dh = numpy.zeros(numPts)        
        dSFC_dh = mission.get_sfc_d(self.numElem)

        dh = dp('h')
        dSFCSL = dp('SFCSL')

        if self.mode == 'fwd':
            dg('SFC')[:] = 0.0
            if self.get_id('h') in arguments:
                dg('SFC')[:] += (dSFC_dh * dh) * 1e3/1e-6
            if self.get_id('SFCSL') in arguments:
                dg('SFC')[:] += 1.0

        if self.mode == 'rev':
            dh[:] = 0.0
            dSFCSL[:] = 0.0
            du('SFC')[:] = 0.0
            if self.get_id('h') in arguments:
                dh[:] = dSFC_dh * dg('SFC')
            if self.get_id('SFCSL') in arguments:
                dSFCSL[:] = 1.0
    
class Sys_gamma(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        gamma_IC = self.kwargs['gamma_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('gamma', size=numPts, val=gamma_IC)
        self._declare_argument('h', indices=iPts)
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        
    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p('h') * 1e3
        x_ends = p(['x',0]) * 1e6
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)
        gamma = mission.get_gamma(self.numElem, h, x_int)
        u(['gamma',-1])[:] = gamma / 1e-1
    
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    
    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        numPts = self.numElem+1
        dgamma_dh1 = numpy.zeros(numPts)
        dgamma_dh2 = numpy.zeros(numPts)
        dgamma_dh3 = numpy.zeros(numPts)
        dgamma_dh4 = numpy.zeros(numPts)
        dgamma_dh5 = numpy.zeros(numPts)
        dgamma_dx = numpy.zeros(numPts)
        h = p('h') * 1e3
        x_ends = p(['x',0]) * 1e6
        x = numpy.linspace(x_ends[0], x_ends[1], numPts)

        dh = dp('h')
        dx = dp(['x',0])

        [dgamma_dh1, dgamma_dh2, dgamma_dh3, dgamma_dh4, dgamma_dh5,
         dgamma_dx0, dgamma_dx1] = mission.get_gamma_d(numPts-1, h, x)

        if self.mode == 'fwd':
            dg('gamma')[:] = 0.0
            if self.get_id('h') in arguments:
                dg('gamma')[:] += (dgamma_dh1*dh + 
                                   dgamma_dh2*dh + 
                                   dgamma_dh3*dh + 
                                   dgamma_dh4*dh + 
                                   dgamma_dh5*dh) * \
                                   1e3/1e-1
            if self.get_id('x') in arguments:
                dg('gamma')[:] += (dgamma_dx0*dx[0] + 
                                   dgamma_dx1*dx[1])  * \
                                   1e6/1e-1

        if self.mode == 'rev':
            dh[:] = 0.0
            dx[:] = 0.0
            du('gamma')[:] = 0.0
            if self.get_id('h') in arguments:
                dh[:] = (dgamma_dh1 + dgamma_dh2 + dgamma_dh3 +
                         dgamma_dh4 + dgamma_dh5)*dg('gamma')
            if self.get_id('x') in arguments:
                dx[0] = dgamma_dx0 * dg('gamma')
                dx[1] = dgamma_dx1 * dg('gamma')
            if self.get_id('gamma') in arguments:
                du('gamma')[:] = 0.0
   
class Sys_Temp(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        Temp_IC = self.kwargs['Temp_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('Temp', size=numPts, val=Temp_IC,
                               lower=0.001)
        self._declare_argument('h', indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p('h') * 1e3
        #print  " current h-temp: ", h
        Temp = mission.get_temp(self.numElem, h)
        #print " current Temp: ", Temp
        u('Temp')[:] = Temp / 1e2
    
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    
    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        numPts = self.numElem+1
        dTemp_dh = numpy.zeros(numPts)
        h = p('h') * 1e3

        dh = dp('h')
        
        dTemp_dh = mission.get_temp_d(self.numElem, h)
        
        if self.mode == 'fwd':
            dg('Temp')[:] = 0.0
            if self.get_id('h') in arguments:
                dg('Temp')[:] += (dTemp_dh * dh) * \
                    1e3/1e2
        if self.mode == 'rev':
            dh[:] = 0.0
            du('Temp')[:] = 0.0
            if self.get_id('h') in arguments:
                dh[:] = dTemp_dh * dg('Temp')
            if self.get_id('Temp') in arguments:
                du('Temp')[:] = 0.0
    
class Sys_rho(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        rho_IC = self.kwargs['rho_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('rho', size=numPts, val=rho_IC,
                               lower=0.001)
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument('Temp', indices=iPts)
        
    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        g = p(['g',0])
        Temp = p('Temp') * 1e2
        rho = mission.get_rho(self.numElem, g, Temp)

        #print "current rho: ", rho
        u('rho')[:] = rho
    
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    
    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        g = 9.81
        Temp = p('Temp') * 1e2
        numPts = self.numElem+1
        drho_dTemp = numpy.zeros(numPts)

        dTemp = dp('Temp')

        drho_dTemp = mission.get_rho_d(self.numElem, g, Temp)

        if self.mode == 'fwd':
            dg('rho')[:] = 0.0
            if self.get_id('Temp') in arguments:
                dg('rho')[:] += (drho_dTemp * dTemp) * 1e2
        if self.mode == 'rev':
            dTemp[:] = 0.0
            du('rho')[:] = 0.0
            if self.get_id('Temp') in arguments:
                dTemp[:] = drho_dTemp * dg('rho')
            if self.get_id('rho') in arguments:
                du('rho')[:] = 0.0
    
class Sys_v(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        v_IC = self.kwargs['v_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        ##print "v_IC: ", v_IC
        self._declare_variable('v', size=numPts, val=v_IC,
                               lower=0.0)
        self._declare_argument('Temp', indices=iPts)
        self._declare_argument('v_ends', indices=[0,1])
        self._declare_argument('M_ends', indices=[0,1])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        v_ends = p('v_ends')
        M_ends = p('M_ends')
        Temp = p('Temp') * 1e2
        v_copy = copy.copy(v_ends)

        if v_copy[0] != -1:
            v_copy[0] *= 1e2
        if v_copy[1] != -1:
            v_copy[1] *= 1e2

        if v_copy[0] == -1:
            v_copy[0] = M_ends[0]*numpy.sqrt(1.4*288*Temp[0])
        if v_copy[1] == -1:
            v_copy[1] = M_ends[1]*numpy.sqrt(1.4*288*Temp[self.numElem])

        v = mission.get_v(self.numElem, v_copy)
        #print "current v: ", v
        u('v')[:] = v / 1e2
    
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    
    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        numPts = self.numElem+1
        dv_dvEnds0 = numpy.zeros(numPts)
        dv_dvEnds1 = numpy.zeros(numPts)
        dv_dMEnds0 = numpy.zeros(numPts)
        dv_dMEnds1 = numpy.zeros(numPts)
        dv_dTemp = numpy.zeros(numPts)
        M_ends = p('M_ends')
        Temp = p('Temp') * 1e2
        v_ends = p('v_ends')
        v_copy = copy.copy(v_ends)

        if v_copy[0] != -1:
            v_copy[0] *= 1e2
        if v_copy[1] != -1:
            v_copy[1] *= 1e2

        dM = dp('M_ends')
        dv = dp('v_ends')
        dT = dp('Temp')

        [dv_dvEnds0, dv_dvEnds1] = mission.get_v_d(self.numElem, v_copy)

        if v_copy[0] == -1:
            dv_dMEnds0 = dv_dvEnds0*numpy.sqrt(1.4*288*Temp[0])
            dv_dvEnds0 = numpy.zeros(numPts)
        if v_copy[1] == -1:
            dv_dMEnds1 = dv_dvEnds1*numpy.sqrt(1.4*288*Temp[numPts-1])
            dv_dvEnds1 = numpy.zeros(numPts)
        # FIX TEMP DEPENDENCE FOR MACH SPECIFICATION
        if self.mode == 'fwd':
            dg('v')[:] = 0.0
            if self.get_id('v_ends') in arguments:
                dg('v')[:] += (dv_dvEnds0 * dv[0] +
                               dv_dvEnds1 * dv[1])
            if self.get_id('M_ends') in arguments:
                dg('v')[:] += (dv_dMEnds0 * dM[0] +
                               dv_dMEnds1 * dM[1]) * \
                               1/1e2
        if self.mode == 'rev':
            dM[:] = 0.0
            dv[:] = 0.0
            dT[:] = 0.0
            du('v')[:] = 0.0
            if self.get_id('v_ends') in arguments:
                dv[:] = (dv_dvEnds0 + dv_dvEnds1)*dg('v')
            if self.get_id('M_ends') in arguments:
                dM[:] = (dv_dMEnds0 + dv_dMEnds1)*dg('v')
            if self.get_id('v') in arguments:
                du('v')[:] = 0.0

class Sys_CT(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        t_IC = self.kwargs['t_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('CT', size=numPts, val=t_IC)
        self._declare_argument('h', indices=iPts)
        self._declare_argument('tau_init', indices=[0])
        self._declare_argument('rho', indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument('v', indices=iPts)
        self._declare_argument(['cThrustSL',0], indices=[0])

    def apply_G(self):
        climb = self.kwargs['climb']
        
        if climb == -1:
            pass
        else:
            numElem = self.numElem
            p = self.vec['p']
            u = self.vec['u']
            
            tauInit = p('tau_init')
            rho = p('rho')
            v = p('v') * 1e2
            S = p(['S',0]) * 1e2
            cThrustSL = p(['cThrustSL',0]) * 1e6
            h = p('h') * 1e3

            CT = mission.get_ct_init(numElem, tauInit, S, cThrustSL, rho, v, h)

            #print "current CT: ", CT
            u('CT')[:] = CT / 1e-1

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dGdp(self, arguments):
        climb = self.kwargs['climb']

        if climb == -1:
            pass
        else:
            p = self.vec['p']
            dp = self.vec['dp']
            dg = self.vec['dg']
            du = self.vec['du']
        
            tau = p('tau_init')
            S = p(['S',0]) * 1e2
            cThrustSL = p(['cThrustSL',0]) * 1e6
            rho = p('rho')
            v = p('v') * 1e2
            h = p('h') * 1e3

            dtau = dp('tau_init')
            dS = dp(['S',0])
            dcThrustSL = dp(['cThrustSL',0])
            drho = dp('rho')
            dv = dp('v')
            dh = dp('h')
            
            [dCTdcThrustSL, dCTdH, dCTdRho, dCTdV, dCTdS,
             dCTdTau] = mission.get_ct_init_d(self.numElem, tau, S, 
                                              cThrustSL, rho, v, h)
        
            if self.mode == 'fwd':
                dg('CT')[:] = 0.0
                if self.get_id('cThrustSL') in arguments:
                    dg('CT')[:] += (dCTdcThrustSL *
                                    dcThrustSL) * 1e6/1e-1
                if self.get_id('h') in arguments:
                    dg('CT')[:] += (dCTdH * dh) * 1e3/1e-1
                if self.get_id('rho') in arguments:
                    dg('CT')[:] += (dCTdRho * drho) / 1e-1
                if self.get_id('v') in arguments:
                    dg('CT')[:] += (dCTdV * dv) * 1e2/1e-1
                if self.get_id('S') in arguments:
                    dg('CT')[:] += (dCTdS * dS) * 1e2/1e-1
                if self.get_id('tau_init') in arguments:
                    dg('CT')[:] += (dCTdTau * dtau) / 1e-1
            if self.mode == 'rev':
                dtau[:] = 0.0
                dS[:] = 0.0
                dcThrustSL[:] = 0.0
                drho[:] = 0.0
                dv[:] = 0.0
                dh[:] = 0.0
                du('CT')[:] = 0.0
                if self.get_id('cThrustSL') in arguments:
                    dcThrustSL[:] = dCTdcThrustSL * dg('CT')
                if self.get_id('h') in arguments:
                    dh[:] = dCTdH * dg('CT')
                if self.get_id('rho') in arguments:
                    drho[:] = dCTdRho * dg('CT')
                if self.get_id('v') in arguments:
                    dv[:] = dCTdV * dg('CT')
                if self.get_id('S') in arguments:
                    dS[:] = dCTdS * dg('CT')
                if self.get_id('tau_init') in arguments:
                    dtau[:] = dCTdTau * dg('CT')

class Sys_tau(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('tau', size=numPts)
        self._declare_argument('CT', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        self._declare_argument('v', indices=iPts)
        self._declare_argument('h', indices=iPts)
        self._declare_argument(['cThrustSL',0], indices=[0])
        self._declare_argument(['S',0], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        CT = p('CT') * 1e-1
        rho = p('rho')
        v = p('v') * 1e2
        h = p('h') * 1e3
        cThrustSL = p(['cThrustSL',0]) * 1e6
        S = p(['S',0]) * 1e2

        tau = mission.get_tau(self.numElem, cThrustSL, S, CT, rho, v, h)
        #print "current tau: ", tau
        u(['tau',-1])[:] = tau

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']

        cThrustSL = p(['cThrustSL',0]) * 1e6
        S = p(['S',0]) * 1e2
        h = p('h') * 1e3
        CT = p('CT') * 1e-1
        rho = p('rho')
        v = p('v') * 1e2

        dcThrustSL = dp(['cThrustSL', 0])
        dS = dp(['S', 0])
        dh = dp('h')
        dCT = dp('CT')
        drho = dp('rho')
        dv = dp('v')
        
        [dtdcThrustSL, dtdh, dtdCT, dtdRho, dtdV,
         dtdS] = mission.get_tau_d(self.numElem, cThrustSL, S, h, CT,
                                   rho, v)

        if self.mode == 'fwd':
            dg('tau')[:] = 0.0
            if self.get_id('cThrustSL') in arguments:
                dg('tau')[:] += (dtdcThrustSL * 
                                 dcThrustSL) * 1e6
            if self.get_id('h') in arguments:
                dg('tau')[:] += (dtdh * dh) * 1e3
            if self.get_id('CT') in arguments:
                dg('tau')[:] += (dtdCT * dCT) * 1e-1
            if self.get_id('rho') in arguments:
                dg('tau')[:] += dtdRho * drho
            if self.get_id('v') in arguments:
                dg('tau')[:] += (dtdV * dv) * 1e2
            if self.get_id('S') in arguments:
                dg('tau')[:] += (dtdS * dS) * 1e2
        if self.mode == 'rev':
            dcThrustSL[:] = 0.0
            dh[:] = 0.0
            dCT[:] = 0.0
            drho[:] = 0.0
            dv[:] = 0.0
            dS[:] = 0.0
            du('tau')[:] = 0.0
            if self.get_id('cThrustSL') in arguments:
                dcThrustSL[:] = dtdcThrustSL * dg('tau')
            if self.get_id('h') in arguments:
                dh[:] = dtdh * dg('tau')
            if self.get_id('CT') in arguments:
                dCT[:] = dtdCT * dg('tau')
            if self.get_id('rho') in arguments:
                drho[:] = dtdRho * dg('tau')
            if self.get_id('v') in arguments:
                dv[:] = dtdV * dg('tau')
            if self.get_id('S') in arguments:
                dS[:] = dtdS * dg('tau')
            if self.get_id('tau') in arguments:
                du('tau')[:] = 0.0

class Sys_CL(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CL_IC = self.kwargs['CL_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('CL', size=numPts, val=CL_IC)
        self._declare_argument('Wf', indices=iPts)
        self._declare_argument('gamma', indices=iPts)
        self._declare_argument('CT', indices=iPts)
        self._declare_argument('alpha', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        self._declare_argument('v', indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['Wac',0], indices=[0])
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        
    def apply_F(self):
        self.numInt = self.kwargs['numInt']

        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        Wf = p('Wf') * 1e6
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        alpha = p('alpha') * 1e-1
        rho = p('rho')
        v = p('v') * 1e2
        S = p(['S',0]) * 1e2
        Wac = p(['Wac',0]) * 1e6
        g = 9.81
        CL = u('CL')
        x_ends = p(['x',0]) * 1e6
        R = f('CL')
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        
        CLRes = numpy.zeros(self.numElem+1)
        R1 = numpy.linspace(1.0, 0.0, self.numInt)
        R2 = numpy.linspace(0.0, 1.0, self.numInt)

        for i in xrange(self.numElem):
            rhoTemp = numpy.linspace(rho[i], rho[i+1], self.numInt)
            vTemp = numpy.linspace(v[i], v[i+1], self.numInt)
            gammaTemp = numpy.linspace(gamma[i], gamma[i+1], self.numInt)
            CLTemp = numpy.linspace(CL[i], CL[i+1], self.numInt)
            WTemp = numpy.linspace(Wac+Wf[i], Wac+Wf[i+1], self.numInt)
            aTemp = numpy.linspace(alpha[i], alpha[i+1], self.numInt)
            CTTemp = numpy.linspace(CT[i], CT[i+1], self.numInt)

            deltax = (x_int[i+1]-x_int[i])/self.numInt
            QTemp = 0.5*rhoTemp*vTemp*vTemp*S
            dxTemp = numpy.ones(self.numInt) * deltax
            dxTemp[0] = 0.5*deltax
            dxTemp[-1] = 0.5*deltax

            cosGamma = numpy.cos(gammaTemp)
            sinAlpha = numpy.sin(aTemp)

            CLResTemp = -QTemp*CLTemp + WTemp*cosGamma - CTTemp*QTemp*sinAlpha
            CLRes[i] += numpy.sum(CLResTemp * R1 * dxTemp)
            CLRes[i+1] += numpy.sum(CLResTemp * R2 * dxTemp)
        R[:] = CLRes / 1e7
            
        """
        CLRes = mission.get_cl(self.numElem, self.numInt, Wac, S, g, 
                               x_int, v, rho, CL, Wf, gamma, CT, alpha)

        QTemp = numpy.zeros(self.numInt)
        CLRes = numpy.zeros(self.numElem+1)

        R1 = numpy.linspace(1.0, 0.0, self.numInt)
        R2 = numpy.linspace(0.0, 1.0, self.numInt)
        
        for i in xrange(self.numElem):
            rhoTemp = numpy.linspace(rho[i], rho[i+1], self.numInt)
            vTemp = numpy.linspace(v[i], v[i+1], self.numInt)
            deltax = (x_int[i+1]-x_int[i])/self.numInt
            
            for j in xrange(self.numInt):
                QTemp[j] = 0.5*rhoTemp[j]*vTemp[j]*vTemp[j]*S
            xTemp = numpy.ones(self.numInt)*deltax
            xTemp[0] = 0.5*deltax
            xTemp[-1] = 0.5*deltax
            
            gammaTemp = numpy.linspace(gamma[i], gamma[i+1], self.numInt)
            for j in xrange(self.numInt):
                gammaTemp[j] = gamma[i] * R1[j] + gamma[i+1] * R2[j]
            cosGamma = numpy.cos(gammaTemp)

            CLTemp = numpy.linspace(CL[i], CL[i+1], self.numInt)
            WTemp = numpy.linspace(Wf[i], Wf[i+1], self.numInt)
            aTemp = numpy.linspace(alpha[i], alpha[i+1], self.numInt)
            CTTemp = numpy.linspace(CT[i], CT[i+1], self.numInt)
            
            for j in xrange(self.numInt):
                WTemp[j] += Wac
                '''
                CLRes[i] += (-QTemp[j]*CLTemp[j]+WTemp[j]*cosGamma[j] -
                              numpy.sin(aTemp[j])*CTTemp[j]*
                              QTemp[j])*R1[j]*xTemp[j]

                CLRes[i+1] += (-QTemp[j]*CLTemp[j]+WTemp[j]*cosGamma[j] -
                                numpy.sin(aTemp[j])*CTTemp[j]*
                                QTemp[j])*R2[j]*xTemp[j]
                '''
                CLRes[i] += (CLTemp[j] - numpy.sin(vTemp[j])) * R1[j] * xTemp[j]
                CLRes[i+1] += (CLTemp[j] - numpy.sin(vTemp[j])) * R2[j] * xTemp[j]

        #print "current CL: ", CL
        #print "current CLRes: ", CLRes
        f('CL')[:] = CLRes
        """

    def apply_dFdpu1(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dFdpu0(self, args):
        self.numInt = self.kwargs['numInt']

        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']

        Wf = p('Wf') * 1e6
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        alpha = p('alpha') * 1e-1
        rho = p('rho')
        v = p('v') * 1e2
        S = p(['S',0]) * 1e2
        Wac = p(['Wac',0]) * 1e6
        g = 9.81 
        CL = u('CL')
        x_ends = p(['x',0]) * 1e6
        R = f('CL')
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        dWf = dp('Wf') * 1e6 / 1e7
        dgamma = dp('gamma') * 1e-1 / 1e7
        dCT = dp('CT') * 1e-1 / 1e7
        dalpha = dp('alpha') * 1e-1 / 1e7
        drho = dp('rho') / 1e7
        dv = dp('v') * 1e2 / 1e7
        dS = dp(['S',0]) * 1e2 / 1e7
        dWac = dp(['Wac',0]) * 1e6 / 1e7
        dCL = du('CL') / 1e7
        dx_ends = dp(['x',0]) * 1e6 / 1e7
        dR = df('CL')
        dx_int = 0

        CLRes = numpy.zeros(self.numElem+1)
        R1 = numpy.linspace(1.0, 0.0, self.numInt)
        R2 = numpy.linspace(0.0, 1.0, self.numInt)

        if self.mode == 'fwd':
            dR[:] = 0.0
        elif self.mode == 'rev':
            jason = 1
            # INITIALIZE dARGS!!!

        for i in xrange(self.numElem):
            rhoTemp = numpy.linspace(rho[i], rho[i+1], self.numInt)
            vTemp = numpy.linspace(v[i], v[i+1], self.numInt)
            gammaTemp = numpy.linspace(gamma[i], gamma[i+1], self.numInt)
            CLTemp = numpy.linspace(CL[i], CL[i+1], self.numInt)
            WTemp = numpy.linspace(Wac+Wf[i], Wac+Wf[i+1], self.numInt)
            aTemp = numpy.linspace(alpha[i], alpha[i+1], self.numInt)
            CTTemp = numpy.linspace(CT[i], CT[i+1], self.numInt)

            deltax = (x_int[i+1]-x_int[i])/self.numInt
            QTemp = 0.5*rhoTemp*vTemp*vTemp*S
            dxTemp = numpy.ones(self.numInt) * deltax
            dxTemp[0] /= 2.0
            dxTemp[-1] /= 2.0

            cosGamma = numpy.cos(gammaTemp)
            sinAlpha = numpy.sin(aTemp)

            dQ_drho = 0.5*vTemp*vTemp*S
            dQ_dv = rhoTemp*vTemp*S
            dQ_dS = 0.5*rhoTemp*vTemp*vTemp

            dcos_dgamma = -numpy.sin(gammaTemp)
            dsin_dalpha = numpy.cos(aTemp)

            dR_dQ = -CLTemp - CTTemp*sinAlpha
            dR_dCL = -QTemp
            dR_dW = numpy.array(cosGamma)
            dR_dcos = numpy.array(WTemp)
            dR_dCT = -QTemp*sinAlpha
            dR_dsin = -CTTemp*QTemp

            dR_drho = dR_dQ * dQ_drho
            dR_dv = dR_dQ * dQ_dv
            dR_dS = dR_dQ * dQ_dS
            dR_dWac = numpy.array(dR_dW)
            dR_dWf = numpy.array(dR_dW)
            dR_dgamma = dR_dcos * dcos_dgamma
            dR_dalpha = dR_dsin * dsin_dalpha

            if self.mode == 'fwd':
                if self.get_id('rho') in args:
                    dR[i] += numpy.sum(dR_drho * R1 * dxTemp * R1) * drho[i]
                    dR[i+1] += numpy.sum(dR_drho * R2 * dxTemp * R2) * drho[i+1]
                    dR[i] += numpy.sum(dR_drho * R1 * dxTemp * R2) * drho[i+1]
                    dR[i+1] += numpy.sum(dR_drho * R2 * dxTemp * R1) * drho[i]
                if self.get_id('v') in args:
                    dR[i] += numpy.sum(dR_dv * R1 * dxTemp * R1) * dv[i]
                    dR[i+1] += numpy.sum(dR_dv * R2 * dxTemp * R2) * dv[i+1]
                    dR[i] += numpy.sum(dR_dv * R1 * dxTemp * R2) * dv[i+1]
                    dR[i+1] += numpy.sum(dR_dv * R2 * dxTemp * R1) * dv[i]
                if self.get_id('S') in args:
                    dR[i] += numpy.sum(dR_dS * R1 * dxTemp) * dS[0]
                    dR[i+1] += numpy.sum(dR_dS * R2 * dxTemp) * dS[0]
                if self.get_id('Wac') in args:
                    dR[i] += numpy.sum(dR_dWac * R1 * dxTemp) * dWac[0]
                    dR[i+1] += numpy.sum(dR_dWac * R2 * dxTemp) * dWac[0]
                if self.get_id('Wf') in args:
                    dR[i] += numpy.sum(dR_dWf * R1 * dxTemp * R1) * dWf[i]
                    dR[i+1] += numpy.sum(dR_dWf * R2 * dxTemp * R2) * dWf[i+1]
                    dR[i] += numpy.sum(dR_dWf * R1 * dxTemp * R2) * dWf[i+1]
                    dR[i+1] += numpy.sum(dR_dWf * R2 * dxTemp * R1) * dWf[i]
                if self.get_id('gamma') in args:
                    dR[i] += numpy.sum(dR_dgamma * R1 * dxTemp * R1) * dgamma[i]
                    dR[i+1] += numpy.sum(dR_dgamma * R2 * dxTemp * R2) * dgamma[i+1]
                    dR[i] += numpy.sum(dR_dgamma * R1 * dxTemp * R2) * dgamma[i+1]
                    dR[i+1] += numpy.sum(dR_dgamma * R2 * dxTemp * R1) * dgamma[i]
                if self.get_id('CT') in args:
                    dR[i] += numpy.sum(dR_dCT * R1 * dxTemp * R1) * dCT[i]
                    dR[i+1] += numpy.sum(dR_dCT * R2 * dxTemp * R2) * dCT[i+1]
                    dR[i] += numpy.sum(dR_dCT * R1 * dxTemp * R2) * dCT[i+1]
                    dR[i+1] += numpy.sum(dR_dCT * R2 * dxTemp * R1) * dCT[i]
                if self.get_id('alpha') in args:
                    dR[i] += numpy.sum(dR_dalpha * R1 * dxTemp * R1) * dalpha[i]
                    dR[i+1] += numpy.sum(dR_dalpha * R2 * dxTemp * R2) * dalpha[i+1]
                    dR[i] += numpy.sum(dR_dalpha * R1 * dxTemp * R2) * dalpha[i+1]
                    dR[i+1] += numpy.sum(dR_dalpha * R2 * dxTemp * R1) * dalpha[i]
                if self.get_id('CL') in args:
                    dR[i] += numpy.sum(dR_dCL * R1 * dxTemp * R1) * dCL[i]
                    dR[i+1] += numpy.sum(dR_dCL * R2 * dxTemp * R2) * dCL[i+1]
                    dR[i] += numpy.sum(dR_dCL * R1 * dxTemp * R2) * dCL[i+1]
                    dR[i+1] += numpy.sum(dR_dCL * R2 * dxTemp * R1) * dCL[i]

    def apply_dFdpu(self, arguments):
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']

        Wac = p(['Wac',0]) * 1e6
        S = p(['S',0]) * 1e2
        x_ends = p(['x',0]) * 1e6
        v = p('v') * 1e2
        rho = p('rho')
        CL = u('CL')
        Wf = p('Wf') * 1e6
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        alpha = p('alpha') * 1e-1
        g = 9.81
        x = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        dWac = dp(['Wac', 0])
        dS = dp(['S', 0])
        dx = dp(['x', 0])
        dv = dp('v')
        drho = dp('rho')
        dWf = dp('Wf')
        dgamma = dp('gamma')
        dCT = dp('CT')
        dalpha = dp('alpha')

        [dCLdWac, dCLdS, dCLdV1, dCLdV2, dCLdV3, dCLdRho1, dCLdRho2,
         dCLdRho3, dCLdCL1, dCLdCL2, dCLdCL3, dCLdWf1, dCLdWf2, dCLdWf3,
         dCLdGamma1, dCLdGamma2, dCLdGamma3, dCLdThrust1, dCLdThrust2,
         dCLdThrust3, dCLdAlpha1, dCLdAlpha2,
         dCLdAlpha3] = mission.get_cl_d(self.numElem, self.numInt, Wac,
                                        S, g, x, v, rho, CL, Wf,
                                        gamma, CT, alpha)
        '''
        R1 = numpy.linspace(1.0, 0.0, self.numInt)
        R2 = numpy.linspace(0.0, 1.0, self.numInt)

        for i in xrange(self.numElem+1):
            dCLdGamma1[i] = 0.0
            dCLdGamma2[i] = 0.0
            dCLdGamma3[i] = 0.0
        
            dCLdV1[i] = 0.0
            dCLdV2[i] = 0.0
            dCLdV3[i] = 0.0
            dCLdThrust1[i] = 0.0
            dCLdThrust2[i] = 0.0
            dCLdThrust3[i] = 0.0

        QTemp = numpy.zeros(self.numInt)
        dVTemp1 = numpy.zeros(self.numInt)
        dVTemp2 = numpy.zeros(self.numInt)
        dQTempdV1 = numpy.zeros(self.numInt)
        dQTempdV2 = numpy.zeros(self.numInt)
        for i in xrange(self.numElem):
            rhoTemp = numpy.linspace(rho[i], rho[i+1], self.numInt)
            vTemp = numpy.linspace(v[i], v[i+1], self.numInt)
            deltax = (x[i+1]-x[i])/self.numInt
            
            for j in xrange(self.numInt):
                QTemp[j] = 0.5*rhoTemp[j]*vTemp[j]*vTemp[j]*S
            xTemp = numpy.ones(self.numInt)*deltax
            xTemp[0] = 0.5*deltax
            xTemp[-1] = 0.5*deltax

            for j in xrange(self.numInt):
                dVTemp1[j] = R1[j]
                dVTemp2[j] = R2[j]
            
                dQTempdV1[j] = rhoTemp[j]*vTemp[j]*dVTemp1[j]*S
                dQTempdV2[j] = rhoTemp[j]*vTemp[j]*dVTemp2[j]*S
            
            gammaTemp = numpy.linspace(gamma[i], gamma[i+1], self.numInt)
            for j in xrange(self.numInt):
                gammaTemp[j] = gamma[i] * R1[j] + gamma[i+1] * R2[j]
            dGammaTemp1 = R1
            dGammaTemp2 = R2
            cosGamma = numpy.cos(gammaTemp)
            dCosGamma1 = -numpy.sin(gammaTemp)*dGammaTemp1
            dCosGamma2 = -numpy.sin(gammaTemp)*dGammaTemp2
            
            CLTemp = numpy.linspace(CL[i], CL[i+1], self.numInt)
            WTemp = numpy.linspace(Wf[i], Wf[i+1], self.numInt)
            aTemp = numpy.linspace(alpha[i], alpha[i+1], self.numInt)
            CTTemp = numpy.linspace(CT[i], CT[i+1], self.numInt)
            dCTTemp1 = R1
            dCTTemp2 = R2

            for j in xrange(self.numInt):
                WTemp[j] += Wac
                dCLdGamma2[i] += R1[j] * xTemp[j] * WTemp[j] * dCosGamma1[j]
                dCLdGamma3[i] += R1[j] * xTemp[j] * WTemp[j] * dCosGamma2[j]
                dCLdGamma1[i+1] += R2[j] * xTemp[j] * WTemp[j] * dCosGamma1[j]
                dCLdGamma2[i+1] += R2[j] * xTemp[j] * WTemp[j] * dCosGamma2[j]
                
                dCLdV2[i] -= R1[j] * xTemp[j] * (CLTemp[j] * dQTempdV1[j]
                                                 + numpy.sin(aTemp[j]) *
                                                 CTTemp[j] * dQTempdV1[j])
                dCLdV3[i] -= R1[j] * xTemp[j] * (CLTemp[j] * dQTempdV2[j]
                                                 + numpy.sin(aTemp[j]) * 
                                                 CTTemp[j] * dQTempdV2[j])
                dCLdV1[i+1] -= R2[j] * xTemp[j] * (CLTemp[j] * dQTempdV1[j]
                                                   + numpy.sin(aTemp[j]) * 
                                                   CTTemp[j] * dQTempdV1[j])
                dCLdV2[i+1] -= R2[j] * xTemp[j] * (CLTemp[j] * dQTempdV2[j]
                                                   + numpy.sin(aTemp[j]) * 
                                                   CTTemp[j] * dQTempdV2[j])
                
                dCLdV2[i] -= numpy.cos(vTemp[j])*R1[j]*R1[j]*xTemp[j]
                dCLdV3[i] -= numpy.cos(vTemp[j])*R2[j]*R1[j]*xTemp[j]
                dCLdV1[i+1] -= numpy.cos(vTemp[j])*R1[j]*R2[j]*xTemp[j]
                dCLdV2[i+1] -= numpy.cos(vTemp[j])*R2[j]*R2[j]*xTemp[j]
                
                dCLdThrust2[i] -= (R1[j] * xTemp[j] * numpy.sin(aTemp[j]) *
                                   QTemp[j] * dCTTemp1[j])
                dCLdThrust3[i] -= (R1[j] * xTemp[j] * numpy.sin(aTemp[j]) *
                                   QTemp[j] * dCTTemp2[j])
                dCLdThrust1[i+1] -= (R2[j] * xTemp[j] * numpy.sin(aTemp[j])
                                     * QTemp[j] * dCTTemp1[j])
                dCLdThrust2[i+1] -= (R2[j] * xTemp[j] * numpy.sin(aTemp[j])
                                     * QTemp[j] * dCTTemp2[j])
                '''
        if self.mode == 'fwd':
            df('CL')[:] = 0.0
            if self.get_id('Wac') in arguments:
                df('CL')[:] += (dCLdWac * dWac) * 1e6 / 1e7
            if self.get_id('S') in arguments:
                df('CL')[:] += dCLdS * dS * 1e2 / 1e7
            if self.get_id('v') in arguments:
                df('CL')[:] += ((dCLdV1 + dCLdV2 + dCLdV3) * 
                                dv) * 1e2 / 1e7
            if self.get_id('rho') in arguments:
                df('CL')[:] += ((dCLdRho1 + dCLdRho2 + dCLdRho3) *
                                drho) / 1e7
            if self.get_id('Wf') in arguments:
                df('CL')[:] += ((dCLdWf1 + dCLdWf2 + dCLdWf3) *
                                dWf) * 1e6 / 1e7
            if self.get_id('gamma') in arguments:
                df('CL')[:] += ((dCLdGamma1 + dCLdGamma2 + dCLdGamma3) *
                                dgamma) * 1e-1 / 1e7
            if self.get_id('CT') in arguments:
                df('CL')[:] += (dCLdThrust1 + dCLdThrust2 +
                                dCLdThrust3) * dCT * 1e-1 / 1e7
            if self.get_id('alpha') in arguments:
                df('CL')[:] += (dCLdAlpha1 + dCLdAlpha2 +
                                dCLdAlpha3) * dalpha * 1e-1 / 1e7
            if self.get_id('CL') in arguments:
                df('CL')[:] += ((dCLdCL1 + dCLdCL2 + dCLdCL3) *
                                du('CL'))  / 1e7
        if self.mode == 'rev':
            dWac[:] = 0.0
            dS[:] = 0.0
            dv[:] = 0.0
            drho[:] = 0.0
            dWf[:] = 0.0
            dgamma[:] = 0.0
            dCT[:] = 0.0
            dalpha[:] = 0.0
            du('CL')[:] = 0.0
            if self.get_id('Wac') in arguments:
                dWac[:] = dCLdWac * df('CL')
            if self.get_id('S') in arguments:
                dS[:] = dCLdS * df('CL')
            if self.get_id('v') in arguments:
                dv[:] = (dCLdV1 + dCLdV2 + dCLdV3) * df('CL')
            if self.get_id('rho') in arguments:
                drho[:] = (dCLdRho1 + dCLdRho2 + dCLdRho3) * \
                    df('CL')
            if self.get_id('CL') in arguments:
                du('CL')[:] = (dCLdCL1 + dCLdCL2 + dCLdCL3) * \
                    df('CL')
            if self.get_id('Wf') in arguments:
                dWf[:] = (dCLdWf1 + dCLdWf2 + dCLdWf3) * \
                    df('CL')
            if self.get_id('gamma') in arguments:
                dgamma[:] = (dCLdGamma1 + dCLdGamma2 +
                             dCLdGamma3) * df('CL')
            if self.get_id('CT') in arguments:
                dCT[:] = (dCLdThrust1 + dCLdThrust2 + 
                          dCLdThrust3) * df('CL')
            if self.get_id('alpha') in arguments:
                dalpha[:] = (dCLdAlpha1 + dCLdAlpha2 +
                             dCLdAlpha3) * df('CL')

class Sys_alpha(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        a_IC = self.kwargs['a_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('alpha', size=numPts, val=a_IC)
        self._declare_argument('CL', indices=iPts)
        self._declare_argument('eta', indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']

        CL = p('CL')
        eta = p('eta') * 1e-1
        alpha = mission.get_alpha(self.numElem, CL, eta)
        #print "current alpha: ", alpha
        u('alpha')[:] = alpha / 1e-1

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']

        CL = p('CL')
        eta = p('eta') * 1e-1

        dCL = dp('CL')
        deta = dp('eta')
        
        [dAdEta, dAdCL] = mission.get_alpha_d(self.numElem, CL, eta)
        
        if self.mode == 'fwd':
            dg('alpha')[:] = 0.0
            if self.get_id('eta') in arguments:
                dg('alpha')[:] += (dAdEta * deta)
            if self.get_id('CL') in arguments:
                dg('alpha')[:] += dAdCL * dCL / 1e-1
        if self.mode == 'rev':
            dCL[:] = 0.0
            deta[:] = 0.0
            du('alpha')[:] = 0.0
            if self.get_id('eta') in arguments:
                deta[:] = dAdEta * dg('alpha')
            if self.get_id('CL') in arguments:
                dCL[:] = dAdCL * dg('alpha')
            if self.get_id('alpha') in arguments:
                du('alpha')[:] = 0.0

class Sys_CD(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CD_IC = self.kwargs['CD_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('CD', size=numPts, val=CD_IC)
        self._declare_argument('CL', indices=iPts)
        self._declare_argument(['AR',0], indices=[0])
        self._declare_argument(['e',0], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        CL = p('CL')
        AR = p(['AR',0])
        e = p(['e',0])
        CD = mission.get_cd(self.numElem, AR, e, CL)
        #print "current CD: ", CD
        u('CD')[:] = CD / 1e-1

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']

        AR = p(['AR',0])
        e = p(['e',0])
        CL = p('CL')

        dAR = dp(['AR', 0])
        de = dp(['e', 0])
        dCL = dp('CL')
        
        [dCDdAR, dCDde, dCDdCL] = mission.get_cd_d(self.numElem, AR, e, CL)
        
        if self.mode == 'fwd':
            dg('CD')[:] = 0.0
            if self.get_id('AR') in arguments:
                dg('CD')[:] += dCDdAR * dAR / 1e-1
            if self.get_id('e') in arguments:
                dg('CD')[:] += dCDde * de / 1e-1
            if self.get_id('CL') in arguments:
                dg('CD')[:] += dCDdCL * dCL / 1e-1

        if self.mode == 'rev':
            dAR[:] = 0.0
            de[:] = 0.0
            dCL[:] = 0.0
            du('CD')[:] = 0.0
            if self.get_id('AR') in arguments:
                dAR[:] = dCDdAR * dg('CD')
            if self.get_id('e') in arguments:
                de[:] = dCDde * dg('CD')
            if self.get_id('CL') in arguments:
                dCL[:] = dCDdCL * dg('CD')
            if self.get_id('CD') in arguments:
                du('CD')[:] = 0.0

class Sys_Wf(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        Wf_IC = self.kwargs['Wf_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        numSeg = self.kwargs['numSeg']

        self._declare_variable('Wf', size=numPts, val=Wf_IC)
        self._declare_argument('v', indices=iPts)
        self._declare_argument('gamma', indices=iPts)
        self._declare_argument('CT', indices=iPts)
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument('SFC', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        #self._declare_argument(['WfSeg',0], indices=[self.copy+1])
        self._declare_argument(['S',0], indices=[0])
        if self.copy != numSeg-1:
            self._declare_argument(['Wf', self.copy+1], indices=[0])
        
    def apply_F(self):
        self.numInt = self.kwargs['numInt']
        numSeg = self.kwargs['numSeg']
        
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        x = p(['x',0]) * 1e6
        g = p(['g',0])
        v = p('v') * 1e2
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        SFC = p('SFC') * 1e-6
        rho = p('rho')
        WfIn = u('Wf') * 1e6
        #WfSeg = p(['WfSeg',0]) * 1e6
        if self.copy != numSeg-1:
            WfSeg = p(['Wf', self.copy+1]) * 1e6
        else:
            WfSeg = 0.0
        S = p(['S',0]) * 1e2
        xInt = numpy.linspace(x[0], x[1], self.numElem+1)

        Wf = mission.get_wf(self.numInt, self.numElem, xInt, v, gamma,
                            CT, SFC, rho, WfIn, g, WfSeg, S)
        
        #print "current Wf: ", WfIn
        #print "current WfRes: ", Wf-WfIn
        f('Wf')[:] = (WfIn-Wf) / 1e6

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dFdpu(self, arguments):
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']
        numSeg = self.kwargs['numSeg']

        g = 9.81
        if self.copy != numSeg-1:
            WfSeg = p(['Wf', self.copy+1]) * 1e6
        else:
            WfSeg = 0.0
        
        S = p(['S',0]) * 1e2
        x_ends = p(['x',0]) * 1e6
        v = p('v') * 1e2
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        SFC = p('SFC') * 1e-6
        rho = p('rho')
        x = numpy.linspace(x_ends[0],x_ends[1],self.numElem+1)

        #dWfSeg = dp(['WfSeg', 0])
        dS = dp(['S', 0])
        dx = dp(['x', 0])
        dv = dp('v')
        dgamma = dp('gamma')
        dCT = dp('CT')
        dSFC = dp('SFC')
        drho = dp('rho')

        [dWfdCT1, dWfdCT2, dWfdSFC1, dWfdSFC2, dWfdV1, dWfdV2,
         dWfdGamma1, dWfdGamma2, dWfdRho1, dWfdRho2,
         dWfdS] = mission.get_wf_d(self.numElem, self.numInt,
                                   g, WfSeg, S, x, v, gamma,
                                   CT, SFC, rho)

        if self.mode == 'fwd':
            df('Wf')[:] = 0.0
            if self.get_id('CT') in arguments:
                df('Wf')[:] -= ((dWfdCT1 + dWfdCT2) * 
                                dCT) * 1e-1/1e6
            if self.get_id('SFC') in arguments:
                df('Wf')[:] -= ((dWfdSFC1 + dWfdSFC2) * 
                                dSFC) * 1e-6/1e6
            if self.get_id('v') in arguments:
                df('Wf')[:] -= ((dWfdV1 + dWfdV2) * 
                                dv) * 1e2/1e6
            if self.get_id('gamma') in arguments:
                df('Wf')[:] -= ((dWfdGamma1 + dWfdGamma2) * 
                                dgamma) * 1e-1/1e6
            if self.get_id('rho') in arguments:
                df('Wf')[:] -= ((dWfdRho1 + dWfdRho2) * 
                                drho) * 1/1e6
            if self.get_id('S') in arguments:
                df('Wf')[:] -= dWfdS * dS * 1e2/1e6
            #if self.get_id('WfSeg') in arguments:
                #df('Wf')[:] -= (numpy.ones(self.numElem+1) *
                                     #dWfSeg)
            if self.get_id('Wf') in arguments:
                df('Wf')[:] += du('Wf')
        if self.mode == 'rev':
            dCT[:] = 0.0
            dSFC[:] = 0.0
            dv[:] = 0.0
            dgamma[:] = 0.0
            dS[:] = 0.0
            #dWfSeg[:] = 0.0
            du('Wf')[:] = 0.0
            if self.get_id('CT') in arguments:
                dCT[:] = (dWfdCT1 + dWfdCT2) * df('Wf')
            if self.get_id('SFC') in arguments:
                dSFC[:] = (dWfdSFC1 + dWfdSFC2) * df('Wf')
            if self.get_id('v') in arguments:
                dv[:] = (dWfdV1 + dWfdV2) * df('Wf')
            if self.get_id('gamma') in arguments:
                dgamma[:] = (dWfdGamma1 + dWfdGamma2) * \
                    dp('gamma')
            if self.get_id('S') in arguments:
                dS[:] = dWfdS * df('Wf')
            #if self.get_id('WfSeg') in arguments:
                #dWfSeg[:] = numpy.ones(self.numElem+1) * \
                    #df('Wf')
            if self.get_id('Wf') in arguments:
                du('Wf')[:] = df('Wf')

class Sys_CM(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CM_IC = self.kwargs['CM_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('CM', size=numPts, val=CM_IC)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['chord',0], indices=[0])
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        self._declare_argument('v', indices=iPts)
        self._declare_argument('rho', indices=iPts)

    def apply_F(self):
        self.numInt = self.kwargs['numInt']
        
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        S = p(['S',0]) * 1e2
        chord = p(['chord',0])
        x_ends = p(['x',0]) * 1e6
        v = p('v') * 1e2
        rho = p('rho')
        CM = u('CM')
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        CMRes = mission.get_cm(self.numElem, self.numInt, S, chord,
                               x_int, v, rho, CM)
        #print "current CM: ", CM
        #print "current CMRes: ", CMRes

        f('CM')[:] = CMRes

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dFdpu(self, arguments):
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']

        S = p(['S',0]) * 1e2
        chord = p(['chord',0])
        x_ends = p(['x',0]) * 1e6
        v = p('v') * 1e2
        rho = p('rho')
        CM = u('CM')
        x = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        dS = dp(['S', 0])
        dchord = dp(['chord', 0])
        dx = dp(['x', 0])
        dv = dp('v')
        drho = dp('rho')
        
        [dCMdS, dCMdC, dCMdV1, dCMdV2, dCMdV3, dCMdRho1, dCMdRho2,
         dCMdRho3, dCMdCM1, dCMdCM2,
         dCMdCM3] = mission.get_cm_d(self.numElem, self.numInt, S,
                                     chord, x, v, rho, CM)

        if self.mode == 'fwd':
            df('CM')[:] = 0.0
            if self.get_id('S') in arguments:
                df('CM')[:] += (dCMdS * dS) * 1e2
            if self.get_id('chord') in arguments:
                df('CM')[:] += dCMdC * dchord
            if self.get_id('v') in arguments:
                df('CM')[:] += ((dCMdV1 + dCMdV2 + dCMdV3) *
                                dv) * 1e2
            if self.get_id('rho') in arguments:
                df('CM')[:] += ((dCMdRho1 + dCMdRho2 + dCMdRho3) *
                                drho)
            if self.get_id('CM') in arguments:
                df('CM')[:] += (dCMdCM1 + dCMdCM2 + dCMdCM3) * du('CM')
        if self.mode == 'rev':
            dS[:] = 0.0
            dchord[:] = 0.0
            dv[:] = 0.0
            drho[:] = 0.0
            du('CM')[:] = 0.0
            if self.get_id('S') in arguments:
                dS[:] = dCMdS * df('CM')
            if self.get_id('chord') in arguments:
                dchord[:] = dCMdC * df('CM')
            if self.get_id('v') in arguments:
                dv[:] = (dCMdV1 + dCMdV2 + dCMdV3) * \
                    df('CM')
            if self.get_id('rho') in arguments:
                drho[:] = (dCMdRho1 + dCMdRho2 + dCMdRho2) * \
                    df('CM')
            if self.get_id('CM') in arguments:
                du('CM')[:] = (dCMdCM1 + dCMdCM2 + dCMdCM3) * df('CM')

class Sys_eta(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        e_IC = self.kwargs['e_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('eta', size=numPts, val=e_IC)
        self._declare_argument('CM', indices=iPts)
        self._declare_argument('alpha', indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        CM = p('CM')
        alpha = p('alpha') * 1e-1
        eta = mission.get_eta(self.numElem, CM, alpha)
        #print "current eta: ", eta
        u('eta')[:] = eta / 1e-1

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']

        CM = p('CM')
        alpha = p('alpha') * 1e-1

        dCM = dp('CM')
        dalpha = dp('alpha')

        [dEdAlpha, dEdCM] = mission.get_eta_d(self.numElem, alpha, CM)

        if self.mode == 'fwd':
            dg('eta')[:] = 0.0
            if self.get_id('CM') in arguments:
                dg('eta')[:] += dEdCM * dCM / 1e-1
            if self.get_id('alpha') in arguments:
                dg('eta')[:] += dEdAlpha * dalpha

        if self.mode == 'rev':
            dCM[:] = 0.0
            dalpha[:] = 0.0
            du('eta')[:] = 0.0
            if self.get_id('CM') in arguments:
                dCM[:] = dEdCM * dg('eta')
            if self.get_id('alpha') in arguments:
                dalpha[:] = dEdAlpha * dg('eta')
            if self.get_id('eta') in arguments:
                du('eta')[:] = dg('eta')

class Sys_x(ImplicitSystem):
    def _declare(self):
        x_IC = self.kwargs['x_IC']
        self.numSeg = self.kwargs['numSeg']
        numElem = self.kwargs['numElem']

        self._declare_variable('x', size=self.numSeg+1, val=x_IC)
        for i in xrange(self.numSeg):
            self._declare_argument(['h_ends', i], indices=[0,1])
            self._declare_argument(['h', i], indices=[numElem[i]])
        
    def apply_F(self):
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        

'''
class Sys_WfSeg(ExplicitSystem):
    def _declare(self):
        WfSeg_IC = self.kwargs['WfSeg_IC']
        numSeg = self.kwargs['numSeg']

        self._declare_variable('WfSeg', size=numSeg+1, val=WfSeg_IC)
        for i in xrange(numSeg):
            self._declare_argument(['Wf', i], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']

        numSeg = self.kwargs['numSeg']

        Wf_ends = numpy.zeros(numSeg)
        for i in xrange(numSeg):
            Wf_ends[i] = p(['Wf', i]) * 1e6
        WfSeg = u(['WfSeg', 0]) * 1e6
        print 'Wf', Wf_ends

        diff = 0
        for i in xrange(numSeg):
            j = numSeg-i-1
            print 'WfSeg', Wf_ends[j], WfSeg[j]
            diff += Wf_ends[j] - WfSeg[j]
            u(['WfSeg', 0])[j] += diff/1e6

        print 'new WfSeg', u(['WfSeg', 0])*1e6

    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)
'''
class Trajectory(object):

    def __init__(self, numSeg):
        self.numSeg = numSeg
        self.hPts = numpy.zeros(numSeg+1)
        self.xPts = numpy.zeros(numSeg+1)
        self.vPts = numpy.zeros(numSeg+1)
        self.MPts = numpy.zeros(numSeg+1)
        self.vMType = []
        self.tPts = numpy.zeros(numSeg)
        self.hDotPts = numpy.zeros(numSeg)
        self.thType = []
        self.accType = []
        self.opt = []
        self.numElem = numpy.zeros(numSeg, int)
        self.crtPt = 0

    def add_seg_point(self, h, x, v=-1, M=-1, tau=-1, hDot=-1, numElem=-1):
        
        if v == -1 and M == -1:
            raise Exception('Either velocity or Mach number has to be specified at each point')
        if v != -1 and M != -1:
            raise Exception('Only specify velocity OR Mach number at each point')
        if tau == -1 and hDot == -1 and self.crtPt != 0:
            raise Exception('Specify either throttle or descent rate for subsequent points')
        if tau != -1 and hDot != -1:
            raise Exception('Only specify throttle OR descent rate for subsequent points')
        if numElem == -1 and self.crtPt != 0:
            raise Exception('Specify number of elements per segment')
        if self.crtPt > self.numSeg:
            raise Exception('Segment point exceed number of segments')

        self.hPts[self.crtPt] = h/1e3
        self.xPts[self.crtPt] = x/1e6
        if v == -1:
            self.MPts[self.crtPt] = M
            self.vPts[self.crtPt] = -1
            self.vMType.append('M')
        else:
            self.MPts[self.crtPt] = -1
            self.vPts[self.crtPt] = v/1e2
            self.vMType.append('v')
        
        if tau == -1:
            self.tPts[self.crtPt-1] = -1
            self.hDotPts[self.crtPt-1] = hDot
            self.thType.append('h')
        else:
            self.tPts[self.crtPt-1] = tau
            self.hDotPts[self.crtPt-1] = -1
            self.thType.append('t')

        self.opt.append(False)
        self.numElem[self.crtPt-1] = numElem
        self.crtPt += 1

    def set_params(self, kw):
        self.S = kw['S']
        self.Wac = kw['Wac']
        self.cThrustSL = kw['cThrustSL']
        self.SFCSL = kw['SFCSL']
        self.chord = kw['chord']
        self.inertia = kw['inertia']
        self.AR = kw['AR']
        self.e = kw['e']
        self.g = kw['g']

    def set_IC(self, h_IC=None, a_IC=None, t_IC=None, e_IC=None,
                   CL_IC=None, CD_IC=None, CM_IC=None, Wf_IC=None, x_IC=None):
        # initial condition given at all element control points, total
        # of sum(numElem[i]) + 1 initial values for a given mission
        
        totalElem = 0
        for pt in xrange(self.crtPt-1):
            totalElem += self.numElem[pt]

        if h_IC == None:
            self.h_IC = []
            for pt in xrange(self.crtPt-1):
                if self.hDotPts[pt] == -1:
                    self.h_IC.extend(numpy.linspace(self.hPts[pt], self.hPts[pt+1], 
                                                    self.numElem[pt]+1))
                    #self.h_IC.extend(numpy.ones(self.numElem[pt]+1)*self.hPts[pt])
                    self.h_IC.pop()
                else:
                    h_ends = [self.hPts[pt],self.hPts[pt+1]]
                    v_ends = [self.vPts[pt],self.vPts[pt+1]]
                    if v_ends[pt] == -1:
                        v_ends[pt] = self.MPts[pt]*numpy.sqrt(1.4*288*(288.16-(6.5e-3)*self.hPts[pt]))
                    if v_ends[pt+1] == -1:
                        v_ends[pt+1] = self.MPts[pt+1]*numpy.sqrt(1.4*288*(288.16-(6.5e-3)*self.hPts[pt+1]))
                    h_dot = self.hDotPts[pt]
                    x_ends = [self.xPts[pt],self.xPts[pt+1]]
                    self.h_IC.extend(mission.get_h_des(self.numElem[pt],self.numInt, h_dot, x_ends, h_ends, v_ends))
                    self.h_IC.pop()
            #self.h_IC = 5000.0
            self.h_IC.extend(numpy.ones(1)*self.hPts[pt+1])
            #print "IC stage: ", self.h_IC
        else:
            self.h_IC = h_IC/1e3

        if a_IC == None:
            self.a_IC = numpy.ones(totalElem+1)*3.0*numpy.pi/180.0/1e-1
            #self.a_IC = 3.0*numpy.pi/180.0
        else:
            self.a_IC = a_IC/1e-1

        if e_IC == None:
            self.e_IC = numpy.zeros(totalElem+1)
            #self.e_IC = 0.0
        else:
            self.e_IC = e_IC

        if CL_IC == None:
            self.CL_IC = numpy.ones(totalElem+1)
            #self.CL_IC = 1.0
        else:
            self.CL_IC = CL_IC

        if CD_IC == None:
            self.CD_IC = numpy.ones(totalElem+1)*0.01/1e-1
            #self.CD_IC = 0.01
        else:
            self.CD_IC = CD_IC/1e-1

        if CM_IC == None:
            self.CM_IC = numpy.zeros(totalElem+1)
            #self.CM_IC = 0.0
        else:
            self.CM_IC = CM_IC

        if Wf_IC == None:
            self.Wf_IC = numpy.linspace(100000.0,self.Wf_final,totalElem+1)/1e6
            #self.Wf_IC = 100000.0
        else:
            self.Wf_IC = Wf_IC/1e6

        if x_IC == None:
            self.x_IC = numpy.linspace(0.0,self.range,totalElem+1)/1e6
            #self.x_IC = 0.0
        else:
            self.x_IC = x_IC/1e6
        
        h = numpy.array(self.h_IC) * 1e3
        SFC = numpy.array(self.SFCSL) * 1e-6
        x = numpy.array(self.x_IC) * 1e6
        
        self.SFC_IC = mission.get_sfc(totalElem, SFC, h)/1e-6
        # potential bug: must allow diff x spacing in diff segments
        self.gamma_IC = mission.get_gamma(totalElem, h, x)/1e-1
        self.Temp_IC = mission.get_temp(totalElem, h)/1e2
        temp = numpy.array(self.Temp_IC) * 1e2
        self.rho_IC = mission.get_rho(totalElem, self.g, temp)
        
        v_ends = numpy.zeros(self.numSeg+1)
        self.v_IC = []
        for seg in xrange(self.numSeg+1):
            if self.vPts[seg] == -1:
                v_ends[seg] = self.MPts[seg]*numpy.sqrt(1.4*288*self.Temp_IC[seg]*1e2)
            else:
                v_ends[seg] = self.vPts[seg]*1e2

        for seg in xrange(self.numSeg):
            self.v_IC.extend(mission.get_v(self.numElem[seg], [v_ends[seg], v_ends[seg+1]])/1e2)
            self.v_IC.pop()
        self.v_IC.append(v_ends[self.numSeg]/1e2)


        self.t_IC = []
        for seg in xrange(self.numSeg):
            if self.tPts[seg] != -1:
                if t_IC == None:
                    self.t_IC.extend(numpy.ones(self.numElem[seg])*1.0)
                #self.t_IC = 0.5
                else:
                    self.t_IC.extend(t_IC)
                    self.t_IC.pop()
            else:
                self.t_IC.extend(mission.get_ct_init(self.numElem[seg],
                                                     self.tPts[seg],
                                                     self.S,
                                                     self.cThrustSL,
                                                     self.rho_IC,
                                                     self.v_IC))
                self.t_IC.pop()
        self.t_IC.extend(numpy.ones(1)*1.0)

    def set_range(self, mission_range):
        self.range = mission_range

    def set_final_Wf(self, Wf):
        self.Wf_final = Wf

    def set_ingn_intl(self, numInt):
        self.numInt = numInt

    def initialize(self):
        ones = numpy.ones
        ls = numpy.linspace
        ne = self.numElem
        pt1 = 0
        pt2 = 0
        self.segments = []
        for pt in xrange(self.numSeg):
            pt2 += int(self.numElem[pt])
            h_IC = self.h_IC[pt1:pt2+1]
            t_IC = self.t_IC[pt1:pt2+1]
            e_IC = self.e_IC[pt1:pt2+1]
            a_IC = self.a_IC[pt1:pt2+1]
            CL_IC = self.CL_IC[pt1:pt2+1]
            CD_IC = self.CD_IC[pt1:pt2+1]
            CM_IC = self.CM_IC[pt1:pt2+1]
            Wf_IC = self.Wf_IC[pt1:pt2+1]
            SFC_IC = self.SFC_IC[pt1:pt2+1]
            gamma_IC = self.gamma_IC[pt1:pt2+1]
            Temp_IC = self.Temp_IC[pt1:pt2+1]
            rho_IC = self.rho_IC[pt1:pt2+1]
            v_IC = self.v_IC[pt1:pt2+1]
            pt1 += int(self.numElem[pt])
            
            self.segments.append(SerialSystem('segment', pt, 
                                              NL="NLN_GS", 
                                              LN="LIN_GS",
                                              LN_ilimit=1, 
                                              NL_ilimit=1, 
                                              output=True,
                                              subsystems=[
                        SerialSystem('seg_param',pt,
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     LN_ilimit=1, 
                                     NL_ilimit=1, 
                                     output=True,
                                     subsystems=[
                                IndVar('h_ends',pt,val=[self.hPts[pt],self.hPts[pt+1]],size=2),
                                IndVar('v_ends',pt,val=[self.vPts[pt],self.vPts[pt+1]],size=2),
                                IndVar('M_ends',pt,val=[self.MPts[pt],self.MPts[pt+1]],size=2),
                                IndVar('h_dot',pt,val=self.hDotPts[pt],size=1),
                                IndVar('tau_init',pt,val=self.tPts[pt],size=1),
                                ]),
                        SerialSystem('seg_analysis', pt,
                                     NL='NLN_GS', 
                                     LN='KSP_PC',
                                     PC='LIN_GS',
                                     LN_ilimit=100, 
                                     LN_rtol=1e-12, 
                                     LN_atol=1e-14, 
                                     NL_ilimit=100, 
                                     NL_rtol=1e-13,
                                     NL_atol=1e-13,
                                     PC_ilimit=20,
                                     PC_rtol=1e-1, 
                                     PC_atol=1e-4, 
                                     output=True,
                                     subsystems=[
                                Sys_CL('CL',pt,numElem=self.numElem[pt],numInt=self.numInt,CL_IC=CL_IC),
                                Sys_alpha('alpha',pt,numElem=self.numElem[pt],a_IC=a_IC),
                                Sys_CD('CD',pt,numElem=self.numElem[pt],CD_IC=CD_IC),
                                Sys_h('h',pt,climb=self.tPts[pt],h_IC=h_IC,numElem=self.numElem[pt],numInt=self.numInt),
                                Sys_SFC('SFC',pt,numElem=self.numElem[pt],SFC_IC=SFC_IC),
                                Sys_gamma('gamma',pt,numElem=self.numElem[pt],gamma_IC=gamma_IC),
                                Sys_Temp('Temp',pt,numElem=self.numElem[pt],Temp_IC=Temp_IC),
                                Sys_rho('rho',pt,numElem=self.numElem[pt],rho_IC=rho_IC),
                                Sys_v('v',pt,numElem=self.numElem[pt],v_IC=v_IC),
                                Sys_CT('CT',pt,numElem=self.numElem[pt],climb=self.tPts[pt],t_IC=t_IC),
                                Sys_tau('tau',pt,numElem=self.numElem[pt]),
                                Sys_Wf('Wf',pt,numElem=self.numElem[pt],numInt=self.numInt,Wf_IC=Wf_IC, numSeg=self.numSeg),
                                Sys_CM('CM',pt,numElem=self.numElem[pt],numInt=self.numInt,CM_IC=CM_IC),
                                Sys_eta('eta',pt,numElem=self.numElem[pt],e_IC=e_IC),
                                #]),
                                ])
                        ]))
                        
        self.segments.append(IndVar('x',val=self.xPts))
        self.segments.reverse()
        #self.segments.append(Sys_WfSeg('WfSeg', 0, numSeg=self.numSeg, WfSeg_IC=numpy.linspace(100000.0,0.0,self.numSeg+1)/1e6))
        #self.segments.append(IndVar('WfSeg',val=numpy.linspace(100000.0,0.0,self.numSeg+1)/1e6))
        
        self.mainMission = SerialSystem('mission',
                                        NL='NLN_GS', 
                                        LN='LIN_GS',
                                        LN_ilimit=1, 
                                        NL_ilimit=1, 
                                        output=True,
                                        subsystems=[
                SerialSystem('mission_param',
                             NL='NLN_GS',
                             LN='LIN_GS',
                             LN_ilimit=1, 
                             NL_ilimit=1, 
                             output=True,
                             subsystems=[
                        IndVar('S',val=self.S,size=1),
                        IndVar('Wac',val=self.Wac,size=1),
                        IndVar('cThrustSL',val=self.cThrustSL,size=1),
                        IndVar('SFCSL',val=self.SFCSL,size=1),
                        IndVar('chord',val=self.chord,size=1),
                        IndVar('inertia',val=self.inertia,size=1),
                        IndVar('AR',val=self.AR,size=1),
                        IndVar('e',val=self.e,size=1),
                        IndVar('g',val=self.g,size=1),
                        ]),
                SerialSystem('mission_analysis',
                             NL='NLN_GS', 
                             LN='LIN_GS',
                             LN_ilimit=1, 
                             NL_ilimit=1, 
                             output=True,
                             subsystems=self.segments,
                             #self.segments.reverse(),
                             #Sys_xRes(x_IC=self.x_IC),
                             )#NL='NLN_JC')
                ]).setup()
        '''
        self.mainMission = SerialSystem('mission',subsystems=[
                SerialSystem('global_ind',subsystems=[
                        IndVar('S',val=self.S,size=1), 
                        IndVar('Wac',val=self.Wac,size=1),
                        IndVar('cThrustSL',val=self.cThrustSL,size=1),
                        IndVar('SFCSL',val=self.SFCSL,size=1),
                        IndVar('chord',val=self.chord,size=1),
                        IndVar('inertia',val=self.inertia,size=1),
                        IndVar('AR',val=self.AR,size=1),
                        IndVar('e',val=self.e,size=1),
                        IndVar('g',val=self.g,size=1),
                        IndVar('x',val=self.xPts), 
                        IndVar('h_ends',pt,
                               val=[self.hPts[pt],self.hPts[pt+1]],
                               size=2),
                        IndVar('v_ends',pt,
                               val=[self.vPts[pt],self.vPts[pt+1]],
                               size=2),
                        IndVar('M_ends',pt,
                               val=[self.MPts[pt],self.MPts[pt+1]],size=2),
                        IndVar('h_dot',pt,val=self.hDotPts[pt],size=1),
                        IndVar('tau_init',pt,val=self.tPts[pt],size=1),
                        IndVar('numElem',pt,val=ne[pt],size=1),
                        ], NL="NLN_GS"),
                Sys_CL('CL',pt,numElem=self.numElem,numInt=self.numInt,CL_IC=self.CL_IC),
                Sys_alpha('alpha',pt,numElem=self.numElem,a_IC=self.a_IC),
                Sys_CD('CD',pt,numElem=self.numElem,CD_IC=self.CD_IC),
                Sys_h('h',pt,climb=self.tPts[pt],h_IC=self.h_IC,numElem=self.numElem[pt],numInt=self.numInt),
                Sys_SFC('SFC', pt,numElem=self.numElem[pt],SFC_IC=self.SFC_IC),
                Sys_gamma('gamma',pt,numElem=self.numElem[pt],gamma_IC=self.gamma_IC),
                Sys_Temp('Temp',pt,numElem=self.numElem[pt],Temp_IC=self.Temp_IC),
                Sys_rho('rho',pt,numElem=self.numElem[pt],rho_IC=self.rho_IC,),
                Sys_v('v',pt,numElem=self.numElem[pt],v_IC=self.v_IC),
                Sys_CT('CT',pt,numElem=self.numElem[pt],climb=self.tPts[pt],t_IC=t_IC),
                Sys_tau('tau',pt,numElem=self.numElem[pt]),
                Sys_Wf('Wf',pt,numElem=self.numElem[pt],numInt=self.numInt,Wf_IC=Wf_IC),
                Sys_CM('CM',pt,numElem=self.numElem[pt],numInt=self.numInt,CM_IC=self.CM_IC),
                Sys_eta('eta',pt,numElem=self.numElem[pt],e_IC=self.e_IC),
                IndVar('WfSeg',val=numpy.linspace(100000.0,0.0,len(self.segments)))
                ], NL='NEWTON', PC='LIN_GS', LN_ilimit=100, LN_rtol=1e-6, NL_ilimit=100, NL_rtol=1e-16, output=True).setup()
        '''        

        return self.mainMission

params = {
    'S': 4.278, #427.8,                  e2
    'Wac': 1.82419893, #185953*9.81,     e6
    'cThrustSL': 1.02, #1020000.0,       e6
    'SFCSL': 8.951, #8.951e-6,           e-6
    'chord': 8.15,
    'inertia': 4.1, #4.1e7,              e7
    'AR': 8.68,
    'e': 0.8,
    'g': 9.81,
    }

n = 1
missionProblem = Trajectory(1)
missionProblem.add_seg_point(0.0,0.0,v=150.0)
#for i in xrange(n):
missionProblem.add_seg_point(8000.0,20.0e3,v=270.0,tau=1.0,numElem=10)
#missionProblem.add_seg_point(11000.0,40.0e3,M=0.82,tau=1.0,numElem=10)
#missionProblem.add_seg_point(15000.0,60.0e3,M=0.82,tau=1.0,numElem=10)
#missionProblem.add_seg_point(8000.0,4500.0,M=0.82,hDot=-15.3,numElem=5)
#missionProblem.add_seg_point(0.0,5000.0,v=150.0,hDot=-15.3,numElem=5)
missionProblem.set_ingn_intl(100)
missionProblem.set_params(params)
missionProblem.set_range(1000.0e3)
missionProblem.set_final_Wf(0.0)
missionProblem.set_IC()

problemPtr = missionProblem.initialize()
#print 'h_IC: ', problemPtr.vec['u']['h',0]
#print 'SFC_IC: ', problemPtr.vec['u']['SFC',0]
#print 'gamma_IC: ', problemPtr.vec['u']['gamma',0]
#print 'Temp_IC: ', problemPtr.vec['u']['Temp',0]
#print 'rho_IC: ', problemPtr.vec['u']['rho',0]

problemPtr.compute(True).array
print
print 'Computing derivatives'
#print problemPtr.compute_derivatives('fwd', 'v_ends', output=True)#.array
#exit()
'''
problemPtr('segment').kwargs['NL'] = 'NLN_GS'
problemPtr.compute(True).array
'''
problemPtr.check_derivatives_all(fwd=True, rev=False)
"""
#A = numpy.array(problemPtr.compute_derivatives('fwd', 'S', output=True).array)
print problemPtr.compute_derivatives('fwd', 'S', output=True)
h = 1e-3
f0 = numpy.array(problemPtr.vec['u'].array)
problemPtr('S').value += h
problemPtr.vec['du'].array[:] = 1.0
problemPtr.compute(False)
f = numpy.array(problemPtr.vec['u'].array)
problemPtr.vec['du'].array[:] = (f-f0)/h
#print '--------------------NOW COMPUTING DERIVATIVES'
#B = numpy.array(problemPtr.vec['du'].array)
print problemPtr.vec['du']
#print numpy.vstack([A,B]).T
"""
exit()

'''
analytic_array = problemPtr.compute_derivatives('fwd','SFCSL').array[:]
problemPtr.vec['u'](['SFCSL',0])[:] += 1e-5
new_array = problemPtr.compute(True).array
diff_array = (new_array-old_array)/1e-10
print analytic_array
print diff_array
print analytic_array-diff_array
'''
#problemPtr.vec['u'].array[:] *= problemPtr.vec['u0'].array[:]
#print problemPtr.vec['u']
#print problemPtr.vec['u0']
'''
print "---------------------------------------"
for name in ['SFC', 'gamma', 'Temp', 'h', 'rho', 'CT', 'v', 'tau', 'CL',
             'alpha', 'CD', 'Wf', 'CM', 'eta']:
    result = problemPtr([name,0]).check_derivatives(problemPtr.variables.keys())
    print 'check ' + name
    print result[0]
    print result[1]
'''

totalElem = 0
total_x = numpy.array([])
total_h = numpy.array([])
total_v = numpy.array([])
total_a = numpy.array([])
total_t = numpy.array([])
total_e = numpy.array([])
total_w = numpy.array([])
total_p = numpy.array([])
for i in xrange(missionProblem.numSeg):
    totalElem += missionProblem.numElem[i]
    x_ends = problemPtr.vec['u'](['x', 0]) * 1e6
    x_int = numpy.linspace(x_ends[i], x_ends[i+1], missionProblem.numElem[i]+1)
    total_x = numpy.append(total_x, x_int)
    total_h = numpy.append(total_h, problemPtr.vec['u'](['h', i]))
    total_v = numpy.append(total_v, problemPtr.vec['u'](['v', i]))
    total_a = numpy.append(total_a, problemPtr.vec['u'](['alpha', i]))
    total_t = numpy.append(total_t, problemPtr.vec['u'](['tau', i]))
    total_e = numpy.append(total_e, problemPtr.vec['u'](['eta', i]))
    total_w = numpy.append(total_w, problemPtr.vec['u'](['Wf', i]))
    total_p = numpy.append(total_p, problemPtr.vec['u'](['rho', i]))
    numpy.delete(total_x, totalElem)
    numpy.delete(total_h, totalElem)
    numpy.delete(total_v, totalElem)
    numpy.delete(total_a, totalElem)
    numpy.delete(total_t, totalElem)
    numpy.delete(total_e, totalElem)
    numpy.delete(total_w, totalElem)
    numpy.delete(total_p, totalElem)
numSeg = missionProblem.numSeg
numpy.append(total_x, numpy.array([x_ends[numSeg]]))
numpy.append(total_h, numpy.array([missionProblem.hPts[numSeg]]))
numpy.append(total_v, numpy.array([problemPtr.vec['u'](['v', numSeg-1])[-1]]))
numpy.append(total_a, numpy.array([problemPtr.vec['u'](['alpha', numSeg-1])[-1]]))
numpy.append(total_t, numpy.array([problemPtr.vec['u'](['tau', numSeg-1])[-1]]))
numpy.append(total_e, numpy.array([problemPtr.vec['u'](['eta', numSeg-1])[-1]]))
numpy.append(total_w, numpy.array([problemPtr.vec['u'](['Wf', numSeg-1])[-1]]))
numpy.append(total_p, numpy.array([problemPtr.vec['u'](['rho', numSeg-1])[-1]]))
    
pylab.figure
pylab.subplot(711)
pylab.plot(total_x/1000.0, total_h)
pylab.ylabel('Altitude (km)')
pylab.subplot(712)
pylab.plot(total_x/1000.0, total_v*1e2)
pylab.ylabel('Velocity (m/s)')
pylab.subplot(713)
pylab.plot(total_x/1000.0, total_a*180.0e-1/numpy.pi)
pylab.ylabel('AoA (deg)')
pylab.subplot(714)
pylab.plot(total_x/1000.0, total_t)
pylab.ylabel('Throttle')
pylab.subplot(715)
pylab.plot(total_x/1000.0, total_e*180.0e-1/numpy.pi)
pylab.ylabel('Trim Angle (deg)')
pylab.subplot(716)
pylab.plot(total_x/1000.0, total_w*1e6/(9.81*0.804))
pylab.ylabel('Fuel (L)')
pylab.subplot(717)
pylab.plot(total_x/1000.0, total_p)
pylab.ylabel('rho')
pylab.xlabel('Distance (km)')
pylab.show()

'''
name = 'v'
result = problemPtr(['CL',0]).check_derivatives()
print 'checking CL-' + name
print result[0]
print result[1]
'''
"""
print problemPtr.vec['u']['h',0]
print problemPtr.vec['p']['SFC',0]['h',0]

main = problemPtr
if main(['h2',0]).comm is not None:
    print 'h2-derv:', main(['h2',0]).check_derivatives('fwd')
if main(['SFC',0]).comm is not None:
    print 'SFC-derv:'
    print main(['SFC',0]).check_derivatives('fwd')
"""
