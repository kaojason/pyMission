from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
import numpy
import copy
from framework import *
from optimization import *
import mission
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab
import MBI, scipy.sparse


class Sys_v_bspline(ExplicitSystem):

    def _declare(self):
        self.numPts = self.kwargs['numPts']
        self.numCP = self.kwargs['numCP']
        self.jac = self.kwargs['jac']

        self._declare_variable('v', size=self.numPts)
        self._declare_argument('v_CP', indices=range(self.numCP))

    def apply_G(self):
        u = self.vec['u']('v')
        p = self.vec['p']('v_CP') * 1e2
        u[:] = self.jac.dot(p[:]) / 1e2

    def apply_dGdp(self, args):
        dg = self.vec['dg']('v')
        dp = self.vec['dp']('v_CP')

        if self.mode == 'fwd':
            dg[:] = 0.0
            if self.get_id('v_CP') in args:
                dg[:] += self.jac.dot(dp[:])
        if self.mode == 'rev':
            dp[:] = 0.0
            if self.get_id('v_CP') in args:
                dp[:] += self.jac.T.dot(dg[:])


class Sys_h_bspline(ExplicitSystem):

    def _declare(self):
        self.numPts = self.kwargs['numPts']
        self.numCP = self.kwargs['numCP']
        self.jac = self.kwargs['jac']

        self._declare_variable('h', size=self.numPts)
        self._declare_argument('h_CP', indices=range(self.numCP))

    def apply_G(self):
        u = self.vec['u']('h')
        p = self.vec['p']('h_CP')
        u[:] = self.jac.dot(p[:])

    def apply_dGdp(self, args):
        dg = self.vec['dg']('h')
        dp = self.vec['dp']('h_CP')

        if self.mode == 'fwd':
            dg[:] = 0.0
            if self.get_id('h_CP') in args:
                dg[:] += self.jac.dot(dp[:])
        if self.mode == 'rev':
            dp[:] = 0.0
            if self.get_id('h_CP') in args:
                dp[:] += self.jac.T.dot(dg[:])


class Sys_x_bspline(ExplicitSystem):

    def _declare(self):
        self.numPts = self.kwargs['numPts']
        self.numCP = self.kwargs['numCP']
        self.jac = self.kwargs['jac']

        self._declare_variable('x', size=self.numPts)
        self._declare_argument('x_CP', indices=range(self.numCP))

    def apply_G(self):
        u = self.vec['u']('x')
        p = self.vec['p']('x_CP')
        u[:] = self.jac.dot(p[:])

    def apply_dGdp(self, args):
        dg = self.vec['dg']('x')
        dp = self.vec['dp']('x_CP')

        if self.mode == 'fwd':
            dg[:] = 0.0
            if self.get_id('x_CP') in args:
                dg[:] += self.jac.dot(dp[:])
        if self.mode == 'rev':
            dp[:] = 0.0
            if self.get_id('x_CP') in args:
                dp[:] += self.jac.T.dot(dg[:])


class Sys_gamma_bspline(ExplicitSystem):

    def _declare(self):
        self.numPts = self.kwargs['numPts']
        self.numCP = self.kwargs['numCP']
        self.jac = self.kwargs['jac']

        self._declare_variable('gamma', size=self.numPts)
        self._declare_argument('h_CP', indices=range(self.numCP))

    def apply_G(self):
        u = self.vec['u']('gamma')
        p = self.vec['p']('h_CP')
        u[:] = self.jac.dot(p[:])*1e6/1e3 / 1e-1

    def apply_dGdp(self, args):
        dg = self.vec['dg']('gamma')
        dp = self.vec['dp']('h_CP')

        if self.mode == 'fwd':
            dg[:] = 0.0
            if self.get_id('h_CP') in args:
                dg[:] += self.jac.dot(dp[:]) * 1e3/1e-1
        if self.mode == 'rev':
            dp[:] = 0.0
            if self.get_id('h_CP') in args:
                dp[:] += self.jac.T.dot(dg[:]) * 1e3/1e-1


class Sys_h(ImplicitSystem):
    def _declare(self):
        h_IC = self.kwargs['h_IC']
        self.climb = self.kwargs['climb']
        self.numPts = self.kwargs['numElem']+1
        self.numInt = self.kwargs['numInt']
        lower = 0.0
        iPts = range(self.numPts)

        self._declare_variable('h', size=self.numPts, 
                               val=h_IC, lower=lower)#, 
                               #u_scal=1e3, f_scal=1e3)
        self._declare_argument('h_ends', indices=[0,1])
        #self._declare_argument(['v_ends',-1], indices=[0,1])
        self._declare_argument('x', indices=iPts)
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
        self._nln_init()
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
            x_ends = p('x') * 1e6
            Wf = p('Wf') * 1e6
            CT = p('CT') * 1e-1
            alpha = p('alpha') * 1e-1
            CD = p('CD') * 1e-1
            rho = p('rho')
            v = p('v') * 1e2
            S = p(['S',0]) * 1e2
            h = u('h') * 1e3
            Wac = p(['Wac',0]) * 1e6
            
            h_new = numpy.zeros(self.numPts)

            h_new = mission.get_h(self.numPts-1, self.numInt, S, Wac, x_ends, 
                                  h_ends, Wf, CT, alpha, CD, rho, v)

            #print "h_new", h_new
            #print "h", h

            #for i in xrange(self.numPts):
            #    if h_new[i] < 0:
            #        h_new[i] = 0
                    #raise Exception('negative altitude!')
            
            #print "current h: ", h
            f('h')[:] = (h - h_new) / 1e3
        self._nln_final()
    
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
            x_ends = p('x') * 1e6
            h_ends = p('h_ends') * 1e3
            Wf = p('Wf') * 1e6
            CT = p('CT') * 1e-1
            alpha = p('alpha') * 1e-1
            CD = p('CD') * 1e-1
            rho = p('rho')
            v = p('v') * 1e2
            Wac = p(['Wac',0]) * 1e6

            dS = dp(['S',0])
            dx_ends = dp('x')
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
             
            dSTemp = numpy.zeros(numElem+1)
            dWacTemp = numpy.zeros(numElem+1)
            dxTemp0 = numpy.zeros(numElem+1)
            dxTemp1 = numpy.zeros(numElem+1)
            dWfTemp = numpy.zeros(numElem+1)
            drhoTemp = numpy.zeros(numElem+1)
            dvTemp = numpy.zeros(numElem+1)
            dCTTemp = numpy.zeros(numElem+1)
            daTemp = numpy.zeros(numElem+1)
            dCDTemp = numpy.zeros(numElem+1)
            
            if self.mode == 'fwd':
                df('h')[:] = 0.0
                for i in xrange(numElem):
                    if self.get_id('S') in arguments:
                        dSTemp[i+1] += dSTemp[i]
                        dSTemp[i+1] -= dhdS[i+1] * dS
                        df('h')[i+1] += dSTemp[i+1] * 1e2/1e3
                    if self.get_id('Wac') in arguments:
                        dWacTemp[i+1] += dWacTemp[i]
                        dWacTemp[i+1] -= dhdWac[i+1] * dWac
                        df('h')[i+1] += dWacTemp[i+1] * 1e6/1e3
                    #FIX X DERIVATIVE
                    if self.get_id('x') in arguments:
                        dxTemp0[i+1] += dxTemp0[i]
                        dxTemp0[i+1] -= (dhdx0[i+1] * dx_ends[0] +
                                         dhdx1[i+1] * dx_ends[1])
                        df('h')[i+1] += dxTemp0[i+1] * 1e6/1e3
                    if self.get_id('Wf') in arguments:
                        dWfTemp[i+1] += dWfTemp[i]
                        dWfTemp[i+1] -= (dhdWf1[i+1] * dWf[i] + 
                                         dhdWf2[i+1] * dWf[i+1])
                        df('h')[i+1] += dWfTemp[i+1] * 1e6/1e3
                    if self.get_id('rho') in arguments:
                        drhoTemp[i+1] += drhoTemp[i]
                        drhoTemp[i+1] -= (dhdRho1[i+1] * drho[i] +
                                          dhdRho2[i+1] * drho[i+1])
                        df('h')[i+1] += drhoTemp[i+1] /1e3
                    if self.get_id('v') in arguments:
                        dvTemp[i+1] += dvTemp[i]
                        dvTemp[i+1] -= (dhdV1[i+1] * dv[i] +
                                        dhdV2[i+1] * dv[i+1])
                        df('h')[i+1] += dvTemp[i+1] * 1e2/1e3
                    if self.get_id('CT') in arguments:
                        dCTTemp[i+1] += dCTTemp[i]
                        dCTTemp[i+1] -= (dhdCT1[i+1] * dCT[i] +
                                         dhdCT2[i+1] * dCT[i+1])
                        df('h')[i+1] += dCTTemp[i+1] * 1e-1/1e3
                    if self.get_id('alpha') in arguments:
                        daTemp[i+1] += daTemp[i]
                        daTemp[i+1] -= (dhdA1[i+1] * dalpha[i] +
                                        dhdA2[i+1] * dalpha[i+1])
                        df('h')[i+1] += daTemp[i+1] * 1e-1/1e3
                    if self.get_id('CD') in arguments:
                        dCDTemp[i+1] += dCDTemp[i]
                        dCDTemp[i+1] -= (dhdCD1[i+1] * dCD[i] + 
                                         dhdCD2[i+1] * dCD[i+1])
                        df('h')[i+1] += dCDTemp[i+1] * 1e-1/1e3
                if self.get_id('h_ends') in arguments:
                    df('h')[:] -= dh_ends[0]
                if self.get_id('h') in arguments:
                    df('h')[:] += du('h')[:]

            if self.mode == 'rev':
                dx_ends[:] = 0.0
                dWf[:] = 0.0
                drho[:] = 0.0
                dv[:] = 0.0
                dalpha[:] = 0.0
                dCT[:] = 0.0
                dCD[:] = 0.0
                dS[:] = 0.0
                dWac[:] = 0.0
                du('h')[:] = 0.0

                dhSum = numpy.sum(df('h'))
                dhSum -= df('h')[0]
                
                for i in xrange(numElem):
                    if self.get_id('S') in arguments:
                        dSTemp[i+1] += dSTemp[i]
                        dSTemp[i+1] -= dhdS[i+1]
                        dS += dSTemp[i+1] * df('h')[i+1] * 1e2/1e3
                    if self.get_id('Wac') in arguments:
                        dWacTemp[i+1] += dWacTemp[i]
                        dWacTemp[i+1] -= dhdWac[i+1]
                        dWac += dWacTemp[i+1] * df('h')[i+1] * 1e6/1e3
                    #FIX X DERIVATIVE
                    if self.get_id('x') in arguments:
                        dxTemp0[i+1] += dxTemp0[i]
                        dxTemp1[i+1] += dxTemp1[i]
                        dxTemp0[i+1] -= dhdx0[i+1]
                        dxTemp1[i+1] -= dhdx1[i+1]
                        dx_ends[0] += dxTemp0[i+1] * df('h')[i+1] * 1e6/1e3
                        dx_ends[1] += dxTemp1[i+1] * df('h')[i+1] * 1e6/1e3
                    if self.get_id('Wf') in arguments:
                        dWf[i] -= dhdWf1[i+1] * dhSum * 1e6/1e3
                        dWf[i+1] -= dhdWf2[i+1] * dhSum * 1e6/1e3
                    if self.get_id('rho') in arguments:
                        drho[i] -= dhdRho1[i+1] * dhSum /1e3
                        drho[i+1] -= dhdRho2[i+1] * dhSum /1e3
                    if self.get_id('v') in arguments:
                        dv[i] -= dhdV1[i+1] * dhSum * 1e2/1e3
                        dv[i+1] -= dhdV2[i+1] * dhSum * 1e2/1e3
                    if self.get_id('CT') in arguments:
                        dCT[i] -= dhdCT1[i+1] * dhSum * 1e-1/1e3
                        dCT[i+1] -= dhdCT2[i+1] * dhSum * 1e-1/1e3
                    if self.get_id('alpha') in arguments:
                        dalpha[i] -= dhdA1[i+1] * dhSum * 1e-1/1e3
                        dalpha[i+1] -= dhdA2[i+1] * dhSum * 1e-1/1e3
                    if self.get_id('CD') in arguments:
                        dCD[i] -= dhdCD1[i+1] * dhSum * 1e-1/1e3
                        dCD[i+1] -= dhdCD2[i+1] * dhSum * 1e-1/1e3
                    dhSum -= df('h')[i+1]
                if self.get_id('h_ends') in arguments:
                    dh_ends[0] = -numpy.sum(df('h')) 
                    dh_ends[1] = 0.0
                if self.get_id('h') in arguments:
                    du('h')[:] += df('h')[:]

        self._lin_final()
            
class Sys_SFC(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        SFC_IC = self.kwargs['SFC_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('SFC', size=numPts, val=SFC_IC)#,
                               #u_scal=1e-6, f_scal=1e-6) 
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
                dg('SFC')[:] += dSFCSL

        if self.mode == 'rev':
            dh[:] = 0.0
            dSFCSL[:] = 0.0
            du('SFC')[:] = 0.0
            
            if self.get_id('h') in arguments:
                dh[:] += dSFC_dh * dg('SFC') * 1e3/1e-6
            if self.get_id('SFCSL') in arguments:
                dSFCSL[:] += numpy.sum(dg('SFC'))

class Sys_Temp(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        Temp_IC = self.kwargs['Temp_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('Temp', size=numPts, val=Temp_IC,
                               lower=0.001)#, u_scal=1e2, f_scal=1e2)
        self._declare_argument('h', indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p('h') * 1e3

        Temp = mission.get_temp(self.numElem, h)
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
                dh[:] = dTemp_dh * dg('Temp') * 1e3/1e2
    
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
        self._nln_init()
        p = self.vec['p']
        u = self.vec['u']
        g = p(['g',0])
        Temp = p('Temp') * 1e2

        rho = mission.get_rho(self.numElem, g, Temp)

        u('rho')[:] = rho
        self._nln_final()

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
                dTemp[:] = drho_dTemp * dg('rho') * 1e2
    
class Sys_v(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        v_IC = self.kwargs['v_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        ##print "v_IC: ", v_IC
        self._declare_variable('v', size=numPts, val=v_IC,
                               lower=0.0)#, u_scal=1e2, f_scal=1e2)
        self._declare_argument('Temp', indices=iPts)
        self._declare_argument('v_ends', indices=[0,1])
        self._declare_argument('M_ends', indices=[0,1])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        v_ends = p('v_ends') * 1e2
        M_ends = p('M_ends')
        Temp = p('Temp') * 1e2

        v_copy = copy.copy(v_ends)
        '''
        if v_copy[0] != -1:
            v_copy[0] *= 1e2
        if v_copy[1] != -1:
            v_copy[1] *= 1e2
        '''
        if v_copy[0] == -1:
            v_copy[0] = M_ends[0]*numpy.sqrt(1.4*288*Temp[0])
        if v_copy[1] == -1:
            v_copy[1] = M_ends[1]*numpy.sqrt(1.4*288*Temp[self.numElem])

        #v = mission.get_v(self.numElem, v_copy**2)
        v = mission.get_v(self.numElem, v_copy)
        #print "current v: ", v
        u('v')[:] = v / 1e2
        
        #print 'v-v_ends', v_ends
        #print 'Temp', Temp
        #print 'v-v', v

    
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
        v_ends = p('v_ends') * 1e2
        v_copy = copy.copy(v_ends)

        '''
        if v_copy[0] != -1:
            v_copy[0] *= 1e2
        if v_copy[1] != -1:
            v_copy[1] *= 1e2
        '''
        dM = dp('M_ends')
        dv = dp('v_ends')
        dT = dp('Temp')
        v = p('v_ends')

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
                dv[0] = numpy.sum(dv_dvEnds0 * dg('v'))
                dv[1] = numpy.sum(dv_dvEnds1 * dg('v'))
            if self.get_id('M_ends') in arguments:
                dM[0] = numpy.sum(dv_dMEnds0 * dg('v'))
                dM[1] = numpy.sum(dv_dMEnds1 * dg('v'))

class Sys_CT(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        t_IC = self.kwargs['t_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('CT', size=numPts, val=t_IC)#,
                               #u_scal=1e-1, f_scal=1e-1)
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

            #print 'CT-v', v

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
                    dcThrustSL[:] = numpy.sum(dCTdcThrustSL * dg('CT')) * 1e6/1e-1
                if self.get_id('h') in arguments:
                    dh[:] = dCTdH * dg('CT') * 1e3/1e-1
                if self.get_id('rho') in arguments:
                    drho[:] = dCTdRho * dg('CT') /1e-1
                if self.get_id('v') in arguments:
                    dv[:] = dCTdV * dg('CT') * 1e2/1e-1
                if self.get_id('S') in arguments:
                    dS[:] = numpy.sum(dCTdS * dg('CT')) * 1e2/1e-1
                if self.get_id('tau_init') in arguments:
                    dtau[:] = numpy.sum(dCTdTau * dg('CT')) /1e-1
                
class Sys_CT_opt(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        self.numInt = self.kwargs['numInt']
        CT_IC = self.kwargs['CT_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('CT', size=numPts, val=CT_IC)
        self._declare_argument('alpha', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        self._declare_argument('v', indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument('CD', indices=iPts)
        self._declare_argument(['Wac',0], indices=[0])
        self._declare_argument('Wf', indices=iPts)
        self._declare_argument('gamma', indices=iPts)
        self._declare_argument('x', indices=iPts)

    def apply_F(self):
        self._nln_init()
        numElem = self.numElem
        numInt = self.numInt
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']
        
        alpha = p('alpha') * 1e-1
        rho = p('rho')
        v = p('v') * 1e2
        S = p(['S',0]) * 1e2
        CD = p('CD') * 1e-1
        Wac = p(['Wac',0]) * 1e6
        Wf = p('Wf') * 1e6
        gamma = p('gamma') * 1e-1
        CT = u('CT') * 1e-1
        x = p('x') * 1e6

        CTRes = mission.get_ct(numElem, numInt, S, Wac, x, alpha,
                               rho, v, CD, Wf, gamma, CT)

        f('CT')[:] = CTRes / 1e7
        self._nln_final()

    def apply_dFdpu(self, arguments):
        self._lin_init()
        numElem = self.numElem
        numInt = self.numInt
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']
        
        alpha = p('alpha') * 1e-1
        rho = p('rho')
        v = p('v') * 1e2
        S = p('S') * 1e2
        CD = p('CD') * 1e-1
        Wac = p(['Wac',0]) * 1e6
        Wf = p('Wf') * 1e6
        gamma = p('gamma') * 1e-1
        CT = u('CT') * 1e-1
        x = p('x') * 1e6

        dalpha = dp('alpha')
        drho = dp('rho')
        dv = dp('v')
        dS = dp('S')
        dCD = dp('CD')
        dWac = dp('Wac')
        dWf = dp('Wf')
        dgamma = dp('gamma')
        dCT = du('CT')
        
        [dCTdS, dCTdWac, dCTdAlpha1, dCTdAlpha2, dCTdAlpha3, dCTdRho1,
         dCTdRho2, dCTdRho3, dCTdV1, dCTdV2, dCTdV3, dCTdCD1, dCTdCD2,
         dCTdCD3, dCTdWf1, dCTdWf2, dCTdWf3, dCTdGamma1,
         dCTdGamma2, dCTdGamma3, dCTdCT1, dCTdCT2,
         dCTdCT3] = mission.get_ct_d(numElem, numInt, S, Wac, x, alpha,
                                     rho, v, CD, Wf, gamma, CT)

        if self.mode == 'fwd':
            df('CT')[:] = 0.0
            for i in xrange(self.numElem+1):
                if i == 0:
                    if self.get_id('S') in arguments:
                        df('CT')[i] += dCTdS[i] * dS * 1e2/1e7
                    if self.get_id('Wac') in arguments:
                        df('CT')[i] += dCTdWac[i] * dWac * 1e6/1e7
                    if self.get_id('alpha') in arguments:
                        df('CT')[i] += (dCTdAlpha2[i] * dalpha[i] +
                                        dCTdAlpha3[i] * dalpha[i+1]) * 1e-1/1e7
                    if self.get_id('rho') in arguments:
                        df('CT')[i] += (dCTdRho2[i] * drho[i] + 
                                        dCTdRho3[i] * drho[i+1]) /1e7
                    if self.get_id('v') in arguments:
                        df('CT')[i] += (dCTdV2[i] * dv[i] + 
                                        dCTdV3[i] * dv[i+1]) * 1e2/1e7
                    if self.get_id('CD') in arguments:
                        df('CT')[i] += (dCTdCD2[i] * dCD[i] +
                                        dCTdCD3[i] * dCD[i+1]) * 1e-1/1e7
                    if self.get_id('Wf') in arguments:
                        df('CT')[i] += (dCTdWf2[i] * dWf[i] + 
                                        dCTdWf3[i] * dWf[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        df('CT')[i] += (dCTdGamma2[i] * dgamma[i] + 
                                        dCTdGamma3[i] * dgamma[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        df('CT')[i] += (dCTdCT2[i] * dCT[i] + 
                                        dCTdCT3[i] * dCT[i+1]) * 1e-1/1e7
                elif i == self.numElem:
                    if self.get_id('S') in arguments:
                        df('CT')[i] += dCTdS[i] * dS * 1e2/1e7
                    if self.get_id('Wac') in arguments:
                        df('CT')[i] += dCTdWac[i] * dWac * 1e6/1e7
                    if self.get_id('alpha') in arguments:
                        df('CT')[i] += (dCTdAlpha1[i] * dalpha[i-1] + 
                                        dCTdAlpha2[i] * dalpha[i]) * 1e-1/1e7
                    if self.get_id('rho') in arguments:
                        df('CT')[i] += (dCTdRho1[i] * drho[i-1] + 
                                        dCTdRho2[i] * drho[i]) /1e7
                    if self.get_id('v') in arguments:
                        df('CT')[i] += (dCTdV1[i] * dv[i-1] + 
                                        dCTdV2[i] * dv[i]) * 1e2/1e7
                    if self.get_id('CD') in arguments:
                        df('CT')[i] += (dCTdCD1[i] * dCD[i-1] + 
                                        dCTdCD2[i] * dCD[i]) * 1e-1/1e7
                    if self.get_id('Wf') in arguments:
                        df('CT')[i] += (dCTdWf1[i] * dWf[i-1] + 
                                        dCTdWf2[i] * dWf[i]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        df('CT')[i] += (dCTdGamma1[i] * dgamma[i-1] +
                                        dCTdGamma2[i] * dgamma[i]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        df('CT')[i] += (dCTdCT1[i] * dCT[i-1] +
                                        dCTdCT2[i] * dCT[i]) * 1e-1/1e7
                else:
                    if self.get_id('S') in arguments:
                        df('CT')[i] += dCTdS[i] * dS * 1e2/1e7
                    if self.get_id('Wac') in arguments:
                        df('CT')[i] += dCTdWac[i] * dWac * 1e6/1e7
                    if self.get_id('alpha') in arguments:
                        df('CT')[i] += (dCTdAlpha1[i] * dalpha[i-1] +
                                        dCTdAlpha2[i] * dalpha[i] + 
                                        dCTdAlpha3[i] * dalpha[i+1]) * 1e-1/1e7
                    if self.get_id('rho') in arguments:
                        df('CT')[i] += (dCTdRho1[i] * drho[i-1] +
                                        dCTdRho2[i] * drho[i] + 
                                        dCTdRho3[i] * drho[i+1]) / 1e7
                    if self.get_id('v') in arguments:
                        df('CT')[i] += (dCTdV1[i] * dv[i-1] + 
                                        dCTdV2[i] * dv[i] + 
                                        dCTdV3[i] * dv[i+1]) * 1e2/1e7
                    if self.get_id('CD') in arguments:
                        df('CT')[i] += (dCTdCD1[i] * dCD[i-1] + 
                                        dCTdCD2[i] * dCD[i] + 
                                        dCTdCD3[i] * dCD[i+1]) * 1e-1/1e7
                    if self.get_id('Wf') in arguments:
                        df('CT')[i] += (dCTdWf1[i] * dWf[i-1] + 
                                        dCTdWf2[i] * dWf[i] + 
                                        dCTdWf3[i] * dWf[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        df('CT')[i] += (dCTdGamma1[i] * dgamma[i-1] + 
                                        dCTdGamma2[i] * dgamma[i] + 
                                        dCTdGamma3[i] * dgamma[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        df('CT')[i] += (dCTdCT1[i] * dCT[i-1] + 
                                        dCTdCT2[i] * dCT[i] + 
                                        dCTdCT3[i] * dCT[i+1]) * 1e-1/1e7
        if self.mode == 'rev':
            dS[:] = 0.0
            dWac[:] = 0.0
            dalpha[:] = 0.0
            drho[:] = 0.0
            dv[:] = 0.0
            dCD[:] = 0.0
            dWf[:] = 0.0
            dgamma[:] = 0.0
            dCT[:] = 0.0

            for i in xrange(self.numElem+1):
                if i == 0:
                    if self.get_id('S') in arguments:
                        dS += dCTdS[i] * df('CT')[i] * 1e2/1e7
                    if self.get_id('Wac') in arguments:
                        dWac += dCTdWac[i] * df('CT')[i] * 1e6/1e7
                    if self.get_id('alpha') in arguments:
                        dalpha[i] += (dCTdAlpha2[i] * df('CT')[i] +
                                      dCTdAlpha1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] += (dCTdRho2[i] * df('CT')[i] +
                                    dCTdRho1[i+1] * df('CT')[i+1]) /1e7
                    if self.get_id('v') in arguments:
                        dv[i] += (dCTdV2[i] * df('CT')[i] +
                                  dCTdV1[i+1] * df('CT')[i+1]) * 1e2/1e7
                    if self.get_id('CD') in arguments:
                        dCD[i] += (dCTdCD2[i] * df('CT')[i] + 
                                   dCTdCD1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                    if self.get_id('Wf') in arguments:
                        dWf[i] += (dCTdWf2[i] * df('CT')[i] +
                                   dCTdWf1[i+1] * df('CT')[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        dgamma[i] += (dCTdGamma2[i] * df('CT')[i] +
                                      dCTdGamma1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        dCT[i] += (dCTdCT2[i] * df('CT')[i] + 
                                   dCTdCT1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                elif i == self.numElem:
                    if self.get_id('S') in arguments:
                        dS += dCTdS[i] * df('CT')[i] * 1e2/1e7
                    if self.get_id('Wac') in arguments:
                        dWac += dCTdWac[i] * df('CT')[i] * 1e6/1e7
                    if self.get_id('alpha') in arguments:
                        dalpha[i] += (dCTdAlpha3[i-1] * df('CT')[i-1] +
                                   dCTdAlpha2[i] * df('CT')[i]) * 1e-1/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] += (dCTdRho3[i-1] * df('CT')[i-1] + 
                                 dCTdRho2[i] * df('CT')[i]) /1e7
                    if self.get_id('v') in arguments:
                        dv[i] += (dCTdV3[i-1] * df('CT')[i-1] +
                               dCTdV2[i] * df('CT')[i]) * 1e2/1e7
                    if self.get_id('CD') in arguments:
                        dCD[i] += (dCTdCD3[i-1] * df('CT')[i-1] + 
                                dCTdCD2[i] * df('CT')[i]) * 1e-1/1e7
                    if self.get_id('Wf') in arguments:
                        dWf[i] += (dCTdWf3[i-1] * df('CT')[i-1] + 
                                dCTdWf2[i] * df('CT')[i]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        dgamma[i] += (dCTdGamma3[i-1] * df('CT')[i-1] + 
                                   dCTdGamma2[i] * df('CT')[i]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        dCT[i] += (dCTdCT3[i-1] * df('CT')[i-1] + 
                                dCTdCT2[i] * df('CT')[i]) * 1e-1/1e7
                else:
                    if self.get_id('S') in arguments:
                        dS += dCTdS[i] * df('CT')[i] * 1e2/1e7
                    if self.get_id('Wac') in arguments:
                        dWac += dCTdWac[i] * df('CT')[i] * 1e6/1e7
                    if self.get_id('alpha') in arguments:
                        dalpha[i] += (dCTdAlpha3[i-1] * df('CT')[i-1] +
                                   dCTdAlpha2[i] * df('CT')[i] +
                                   dCTdAlpha1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] += (dCTdRho3[i-1] * df('CT')[i-1] + 
                                 dCTdRho2[i] * df('CT')[i] +
                                 dCTdRho1[i+1] * df('CT')[i+1]) /1e7
                    if self.get_id('v') in arguments:
                        dv[i] += (dCTdV3[i-1] * df('CT')[i-1] +
                               dCTdV2[i] * df('CT')[i] +
                               dCTdV1[i+1] * df('CT')[i+1]) * 1e2/1e7
                    if self.get_id('CD') in arguments:
                        dCD[i] += (dCTdCD3[i-1] * df('CT')[i-1] + 
                                dCTdCD2[i] * df('CT')[i] +
                                dCTdCD1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                    if self.get_id('Wf') in arguments:
                        dWf[i] += (dCTdWf3[i-1] * df('CT')[i-1] + 
                                dCTdWf2[i] * df('CT')[i] +
                                dCTdWf1[i+1] * df('CT')[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        dgamma[i] += (dCTdGamma3[i-1] * df('CT')[i-1] + 
                                   dCTdGamma2[i] * df('CT')[i] +
                                   dCTdGamma1[i+1] * df('CT')[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        dCT[i] += (dCTdCT3[i-1] * df('CT')[i-1] + 
                                dCTdCT2[i] * df('CT')[i] +
                                dCTdCT1[i+1] * df('CT')[i+1]) * 1e-1/1e7
        self._lin_final()

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
                dcThrustSL[:] = numpy.sum(dtdcThrustSL * dg('tau')) * 1e6
            if self.get_id('h') in arguments:
                dh[:] = dtdh * dg('tau') * 1e3
            if self.get_id('CT') in arguments:
                dCT[:] = dtdCT * dg('tau') * 1e-1
            if self.get_id('rho') in arguments:
                drho[:] = dtdRho * dg('tau')
            if self.get_id('v') in arguments:
                dv[:] = dtdV * dg('tau') * 1e2
            if self.get_id('S') in arguments:
                dS[:] = numpy.sum(dtdS * dg('tau')) * 1e2

class Sys_CL(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        self.numInt = self.kwargs['numInt']
        CL_IC = self.kwargs['CL_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('CL', size=numPts, val=CL_IC)#,
                               #f_scal=1e7)
        self._declare_argument('Wf', indices=iPts)
        self._declare_argument('gamma', indices=iPts)
        self._declare_argument('CT', indices=iPts)
        self._declare_argument('alpha', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        self._declare_argument('v', indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['Wac',0], indices=[0])
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument('x', indices=iPts)
        
    def apply_F(self):
        self._nln_init()

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
        x = p('x') * 1e6
        R = f('CL')
        

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

            deltax = (x[i+1]-x[i])/self.numInt
            QTemp = 0.5*rhoTemp*vTemp*vTemp*S
            dxTemp = numpy.ones(self.numInt) * deltax
            dxTemp[0] = 0.5*deltax
            dxTemp[-1] = 0.5*deltax

            cosGamma = numpy.cos(gammaTemp)
            sinAlpha = numpy.sin(aTemp)

            CLResTemp = -QTemp*CLTemp + WTemp*cosGamma - CTTemp*sinAlpha*QTemp
            CLRes[i] += numpy.sum(CLResTemp * R1 * dxTemp)
            CLRes[i+1] += numpy.sum(CLResTemp * R2 * dxTemp)
        R[:] = CLRes / 1e7
            
        """
        CLRes = mission.get_cl(self.numElem, self.numInt, Wac, S, g, 
                               x_int, v, rho, CL, Wf, gamma, CT, alpha)

        #print "current CL: ", CL
        #print "current CLRes: ", CLRes
        f('CL')[:] = CLRes
        """
        self._nln_final()

    def apply_dFdpu(self, arguments):
        self._lin_init()
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']

        Wac = p(['Wac',0]) * 1e6
        S = p(['S',0]) * 1e2
        x = p('x') * 1e6
        v = p('v') * 1e2
        rho = p('rho')
        CL = u('CL')
        Wf = p('Wf') * 1e6
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        alpha = p('alpha') * 1e-1
        g = 9.81

        dWac = dp(['Wac', 0])
        dS = dp(['S', 0])
        dx = dp('x')
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

        if self.mode == 'fwd':
            df('CL')[:] = 0.0
            for i in xrange(self.numElem+1):
                if i == 0:
                    if self.get_id('Wac') in arguments:
                        df('CL')[0] += (dCLdWac[0] * dWac) * 1e6 / 1e7
                    if self.get_id('S') in arguments:
                        df('CL')[0] += dCLdS[0] * dS * 1e2 / 1e7
                    if self.get_id('v') in arguments:
                        df('CL')[0] += (dCLdV2[0] * dv[0] + 
                                        dCLdV3[0] * dv[1]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        df('CL')[0] += (dCLdRho2[0] * drho[0] +
                                        dCLdRho3[0] * drho[1]) / 1e7
                    if self.get_id('Wf') in arguments:
                        df('CL')[0] += (dCLdWf2[0] * dWf[0] + 
                                        dCLdWf3[0] * dWf[1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        df('CL')[0] += (dCLdGamma2[0] * dgamma[0] +
                                        dCLdGamma3[0] * dgamma[1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        df('CL')[0] += (dCLdThrust2[0] * dCT[0] +
                                        dCLdThrust3[0] * dCT[1]) * 1e-1/1e7
                    if self.get_id('alpha') in arguments:
                        df('CL')[0] += (dCLdAlpha2[0] * dalpha[0] +
                                        dCLdAlpha3[0] * dalpha[1]) * 1e-1/1e7
                    if self.get_id('CL') in arguments:
                        df('CL')[0] += (dCLdCL2[0] * du('CL')[0] + 
                                        dCLdCL3[0] * du('CL')[1]) /1e7
                elif i == self.numElem:
                    if self.get_id('Wac') in arguments:
                        df('CL')[i] += (dCLdWac[i] * dWac) * 1e6 / 1e7
                    if self.get_id('S') in arguments:
                        df('CL')[i] += dCLdS[i] * dS * 1e2 / 1e7
                    if self.get_id('v') in arguments:
                        df('CL')[i] += (dCLdV1[i] * dv[i-1] + 
                                        dCLdV2[i] * dv[i]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        df('CL')[i] += (dCLdRho1[i] * drho[i-1] + 
                                        dCLdRho2[i] * drho[i]) /1e7
                    if self.get_id('Wf') in arguments:
                        df('CL')[i] += (dCLdWf1[i] * dWf[i-1] + 
                                        dCLdWf2[i] * dWf[i]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        df('CL')[i] += (dCLdGamma1[i] * dgamma[i-1] + 
                                        dCLdGamma2[i] * dgamma[i]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        df('CL')[i] += (dCLdThrust1[i] * dCT[i-1] + 
                                        dCLdThrust2[i] * dCT[i]) * 1e-1/1e7
                    if self.get_id('alpha') in arguments:
                        df('CL')[i] += (dCLdAlpha1[i] * dalpha[i-1] +
                                        dCLdAlpha2[i] * dalpha[i]) * 1e-1/1e7
                    if self.get_id('CL') in arguments:
                        df('CL')[i] += (dCLdCL1[i] * du('CL')[i-1] + 
                                        dCLdCL2[i] * du('CL')[i]) /1e7
                else:
                    if self.get_id('Wac') in arguments:
                        df('CL')[i] += (dCLdWac[i] * dWac) * 1e6 / 1e7
                    if self.get_id('S') in arguments:
                        df('CL')[i] += dCLdS[i] * dS * 1e2 / 1e7
                    if self.get_id('v') in arguments:
                        df('CL')[i] += (dCLdV1[i] * dv[i-1] + 
                                        dCLdV2[i] * dv[i] + 
                                        dCLdV3[i] * dv[i+1]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        df('CL')[i] += (dCLdRho1[i] * drho[i-1] + 
                                        dCLdRho2[i] * drho[i] + 
                                        dCLdRho3[i] * drho[i+1]) /1e7
                    if self.get_id('Wf') in arguments:
                        df('CL')[i] += (dCLdWf1[i] * dWf[i-1] + 
                                        dCLdWf2[i] * dWf[i] + 
                                        dCLdWf3[i] * dWf[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        df('CL')[i] += (dCLdGamma1[i] * dgamma[i-1] + 
                                        dCLdGamma2[i] * dgamma[i] + 
                                        dCLdGamma3[i] * dgamma[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        df('CL')[i] += (dCLdThrust1[i] * dCT[i-1] + 
                                        dCLdThrust2[i] * dCT[i] +
                                        dCLdThrust3[i] * dCT[i+1]) * 1e-1/1e7
                    if self.get_id('alpha') in arguments:
                        df('CL')[i] += (dCLdAlpha1[i] * dalpha[i-1] + 
                                        dCLdAlpha2[i] * dalpha[i] +
                                        dCLdAlpha3[i] * dalpha[i+1]) * 1e-1/1e7
                    if self.get_id('CL') in arguments:
                        df('CL')[i] += (dCLdCL1[i] * du('CL')[i-1] + 
                                        dCLdCL2[i] * du('CL')[i] + 
                                        dCLdCL3[i] * du('CL')[i+1]) /1e7
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

            for i in xrange(self.numElem+1):
                if i == 0:
                    if self.get_id('Wac') in arguments:
                        dWac += dCLdWac[i] * df('CL')[i] * 1e6/1e7
                    if self.get_id('S') in arguments:
                        dS += dCLdS[i] * df('CL')[i] * 1e2/1e7
                    if self.get_id('v') in arguments:
                        dv[i] = (dCLdV2[i] * df('CL')[i] +
                                 dCLdV1[i+1] * df('CL')[i+1]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] = (dCLdRho2[i] * df('CL')[i] +
                                   dCLdRho1[i+1] * df('CL')[i+1]) /1e7
                    if self.get_id('Wf') in arguments:
                        dWf[i] = (dCLdWf2[i] * df('CL')[i] +
                                  dCLdWf1[i+1] * df('CL')[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        dgamma[i] = (dCLdGamma2[i] * df('CL')[i] +
                                     dCLdGamma1[i+1] * df('CL')[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        dCT[i] = (dCLdThrust2[i] * df('CL')[i] +
                                  dCLdThrust1[i+1] * df('CL')[i+1]) * 1e-1/1e7
                    if self.get_id('alpha') in arguments:
                        dalpha[i] = (dCLdAlpha2[i] * df('CL')[i] +
                                     dCLdAlpha1[i+1] * df('CL')[i+1]) * 1e-1/1e7
                    if self.get_id('CL') in arguments:
                        du('CL')[i] = (dCLdCL2[i] * df('CL')[i] +
                                       dCLdCL1[i+1] * df('CL')[i+1]) /1e7
                elif i == self.numElem:
                    if self.get_id('Wac') in arguments:
                        dWac += dCLdWac[i] * df('CL')[i] * 1e6/1e7
                    if self.get_id('S') in arguments:
                        dS += dCLdS[i] * df('CL')[i] * 1e2/1e7
                    if self.get_id('v') in arguments:
                        dv[i] = (dCLdV3[i-1] * df('CL')[i-1] +
                                 dCLdV2[i] * df('CL')[i]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] = (dCLdRho3[i-1] * df('CL')[i-1] +
                                   dCLdRho2[i] * df('CL')[i]) /1e7
                    if self.get_id('Wf') in arguments:
                        dWf[i] = (dCLdWf3[i-1] * df('CL')[i-1] +
                                  dCLdWf2[i] * df('CL')[i]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        dgamma[i] = (dCLdGamma3[i-1] * df('CL')[i-1] +
                                     dCLdGamma2[i] * df('CL')[i]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        dCT[i] = (dCLdThrust3[i-1] * df('CL')[i-1] +
                                  dCLdThrust2[i] * df('CL')[i]) * 1e-1/1e7
                    if self.get_id('alpha') in arguments:
                        dalpha[i] = (dCLdAlpha3[i-1] * df('CL')[i-1] +
                                     dCLdAlpha2[i] * df('CL')[i]) * 1e-1/1e7
                    if self.get_id('CL') in arguments:
                        du('CL')[i] = (dCLdCL3[i-1] * df('CL')[i-1] +
                                       dCLdCL2[i] * df('CL')[i]) /1e7
                else:
                    if self.get_id('Wac') in arguments:
                        dWac += dCLdWac[i] * df('CL')[i] * 1e6/1e7
                    if self.get_id('S') in arguments:
                        dS += dCLdS[i] * df('CL')[i] * 1e2/1e7
                    if self.get_id('v') in arguments:
                        dv[i] = (dCLdV3[i-1] * df('CL')[i-1] +
                                 dCLdV2[i] * df('CL')[i] +
                                 dCLdV1[i+1] * df('CL')[i+1]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] = (dCLdRho3[i-1] * df('CL')[i-1] +
                                   dCLdRho2[i] * df('CL')[i] +
                                   dCLdRho1[i+1] * df('CL')[i+1]) /1e7
                    if self.get_id('Wf') in arguments:
                        dWf[i] = (dCLdWf3[i-1] * df('CL')[i-1] +
                                  dCLdWf2[i] * df('CL')[i] +
                                  dCLdWf1[i+1] * df('CL')[i+1]) * 1e6/1e7
                    if self.get_id('gamma') in arguments:
                        dgamma[i] = (dCLdGamma3[i-1] * df('CL')[i-1] +
                                     dCLdGamma2[i] * df('CL')[i] +
                                     dCLdGamma1[i+1] * df('CL')[i+1]) * 1e-1/1e7
                    if self.get_id('CT') in arguments:
                        dCT[i] = (dCLdThrust3[i-1] * df('CL')[i-1] +
                                  dCLdThrust2[i] * df('CL')[i] +
                                  dCLdThrust1[i+1] * df('CL')[i+1]) * 1e-1/1e7
                    if self.get_id('alpha') in arguments:
                        dalpha[i] = (dCLdAlpha3[i-1] * df('CL')[i-1] +
                                     dCLdAlpha2[i] * df('CL')[i] +
                                     dCLdAlpha1[i+1] * df('CL')[i+1]) * 1e-1/1e7
                    if self.get_id('CL') in arguments:
                        du('CL')[i] = (dCLdCL3[i-1] * df('CL')[i-1] +
                                       dCLdCL2[i] * df('CL')[i] +
                                       dCLdCL1[i+1] * df('CL')[i+1]) /1e7
                
        self._lin_final()

class Sys_alpha(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        a_IC = self.kwargs['a_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('alpha', size=numPts, val=a_IC)#,
                               #u_scal=1e-1, f_scal=1e-1)
        self._declare_argument('CL', indices=iPts)
        self._declare_argument('eta', indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']

        CL = p('CL')
        eta = p('eta') * 1e-1

        alpha = mission.get_alpha(self.numElem, CL, eta)
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
                dCL[:] = dAdCL * dg('alpha') /1e-1

class Sys_CD(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CD_IC = self.kwargs['CD_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('CD', size=numPts, val=CD_IC)#,
                               #u_scal=1e-1, f_scal=1e-1)
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
                dAR[:] = numpy.sum(dCDdAR * dg('CD')) /1e-1
            if self.get_id('e') in arguments:
                de[:] = numpy.sum(dCDde * dg('CD')) /1e-1
            if self.get_id('CL') in arguments:
                dCL[:] = dCDdCL * dg('CD') /1e-1

class Sys_Wf(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        self.numInt = self.kwargs['numInt']
        Wf_IC = self.kwargs['Wf_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        numSeg = self.kwargs['numSeg']

        self._declare_variable('Wf', size=numPts, val=Wf_IC)#,
                               #u_scal=1e6, f_scal=1e6)
        self._declare_argument('v', indices=iPts)
        self._declare_argument('gamma', indices=iPts)
        self._declare_argument('CT', indices=iPts)
        self._declare_argument('x', indices=iPts)
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument('SFC', indices=iPts)
        self._declare_argument('rho', indices=iPts)
        #self._declare_argument(['WfSeg',0], indices=[self.copy+1])
        self._declare_argument(['S',0], indices=[0])
        if self.copy != numSeg-1:
            self._declare_argument(['Wf', self.copy+1], indices=[0])
        
    def apply_F(self):
        self._nln_init()
        numSeg = self.kwargs['numSeg']
        
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        x = p('x') * 1e6
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

        Wf = mission.get_wf(self.numInt, self.numElem, x, v, gamma,
                            CT, SFC, rho, WfIn, g, WfSeg, S)

        f('Wf')[:] = (WfIn-Wf) / 1e6
        self._nln_final()

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dFdpu(self, arguments):
        self._lin_init()
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
        x = p('x') * 1e6
        v = p('v') * 1e2
        gamma = p('gamma') * 1e-1
        CT = p('CT') * 1e-1
        SFC = p('SFC') * 1e-6
        rho = p('rho')

        #dWfSeg = dp(['WfSeg', 0])
        dS = dp(['S', 0])
        dx = dp('x')
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
        dWdCT = numpy.zeros(self.numElem+1)
        dWdSFC = numpy.zeros(self.numElem+1)
        dWdV = numpy.zeros(self.numElem+1)
        dWdGamma = numpy.zeros(self.numElem+1)
        dWdRho = numpy.zeros(self.numElem+1)
        dWdS = numpy.zeros(self.numElem+1)
        dWdWf = numpy.zeros(self.numElem+1)

        if self.mode == 'fwd':
            df('Wf')[:] = 0.0
            for i in xrange(self.numElem):
                j = self.numElem-1-i
                if self.get_id('CT') in arguments:
                    dWdCT[j] += dWdCT[j+1]
                    dWdCT[j] -= (dWfdCT1[j] * dCT[j] +
                                 dWfdCT2[j] * dCT[j+1])
                    df('Wf')[j] += dWdCT[j] * 1e-1/1e6
                if self.get_id('SFC') in arguments:
                    dWdSFC[j] += dWdSFC[j+1]
                    dWdSFC[j] -= (dWfdSFC1[j] * dSFC[j] + 
                                  dWfdSFC2[j] * dSFC[j+1])
                    df('Wf')[j] += dWdSFC[j] * 1e-6/1e6
                if self.get_id('v') in arguments:
                    dWdV[j] += dWdV[j+1]
                    dWdV[j] -= (dWfdV1[j] * dv[j] + 
                                dWfdV2[j] * dv[j+1])
                    df('Wf')[j] += dWdV[j] * 1e2/1e6
                if self.get_id('gamma') in arguments:
                    dWdGamma[j] += dWdGamma[j+1]
                    dWdGamma[j] -= (dWfdGamma1[j] * dgamma[j] + 
                                    dWfdGamma2[j] * dgamma[j+1])
                    df('Wf')[j] += dWdGamma[j] * 1e-1/1e6
                if self.get_id('rho') in arguments:
                    dWdRho[j] += dWdRho[j+1]
                    dWdRho[j] -= (dWfdRho1[j] * drho[j] + 
                                  dWfdRho2[j] * drho[j+1])
                    df('Wf')[j] += dWdRho[j] /1e6
                if self.get_id('S') in arguments:
                    dWdS[j] += dWdS[j+1]
                    dWdS[j] -= dWfdS[j] * dS
                    df('Wf')[j] += dWdS[j] * 1e2/1e6
            if self.get_id('Wf') in arguments:
                df('Wf')[:] += du('Wf')[:]

        if self.mode == 'rev':
            dCT[:] = 0.0
            dSFC[:] = 0.0
            dv[:] = 0.0
            dgamma[:] = 0.0
            dS[:] = 0.0
            drho[:] = 0.0
            #dWfSeg[:] = 0.0
            du('Wf')[:] = 0.0
            
            dWfSum = numpy.sum(df('Wf'))
            dWfSum -= df('Wf')[self.numElem]

            for i in xrange(self.numElem):
                j = self.numElem-1-i
                if self.get_id('CT') in arguments:
                    dCT[j+1] -= dWfdCT2[j] * dWfSum * 1e-1/1e6
                    dCT[j] -= dWfdCT1[j] * dWfSum * 1e-1/1e6
                if self.get_id('SFC') in arguments:
                    dSFC[j+1] -= dWfdSFC2[j] * dWfSum * 1e-6/1e6
                    dSFC[j] -= dWfdSFC1[j] * dWfSum * 1e-6/1e6
                if self.get_id('v') in arguments:
                    dv[j+1] -= dWfdV2[j] * dWfSum * 1e2/1e6
                    dv[j] -= dWfdV1[j] * dWfSum * 1e2/1e6
                if self.get_id('gamma') in arguments:
                    dgamma[j+1] -= dWfdGamma2[j] * dWfSum * 1e-1/1e6
                    dgamma[j] -= dWfdGamma1[j] * dWfSum * 1e-1/1e6
                if self.get_id('rho') in arguments:
                    drho[j+1] -= dWfdRho2[j] * dWfSum /1e6
                    drho[j] -= dWfdRho1[j] * dWfSum /1e6
                if self.get_id('S') in arguments:
                    dWdS[j] += dWdS[j+1]
                    dWdS[j] -= dWfdS[j]
                    dS += dWdS[j] * df('Wf')[j] * 1e2/1e6
                dWfSum -= df('Wf')[j]
            if self.get_id('Wf') in arguments:
                du('Wf')[:] += df('Wf')[:]

        self._lin_final()

class Sys_CM(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        self.numInt = self.kwargs['numInt']
        CM_IC = self.kwargs['CM_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable('CM', size=numPts, val=CM_IC)#,
                               #f_scal=1e7)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['chord',0], indices=[0])
        self._declare_argument('x', indices=iPts)
        self._declare_argument('v', indices=iPts)
        self._declare_argument('rho', indices=iPts)

    def apply_F(self):
        self._nln_init()
        
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        S = p(['S',0]) * 1e2
        chord = p(['chord',0])
        x = p('x') * 1e6
        v = p('v') * 1e2
        rho = p('rho')
        CM = u('CM')

        #print 'CM-v', v

        CMRes = mission.get_cm(self.numElem, self.numInt, S, chord,
                               x, v, rho, CM)
        #print "current CM: ", CM
        #print "current CMRes: ", CMRes

        f('CM')[:] = CMRes / 1e7
        self._nln_final()

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dFdpu(self, arguments):
        self._lin_init()
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']
        df = self.vec['df']

        S = p(['S',0]) * 1e2
        chord = p(['chord',0])
        x = p('x') * 1e6
        v = p('v') * 1e2
        rho = p('rho')
        CM = u('CM')

        dS = dp(['S', 0])
        dchord = dp(['chord', 0])
        dx = dp('x')
        dv = dp('v')
        drho = dp('rho')
        
        [dCMdS, dCMdC, dCMdV1, dCMdV2, dCMdV3, dCMdRho1, dCMdRho2,
         dCMdRho3, dCMdCM1, dCMdCM2,
         dCMdCM3] = mission.get_cm_d(self.numElem, self.numInt, S,
                                     chord, x, v, rho, CM)

        if self.mode == 'fwd':
            df('CM')[:] = 0.0
            for i in xrange(self.numElem+1):
                if i == 0:
                    if self.get_id('S') in arguments:
                        df('CM')[i] += (dCMdS[i] * dS) * 1e2/1e7
                    if self.get_id('chord') in arguments:
                        df('CM')[i] += dCMdC[i] * dchord /1e7
                    if self.get_id('v') in arguments:
                        df('CM')[i] += (dCMdV2[i] * dv[i] +
                                        dCMdV3[i] * dv[i+1]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        df('CM')[i] += (dCMdRho2[i] * drho[i] +
                                        dCMdRho3[i] * drho[i+1]) /1e7
                    if self.get_id('CM') in arguments:
                        df('CM')[i] += (dCMdCM2[i] * du('CM')[i] + 
                                        dCMdCM3[i] * du('CM')[i+1]) /1e7
                elif i == self.numElem:
                    if self.get_id('S') in arguments:
                        df('CM')[i] += (dCMdS[i] * dS) * 1e2/1e7
                    if self.get_id('chord') in arguments:
                        df('CM')[i] += dCMdC[i] * dchord /1e7
                    if self.get_id('v') in arguments:
                        df('CM')[i] += (dCMdV1[i] * dv[i-1] + 
                                        dCMdV2[i] * dv[i]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        df('CM')[i] += (dCMdRho1[i] * drho[i-1] +
                                        dCMdRho2[i] * drho[i]) /1e7
                    if self.get_id('CM') in arguments:
                        df('CM')[i] += (dCMdCM1[i] * du('CM')[i-1] + 
                                        dCMdCM2[i] * du('CM')[i]) /1e7
                else:
                    if self.get_id('S') in arguments:
                        df('CM')[i] += (dCMdS[i] * dS) * 1e2/1e7
                    if self.get_id('chord') in arguments:
                        df('CM')[i] += dCMdC[i] * dchord /1e7
                    if self.get_id('v') in arguments:
                        df('CM')[i] += (dCMdV1[i] * dv[i-1] + 
                                        dCMdV2[i] * dv[i] +
                                        dCMdV3[i] * dv[i+1]) * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        df('CM')[i] += (dCMdRho1[i] * drho[i-1] +
                                        dCMdRho2[i] * drho[i] +
                                        dCMdRho3[i] * drho[i+1]) /1e7
                    if self.get_id('CM') in arguments:
                        df('CM')[i] += (dCMdCM1[i] * du('CM')[i-1] + 
                                        dCMdCM2[i] * du('CM')[i] +
                                        dCMdCM3[i] * du('CM')[i+1]) /1e7
        if self.mode == 'rev':
            dS[:] = 0.0
            dchord[:] = 0.0
            dv[:] = 0.0
            drho[:] = 0.0
            du('CM')[:] = 0.0

            for i in xrange(self.numElem+1):
                if i == 0:
                    if self.get_id('S') in arguments:
                        dS += (dCMdS[i] * df('CM')[i]) * 1e2/1e7
                    if self.get_id('chord') in arguments:
                        dchord += dCMdC[i] * df('CM')[i] /1e7
                    if self.get_id('v') in arguments:
                        dv[i] += dCMdV2[i] * df('CM')[i] * 1e2/1e7
                        dv[i+1] += dCMdV3[i] * df('CM')[i] * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        drho[i] += dCMdRho2[i] * df('CM')[i] /1e7
                        drho[i+1] += dCMdRho3[i] * df('CM')[i] /1e7
                    if self.get_id('CM') in arguments:
                        du('CM')[i] += dCMdCM2[i] * df('CM')[i] /1e7 
                        du('CM')[i+1] += dCMdCM3[i] * df('CM')[i] /1e7
                elif i == self.numElem:
                    if self.get_id('S') in arguments:
                        dS += (dCMdS[i] * df('CM')[i]) * 1e2/1e7
                    if self.get_id('chord') in arguments:
                        dchord += dCMdC[i] * df('CM')[i] /1e7
                    if self.get_id('v') in arguments:
                        dv[i-1] += dCMdV1[i] * df('CM')[i] * 1e2/1e7 
                        dv[i] += dCMdV2[i] * df('CM')[i] * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        drho[i-1] += dCMdRho1[i] * df('CM')[i] /1e7
                        drho[i] += dCMdRho2[i] * df('CM')[i] /1e7
                    if self.get_id('CM') in arguments:
                        du('CM')[i-1] += dCMdCM1[i] * df('CM')[i] /1e7 
                        du('CM')[i] += dCMdCM2[i] * df('CM')[i] /1e7
                else:
                    if self.get_id('S') in arguments:
                        dS += (dCMdS[i] * df('CM')[i]) * 1e2/1e7
                    if self.get_id('chord') in arguments:
                        dchord += dCMdC[i] * df('CM')[i] /1e7
                    if self.get_id('v') in arguments:
                        dv[i-1] += dCMdV1[i] * df('CM')[i] * 1e2/1e7 
                        dv[i] += dCMdV2[i] * df('CM')[i] * 1e2/1e7
                        dv[i+1] += dCMdV3[i] * df('CM')[i] * 1e2/1e7
                    if self.get_id('rho') in arguments:
                        drho[i-1] += dCMdRho1[i] * df('CM')[i] /1e7
                        drho[i] += dCMdRho2[i] * df('CM')[i] /1e7
                        drho[i+1] += dCMdRho3[i] * df('CM')[i] /1e7
                    if self.get_id('CM') in arguments:
                        du('CM')[i-1] += dCMdCM1[i] * df('CM')[i] /1e7 
                        du('CM')[i] += dCMdCM2[i] * df('CM')[i] /1e7
                        du('CM')[i+1] += dCMdCM3[i] * df('CM')[i] /1e7

        self._lin_final()

class Sys_eta(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        e_IC = self.kwargs['e_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable('eta', size=numPts, val=e_IC)#,
                               #u_scal=1e-1, f_scal=1e-1)
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
                dCM[:] = dEdCM * dg('eta') /1e-1
            if self.get_id('alpha') in arguments:
                dalpha[:] = dEdAlpha * dg('eta')

class Sys_Wf_obj(ExplicitSystem):
    def _declare(self):
        self._declare_variable('Wf_obj', size=1)
        self._declare_argument(['Wf', 0], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        Wf = p('Wf')

        u('Wf_obj')[0] = Wf[0]

    def apply_dGdp(self, arguments):
        dp = self.vec['dp']
        dg = self.vec['dg']
        
        if self.mode == 'fwd':
            dg('Wf_obj')[0] = 0.0
            if self.get_id('Wf') in arguments:
                dg('Wf_obj')[0] += dp('Wf')[0]
        if self.mode == 'rev':
            dp('Wf')[0] = 0.0
            if self.get_id('Wf') in arguments:
                dp('Wf')[0] += dg('Wf_obj')[0]

class Sys_h_i(ExplicitSystem):
    def _declare(self):
        self._declare_variable('h_i')
        self._declare_argument('h', indices=[0])
    
    def apply_G(self):
        u = self.vec['u']('h_i')
        p = self.vec['p']('h')
        
        u[0] = p[0]

    def apply_dGdp(self, args):
        du = self.vec['dg']('h_i')
        dp = self.vec['dp']('h')
        
        if self.mode == 'fwd':
            du[0] = 0.0
            if self.get_id('h') in args:
                du[0] += dp[0]
        if self.mode == 'rev':
            dp[0] = 0.0
            if self.get_id('h') in args:
                dp[0] += du[0]

class Sys_h_f(ExplicitSystem):
    def _declare(self):
        numElem = self.kwargs['numElem']
        self._declare_variable('h_f')
        self._declare_argument('h', indices=[numElem])
    
    def apply_G(self):
        u = self.vec['u']('h_f')
        p = self.vec['p']('h')
        
        u[0] = p[0]

    def apply_dGdp(self, args):
        du = self.vec['dg']('h_f')
        dp = self.vec['dp']('h')
        
        if self.mode == 'fwd':
            du[0] = 0.0
            if self.get_id('h') in args:
                du[0] += dp[0]
        if self.mode == 'rev':
            dp[0] = 0.0
            if self.get_id('h') in args:
                dp[0] += du[0]

class Tmin(ExplicitSystem):

    def _declare(self):
        numElem = self.kwargs['numElem']
        self._declare_variable('Tmin')
        self._declare_argument('tau', indices=range(0,numElem+1))
        self.min = 0.01
        self.rho = 30

    def apply_G(self):
        Tmin = self.vec['u']('Tmin')
        tau = self.vec['p']('tau')

        fmax = numpy.max(self.min - tau)
        Tmin[0] = fmax + 1/self.rho * numpy.log(numpy.sum(numpy.exp(self.rho*(self.min - tau - fmax))))

    def apply_dGdp(self, args):
        Tmin = self.vec['u']('Tmin')
        tau = self.vec['p']('tau')

        dTmin = self.vec['dg']('Tmin')
        dtau = self.vec['dp']('tau')

        ind = numpy.argmax(self.min - tau)
        fmax = self.min - tau[ind]
        dfmax_dtau = numpy.zeros(tau.shape[0])
        dfmax_dtau[ind] = -1.0

        deriv = dfmax_dtau + 1/self.rho * \
            1/numpy.sum(numpy.exp(self.rho*(self.min - tau - fmax))) * \
            numpy.exp(self.rho*(self.min - tau - fmax)) * (-self.rho)
        deriv[ind] -= 1/self.rho * \
            1/numpy.sum(numpy.exp(self.rho*(self.min - tau - fmax))) * \
            numpy.sum(numpy.exp(self.rho*(self.min - tau - fmax))) * (-self.rho)
        
        if self.mode == 'fwd':
            dTmin[0] = 0.0
            if self.get_id('tau') in args:
                dTmin[0] += numpy.sum(deriv * dtau[:])
        if self.mode == 'rev':
            dtau[:] = 0.0
            if self.get_id('tau') in args:
                dtau[:] += deriv * dTmin[0]

class Tmax(ExplicitSystem):

    def _declare(self):
        numElem = self.kwargs['numElem']
        self._declare_variable('Tmax')
        self._declare_argument('tau', indices=range(0,numElem+1))
        self.max = 1.0
        self.rho = 30

    def apply_G(self):
        Tmax = self.vec['u']('Tmax')
        tau = self.vec['p']('tau')

        fmax = numpy.max(tau - self.max)
        Tmax[0] = fmax + 1/self.rho * numpy.log(numpy.sum(numpy.exp(self.rho*(tau - self.max - fmax))))

    def apply_dGdp(self, args):
        Tmax = self.vec['u']('Tmax')
        tau = self.vec['p']('tau')

        dTmax = self.vec['dg']('Tmax')
        dtau = self.vec['dp']('tau')

        ind = numpy.argmax(tau - self.max)
        fmax = tau[ind] - self.max
        dfmax_dtau = numpy.zeros(tau.shape[0])
        dfmax_dtau[ind] = 1.0

        deriv = dfmax_dtau + 1/self.rho * \
            1/numpy.sum(numpy.exp(self.rho*(tau - self.max - fmax))) * \
            numpy.exp(self.rho*(tau - self.max - fmax)) * (self.rho)
        deriv[ind] -= 1/self.rho * \
            1/numpy.sum(numpy.exp(self.rho*(tau - self.max - fmax))) * \
            numpy.sum(numpy.exp(self.rho*(tau - self.max - fmax))) * (self.rho)

        if self.mode == 'fwd':
            dTmax[0] = 0.0
            if self.get_id('tau') in args:
                dTmax[0] += numpy.sum(deriv * dtau[:])
        if self.mode == 'rev':
            dtau[:] = 0.0
            if self.get_id('tau') in args:
                dtau[:] += deriv * dTmax[0]      


class Trajectory(object):

    def __init__(self, numSeg, opt=False):
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
        self.opt = opt
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

    def set_IC(self, h_IC=None, a_IC=None, t_IC=None, e_IC=None, CT_IC=None,
               CL_IC=None, CD_IC=None, CM_IC=None, Wf_IC=None, x_IC=None,
               v_IC=None):
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

        if CT_IC == None:
            self.CT_IC = numpy.ones(self.numElem+1)
        else:
            self.CT_IC = CT_IC/1e-1
        
        h = numpy.array(self.h_IC) * 1e3
        SFC = numpy.array(self.SFCSL) * 1e-6
        x = numpy.array(self.x_IC) * 1e6
        
        self.SFC_IC = mission.get_sfc(totalElem, SFC, h)/1e-6
        # potential bug: must allow diff x spacing in diff segments
        self.gamma_IC = mission.get_gamma(totalElem, h, x)/1e-1
        self.Temp_IC = mission.get_temp(totalElem, h)/1e2
        temp = numpy.array(self.Temp_IC) * 1e2
        self.rho_IC = mission.get_rho(totalElem, self.g, temp)
        
        if v_IC == None:
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
        else:
            self.v_IC = v_IC/1e2


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
                                                     self.v_IC,
                                                     self.h_IC))
                self.t_IC.pop()
        self.t_IC.extend(numpy.ones(1)*1.0)

    def set_range(self, mission_range):
        self.range = mission_range

    def set_final_Wf(self, Wf):
        self.Wf_final = Wf

    def set_ingn_intl(self, numInt):
        self.numInt = numInt

    def set_opt(self, dist, numElem, numCP):
        self.add_seg_point(0.0,0.0,v=150.0)
        self.add_seg_point(0.0,dist,v=150.0,tau=0.5,numElem=numElem)
        self.range = dist
        self.numCP = numCP

    def MBI(self):
        n = self.numElem+1
        m = self.numCP

        h = numpy.linspace(0,13,n)
        x = numpy.linspace(0,self.range,n)/1e6

        a = MBI.MBI(h,[x],[m],[4])
        J = a.getJacobian(0,0)
        Jd = a.getJacobian(1,0)

        Cx = self.range/1e6 * 0.5*(1-numpy.cos(numpy.pi*numpy.linspace(0,1,m)))
        #Cx = self.range/1e6 * numpy.linspace(0,1,m)
        dx = Jd.dot(Cx)*1e6

        lins = numpy.linspace(0, n-1, n).astype(int)
        diag = scipy.sparse.csc_matrix((1.0/dx,
                                        (lins,lins)))
        Je = diag.dot(Jd)

        return J, Je, Cx

    def initialize(self):
        ones = numpy.ones
        ls = numpy.linspace
        ne = self.numElem
        pt1 = 0
        pt2 = 0

        jac_h, jac_gamma, Cx = self.MBI()
        numCP = self.numCP

        h_CP_IC = numpy.ones(numCP)
        h_CP_IC[0] = 0
        h_CP_IC[-1] = 0

        if self.opt == False:
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
                                    IndVar('h_ends',pt,val=numpy.array([self.hPts[pt],self.hPts[pt+1]]),size=2),
                                    #u_scal=1e3, f_scal=1e3),
                                    IndVar('v_ends',pt,val=numpy.array([self.vPts[pt],self.vPts[pt+1]]),size=2),
                                    #u_scal=1e2, f_scal=1e2),
                                    IndVar('M_ends',pt,val=numpy.array([self.MPts[pt],self.MPts[pt+1]]),size=2),
                                    IndVar('h_dot',pt,val=self.hDotPts[pt],size=1),
                                    #u_scal=1e1, f_scal=1e1),
                                    IndVar('tau_init',pt,val=self.tPts[pt],size=1),
                                    ]),
                            SerialSystem('seg_analysis', pt,
                                         #NL='NEWTON',
                                         NL='NLN_GS',
                                         LN='KSP_PC',
                                         #PC='LIN_GS',
                                         LN_ilimit=100,
                                         LN_rtol=1e-10,
                                         LN_atol=1e-14,
                                         NL_ilimit=100,#10,
                                         NL_rtol=1e-13,
                                         NL_atol=1e-13,
                                         PC_ilimit=5,
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
                                    ]),
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
                                IndVar('S',val=self.S,size=1),# u_scal=1e2, f_scal=1e2),
                                IndVar('Wac',val=self.Wac,size=1),# u_scal=1e6, f_scal=1e6),
                                IndVar('cThrustSL',val=self.cThrustSL,size=1),# u_scal=1e6, f_scal=1e6),
                                IndVar('SFCSL',val=self.SFCSL,size=1),# u_scal=1e-6, f_scal=1e-6),
                                IndVar('chord',val=self.chord,size=1),
                                #IndVar('inertia',val=self.inertia,size=1, u_scal=1e7, f_scal=1e7),
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
                                     )
                        ]).setup()
        
        else:
            pt = 0
            h_IC = self.h_IC
            t_IC = self.t_IC
            e_IC = self.e_IC
            a_IC = self.a_IC
            CL_IC = self.CL_IC
            CD_IC = self.CD_IC
            CM_IC = self.CM_IC
            Wf_IC = self.Wf_IC
            SFC_IC = self.SFC_IC
            gamma_IC = self.gamma_IC
            Temp_IC = self.Temp_IC
            rho_IC = self.rho_IC
            v_IC = self.v_IC
            CT_IC = self.CT_IC
            
            self.mainMission = Top('mission',
                                   NL='NLN_GS', 
                                   LN='LIN_GS',
                                   LN_ilimit=5, 
                                   NL_ilimit=5, 
                                   NL_rtol=1e-6,
                                   NL_atol=1e-10,
                                   LN_rtol=1e-6,
                                   LN_atol=1e-10,
                                   output=True,
                                   subsystems=[
                    SerialSystem('mission_param',
                                 NL='NLN_GS',
                                 LN='LIN_GS',
                                 LN_ilimit=5, 
                                 NL_ilimit=5, 
                                 output=True,
                                 subsystems=[
                            IndVar('S',val=self.S,size=1),# u_scal=1e2, f_scal=1e2),
                            IndVar('Wac',val=self.Wac,size=1),# u_scal=1e6, f_scal=1e6),
                            IndVar('cThrustSL',val=self.cThrustSL,size=1),# u_scal=1e6, f_scal=1e6),
                            IndVar('SFCSL',val=self.SFCSL,size=1),# u_scal=1e-6, f_scal=1e-6),
                            IndVar('chord',val=self.chord,size=1),
                            #IndVar('inertia',val=self.inertia,size=1, u_scal=1e7, f_scal=1e7),
                            IndVar('AR',val=self.AR,size=1),
                            IndVar('e',val=self.e,size=1),
                            IndVar('g',val=self.g,size=1),
                            ]),
                    SerialSystem('mission_analysis',
                                 NL='NLN_GS', 
                                 LN='KSP_PC',
                                 PC='LIN_GS',
                                 LN_ilimit=30, 
                                 NL_ilimit=30,
                                 PC_ilimit=5,
                                 #NL_rtol=1e-10,
                                 #NL_atol=1e-10,
                                 #LN_rtol=1e-10,
                                 #LN_atol=1e-10,
                                 PC_rtol=1e-6,
                                 PC_atol=1e-10,
                                 output=True,
                                 subsystems=[
                            #IndVar('x', val=numpy.linspace(0.0,self.range,self.numElem[pt]+1)/1e6),
                            IndVar('x_CP', val=Cx, lower=0),
                            #IndVar('h', val=numpy.zeros(self.numElem+1)),
                            IndVar('h_CP', val=h_CP_IC, lower=0),
                            IndVar('v_CP', val=200/1e2, size=numCP, lower=0),
                            Sys_h_bspline('h', numPts=self.numElem[pt]+1, numCP=numCP, jac=jac_h),
                            Sys_x_bspline('x', numPts=self.numElem[pt]+1, numCP=numCP, jac=jac_h),
                            Sys_v_bspline('v', numPts=self.numElem[pt]+1, numCP=numCP, jac=jac_h),
                            Sys_SFC('SFC', numElem=self.numElem[pt], SFC_IC=SFC_IC),
                            Sys_gamma_bspline('gamma', numPts=self.numElem[pt]+1, numCP=numCP, jac=jac_gamma),
                            #Sys_gamma('gamma', numElem=self.numElem[pt], gamma_IC=gamma_IC),
                            Sys_Temp('Temp', numElem=self.numElem[pt], Temp_IC=Temp_IC),
                            Sys_rho('rho', numElem=self.numElem[pt], rho_IC=rho_IC),
                            Sys_CL('CL', numElem=self.numElem[pt], numInt=self.numInt, CL_IC=CL_IC),
                            Sys_alpha('alpha', numElem=self.numElem[pt], a_IC=a_IC),
                            Sys_CD('CD', numElem=self.numElem[pt], CD_IC=CD_IC),
                            Sys_CT_opt('CT', numElem=self.numElem[pt], numInt=self.numInt, CT_IC=CT_IC),
                            Sys_tau('tau', numElem=self.numElem[pt]),
                            Sys_Wf('Wf', numElem=self.numElem[pt], numInt=self.numInt, numSeg=self.numSeg, Wf_IC=Wf_IC),
                            Sys_CM('CM', numElem=self.numElem[pt], numInt=self.numInt, CM_IC=CM_IC),
                            Sys_eta('eta', numElem=self.numElem[pt], e_IC=e_IC),
                            Sys_Wf_obj('Wf_obj')
                            
                            ]
                                 ),
                    Sys_h_i('h_i'),
                    Sys_h_f('h_f', numElem=self.numElem[pt]),
                    Tmin('Tmin', numElem=self.numElem[pt]),
                    Tmax('Tmax', numElem=self.numElem[pt]),
                    ]).setup()
            self.mainMission.reset_counter()
        
        return self.mainMission


class Top(SerialSystem):

    def reset_counter(self):
        self.counter = 0

    def compute(self, output=False):
        temp, success = super(Top, self).compute(output)
        fig = matplotlib.pylab.figure(figsize=(12.0,10.0))
        v = self.vec['u']
        fig.add_subplot(711).plot(v('x')*1000.0, v('h'))
        fig.add_subplot(711).set_ylabel('Altitude (km)')
        fig.add_subplot(712).plot(v('x')*1000.0, v('v')*1e2)
        fig.add_subplot(712).set_ylabel('Velocity (m/s)')
        fig.add_subplot(713).plot(v('x')*1000.0, v('alpha')*1e-1*180.0/numpy.pi)
        fig.add_subplot(713).set_ylabel('AoA (deg)')
        fig.add_subplot(714).plot(v('x')*1000.0, v('tau'))
        fig.add_subplot(714).set_ylabel('Throttle')
        fig.add_subplot(715).plot(v('x')*1000.0, v('eta')*1e-1*180.0/numpy.pi)
        fig.add_subplot(715).set_ylabel('Trim Angle (deg)')
        fig.add_subplot(716).plot(v('x')*1000.0, v('Wf')*1e6/(9.81*0.804))
        fig.add_subplot(716).set_ylabel('Fuel (L)')
        fig.add_subplot(717).plot(v('x')*1000.0, v('rho'))
        fig.add_subplot(717).set_ylabel('rho')
        fig.add_subplot(717).set_xlabel('Distance (km)')
        fig.savefig("plots/OptFig_%i.pdf"%(self.counter))
        fig.savefig("plots/OptFig_%i.png"%(self.counter))
        self.counter += 1

        return temp, success
