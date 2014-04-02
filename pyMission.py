from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
import numpy
import copy
from framework import *
import mission


class Sys_h(ImplicitSystem):
    def _declare(self):
        h_IC = self.kwargs['h_IC']
        self.numPts = self.kwargs['numElem']+1
        self.numInt = self.kwargs['numInt']
        lower = 0.0
        iPts = range(self.numPts)

        self._declare_variable(['h',-1], size=self.numPts, 
                               val=h_IC, lower=lower)
        self._declare_argument(['h_ends',-1], indices=[0,1])
        #self._declare_argument(['v_ends',-1], indices=[0,1])
        self._declare_argument(['x',0], indices=[self.copy,
                                                 self.copy+1])
        #self._declare_argument(['h_dot',-1], indices=[0])
        self._declare_argument(['Wf',-1], indices=iPts)
        self._declare_argument(['CT',-1], indices=iPts)
        self._declare_argument(['alpha',-1], indices=iPts)
        self._declare_argument(['CD',-1], indices=iPts)
        self._declare_argument(['rho',-1], indices=iPts)
        self._declare_argument(['v',-1], indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['Wac',0], indices=[0])

    def apply_F(self):
        self._nln_init()
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
            h_ends = p(['h_ends',-1])
            #v_ends = p(['v_ends',-1])
            #h_dot = p(['h_dot',-1])
            x_ends = p(['x',0])
            Wf = p(['Wf',-1])
            CT = p(['CT',-1])
            alpha = p(['alpha',-1])
            CD = p(['CD',-1])
            rho = p(['rho',-1])
            v = p(['v',-1])
            S = p(['S',0])
            h = u(['h',-1])
            Wac = p(['Wac',0])

            h_new = mission.get_h(self.numPts-1, self.numInt, S, Wac, x_ends, 
                                  h_ends, Wf, CT, alpha, CD, rho, v)
            
            for i in xrange(self.numPts):
                if h_new[i] < 0:
                    h_new[i] = 0
                    #raise Exception('negative altitude!')
            
            #print "current h: ", h
            f(['h',-1])[:] = h - h_new
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
            S = p(['S',0])
            x_ends = p(['x',0])
            h_ends = p(['h_ends',-1])
            Wf = p(['Wf',-1])
            CT = p(['CT',-1])
            alpha = p(['alpha',-1])
            CD = p(['CD',-1])
            rho = p(['rho',-1])
            v = p(['v',-1])
            Wac = p(['Wac',0])
            
            [dhdx0, dhdx1, dhdWf1, dhdWf2, dhdRho1, dhdRho2, dhdV1, dhdV2,
             dhdCT1, dhdCT2, dhdA1, dhdA2, dhdCD1, 
             dhdCD2] = mission.get_h_d(numElem, numInt, S, Wac,  x_ends, h_ends,
                                       Wf, CT, alpha, CD, rho, v)
             
            if self.mode == 'fwd':
                df(['h',-1])[:] = 0.0
                
                #print 'dpx:',dp(['x',0])
                #print 'dpWf:',dp(['Wf',-1])
                #print 'dpRho:',dp(['rho',-1])
                #print 'dpV:',dp(['v',-1])
                #print 'dpCT:',dp(['CT',-1])
                #print 'dpAlpha:',dp(['alpha',-1])
                #print 'dpCD:',dp(['CD',-1])
                #print 'du:',du(['h',-1])
                
                if self.get_id('x') in arguments:
                    df(['h',-1])[:] += dhdx0*dp(['x',0])[0] + dhdx1*dp(['x',0])[1]
                    #print 'dpx0:',dp(['x',0])[0],'dpx1:',dp(['x',0])[1]
                    #print "dhdx0: ", dhdx0*dp(['x',0])[0], "dhdx1: ", dhdx1*dp(['x',0])[1]
                    #print 'dhdx0', dhdx0, 'dhdx1',dhdx1
                    #print 'dhdx',dhdx0*dp(['x',0])[0] + dhdx1*dp(['x',0])[1]
                if self.get_id('Wf') in arguments:
                    df(['h',-1])[:] -= (dhdWf1+dhdWf2)*dp(['Wf',-1])
                    #print 'dhdWf',(dhdWf1+dhdWf2)*dp(['Wf',-1])
                if self.get_id('rho') in arguments:
                    df(['h',-1])[:] -= (dhdRho1+dhdRho2)*dp(['rho',-1])
                    #print 'dhdRho',(dhdRho1+dhdRho2)*dp(['rho',-1])
                if self.get_id('v') in arguments:
                    df(['h',-1])[:] -= (dhdV1+dhdV2)*dp(['v',-1])
                    #print 'dhdV',(dhdV1+dhdV2)*dp(['v',-1])
                if self.get_id('CT') in arguments:
                    df(['h',-1])[:] -= (dhdCT1+dhdCT2)*dp(['CT',-1])
                    #print 'dhdCT',(dhdCT1+dhdCT2)*dp(['CT',-1])
                if self.get_id('alpha') in arguments:
                    df(['h',-1])[:] -= (dhdA1+dhdA2)*dp(['alpha',-1])
                    #print 'dhdA',(dhdA1+dhdA2)*dp(['alpha',-1])
                if self.get_id('CD') in arguments:
                    df(['h',-1])[:] -= (dhdCD1+dhdCD2)*dp(['CD',-1])
                    #print 'dhdCD',(dhdCD1+dhdCD2)*dp(['CD',-1])
                if self.get_id('h_ends') in arguments:
                    df(['h',-1])[:] -= numpy.ones(self.numPts)*dp(['h_ends',-1])[0]
                    #print 'dhdh', numpy.sum(dhdh*du(['h',-1]),axis=1)
                if self.get_id('h') in arguments:
                    df(['h',-1])[:] += du(['h',-1])
                #print 'DH:',df(['h',-1])

            if self.mode == 'rev':
                if self.get_id('x') in arguments:
                    dp(['x',0])[0] = dhdx0*df(['h',-1])
                    dp(['x',0])[1] = dhdx1*df(['h',-1])
                if self.get_id('Wf') in arguments:
                    dp(['Wf',-1])[:] = [dhdWf1, dhdWf2]*df(['h',-1])
                if self.get_id('rho') in arguments:
                    dp(['rho',-1])[:] = [dhdRho1, dhdRho2]*df(['h',-1])
                if self.get_id('v') in arguments:
                    dp(['v',-1])[:] = [dhdV1, dhdV2]*df(['h',-1])
                if self.get_id('CT') in arguments:
                    dp(['CT',-1])[:] = [dhdCT1, dhdCT2]*df(['h',-1])
                if self.get_id('alpha') in arguments:
                    dp(['alpha',-1])[:] = [dhdA1, dhdA2]*df(['h',-1])
                if self.get_id('CD') in arguments:
                    dp(['CD',-1])[:] = [dhdCD1, dhdCD2]*df(['h',-1])
        self._lin_final()

class Sys_SFC(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        SFC_IC = self.kwargs['SFC_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable(['SFC',-1], size=numPts, val=SFC_IC, 
                               u_scal=1e-5, f_scal=1e-5)
        self._declare_argument(['h',-1], indices=iPts)
        self._declare_argument(['SFCSL',0], indices=[0])

    def apply_G(self):
        self._nln_init()
        p = self.vec['p']
        u = self.vec['u']
        h = p(['h',-1])
        SFCSL = p(['SFCSL',0])
        SFC = mission.get_sfc(self.numElem, SFCSL, h)
        u(['SFC',-1])[:] = SFC
        self._nln_final()
    '''
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    '''
    def apply_dGdp(self, arguments):
        self._lin_init()
        mode = self.mode
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        numPts = self.numElem+1
        dSFC_dh = numpy.zeros(numPts)        
        dSFC_dh = mission.get_sfc_d(self.numElem)

        if self.mode == 'fwd':
            dg(['SFC',-1])[:] = 0.0
            if self.get_id('h') in arguments:
                dg(['SFC',-1])[:] += dSFC_dh * dp(['h',-1])
            if self.get_id('SFCSL') in arguments:
                dg(['SFC',-1])[:] += 1.0

        if self.mode == 'rev':
            if self.get_id('h') in arguments:
                dp(['h',-1])[:] = dSFC_dh * dg(['SFC',-1])
            if self.get_id('SFC') in arguments:
                du(['SFC',-1])[:] = 0.0
        self._lin_final()
    
class Sys_gamma(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        gamma_IC = self.kwargs['gamma_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['gamma',-1], size=numPts, val=gamma_IC)
        self._declare_argument(['h',-1], indices=iPts)
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        
    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p(['h',-1])
        x_ends = p(['x',0])
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)
        gamma = mission.get_gamma(self.numElem, h, x_int)
        u(['gamma',-1])[:] = gamma
    '''
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    '''
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
        h = p(['h',-1])
        x_ends = p(['x',0])
        x = numpy.linspace(x_ends[0], x_ends[1], numPts)

        [dgamma_dh1, dgamma_dh2, dgamma_dh3, dgamma_dh4, dgamma_dh5,
         dgamma_dx0, dgamma_dx1] = mission.get_gamma_d(numPts-1, h, x)

        if self.mode == 'fwd':
            dg(['gamma',-1])[:] = 0.0
            if self.get_id('h') in arguments:
                dg(['gamma',-1])[:] += dgamma_dh1*dp(['h',-1]) + \
                    dgamma_dh2*dp(['h',-1]) + \
                    dgamma_dh3*dp(['h',-1]) + \
                    dgamma_dh4*dp(['h',-1]) + \
                    dgamma_dh5*dp(['h',-1])
            if self.get_id('x') in arguments:
                dg(['gamma',-1])[:] += dgamma_dx0*dp(['x',0])[0] + \
                    dgamma_dx1*dp(['x',0])[1]

        if self.mode == 'rev':
            if self.get_id('h') in arguments:
                dp(['h',-1])[:] = [dgamma_dh1, dgamma_dh2, dgamma_dh3,
                                   dgamma_dh4, dgamma_dh5]*dg(['gamma',-1])
            if self.get_id('x') in arguments:
                dp(['x',0])[0] = dgamma_dx0 * dg(['gamma',-1])
                dp(['x',0])[1] = dgamma_dx1 * dg(['gamma',-1])
            if self.get_id('gamma') in arguments:
                du(['gamma',-1])[:] = 0.0
   
class Sys_Temp(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        Temp_IC = self.kwargs['Temp_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['Temp',-1], size=numPts, val=Temp_IC,
                               lower=0.001)
        self._declare_argument(['h',-1], indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p(['h',-1])
        #print  " current h-temp: ", h
        Temp = mission.get_temp(self.numElem, h)
        #print " current Temp: ", Temp
        u(['Temp',-1])[:] = Temp
    '''
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    '''
    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        numPts = self.numElem+1
        dTemp_dh = numpy.zeros(numPts)
        h = p(['h',-1])
        
        dTemp_dh = mission.get_temp_d(self.numElem, h)
        
        if self.mode == 'fwd':
            dg(['Temp',-1])[:] = 0.0
            if self.get_id('h') in arguments:
                dg(['Temp',-1])[:] = dTemp_dh * dp(['h',-1])
        if self.mode == 'rev':
            if self.get_id('h') in arguments:
                dp(['h',-1])[:] = dTemp_dh * dg(['Temp',-1])
            if self.get_id('Temp') in arguments:
                du(['Temp',-1])[:] = 0.0
    
class Sys_rho(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        rho_IC = self.kwargs['rho_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['rho',-1], size=numPts, val=rho_IC,
                               lower=0.001)
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument(['Temp',-1], indices=iPts)
        
    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        g = p(['g',0])
        Temp = p(['Temp',-1])
        rho = mission.get_rho(self.numElem, g, Temp)

        #print "current rho: ", rho
        u(['rho',-1])[:] = rho
    '''
    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)
    '''
    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        g = p(['g',0])
        Temp = p(['Temp',-1])
        numPts = self.numElem+1
        drho_dTemp = numpy.zeros(numPts)

        drho_dTemp = mission.get_rho_d(self.numElem, g, Temp)

        if self.mode == 'fwd':
            dg(['rho',-1])[:] = 0.0
            if self.get_id('Temp') in arguments:
                dg(['rho',-1])[:] += drho_dTemp * dp(['Temp',-1])
        if self.mode == 'rev':
            if self.get_id('Temp') in arguments:
                dp(['Temp',-1])[:] = drho_dTemp * dg(['rho',-1])
            if self.get_id('rho') in arguments:
                du(['rho',-1])[:] = 0.0
    
'''
class Sys_hDepInit(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        numPts = self.numElem+1
        
        self._declare_variable(['SFC',-1],size=numPts)
        self._declare_variable(['gamma',-1],size=numPts)
        self._declare_variable(['Temp',-1],size=numPts)
        self._declare_variable(['rho',-1],size=numPts)
        self._declare_argument(['h',-1],indices=[range(int(numPts))])
        self._declare_argument(['SFCSL',0],indices=[0])
        self._declare_argument(['g',0],indices=[0])
        
    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        h = p(['h',-1])
        SFCSL = p(['SFCSL',0])
        g = p(['g',0])
        SFC = mission.get_sfc(numElem, SFCSL, h)
        gamma = mission.get_gamma(numElem, h, x)
        Temp = mission.get_tmp(numElem, h)
        rho = mission.get_rho(numElem, g, Temp)
        u(['SFC',-1])[:] = SFC
        u(['gamma',-1])[:] = gamma
        u(['Temp',-1])[:] = Temp
        u(['rho',-1])[:] = rho

    def apply_dGdp(self, mode, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        numPts = self.numElem+1
        dSFC_dh = numpy.zeros(numPts)
        dgamma_dh1 = numpy.zeros(numPts)
        dgamma_dh2 = numpy.zeros(numPts)
        dgamma_dh3 = numpy.zeros(numPts)
        dgamma_dh4 = numpy.zeros(numPts)
        dgamma_dh5 = numpy.zeros(numPts)
        dgamma_dx = numpy.zeros(numPts)
        dTemp_dh = numpy.zeros(numPts)
        drho_dTemp = numpy.zeros(numPts)
        h = p(['h',-1])
        x = p(['x',0])

        [dSFC_dh] = mission.get_sfc_d(numElem)
        [dgamma_dh1, dgamma_dh2, dgamma_dh3, dgamma_dh4, dgamma_dh5,
         dgamma_dx] = mission.get_gamma_d(numElem, h, x)
        [dTemp_dh] = mission.get_tmp_d(numElem)
        [drho_dTemp] = mission.get_rho_d(numElem, g)

        if mode == 'fwd':
            dg(['SFC',-1])[:] = 0.0
            dg(['gamma',-1])[:] = 0.0
            dg(['Temp',-1])[:] = 0.0
            dg(['rho',-1])[:] = 0.0
            if self._get_ID('h') in arguments:
                dg(['SFC',-1])[:] += dSFC_dh * dp(['h',-1])
                dg(['gamma',-1])[:] += dgamma_dh1*dp(['h',-1]) + \
                    dgamma_dh2*dp(['h',-1]) + \
                    dgamma_dh3*dp(['h',-1]) + \
                    dgamma_dh4*dp(['h',-1]) + \
                    dgamma_dh5*dp(['h',-1])
                dg(['Temp',-1])[:] += dTemp_dh * dp(['h',-1])
                dg(['rho',-1])[:] += drho_dTemp*dTemp_dh * dp(['h',-1])
            if self._get_ID('x') in arguments:
                dg(['gamma',-1])[:] += dgamma_dx * dp(['h',-1])
        
        if mode == 'rev':
            dp(['h',-1])[:] = 0.0
            dp(['x',-1])[:] = 0.0
            if self._get_ID('h') in arguments:
                dp(['h',-1])[:] += dSFC_dh * dg(['SFC',-1])
                dp(['h',-1])[:] += [dgamma_dh1,dgamma_dh2,dgamma_dh3,
                                    dgamma_dh4,dgamma_dh5]*dg(['gamma',-1])
                dp(['h',-1])[:] += dTemp_dh * dg(['Temp',-1])
                dp(['h',-1])[:] += drho_dTemp*dTemp_dh * dg(['rho',-1])
            if self._get_ID('x') in arguments:
                dp(['h',-1])[:] += dgamma_dx * dg(['gamma',-1])
        
'''
class Sys_v(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        v_IC = self.kwargs['v_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        ##print "v_IC: ", v_IC
        self._declare_variable(['v',-1], size=numPts, val=v_IC,
                               lower=0.0)
        self._declare_argument(['Temp',-1], indices=iPts)
        self._declare_argument(['v_ends',-1], indices=[0,1])
        self._declare_argument(['M_ends',-1], indices=[0,1])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        v_ends = p(['v_ends',-1])
        M_ends = p(['M_ends',-1])
        Temp = p(['Temp',-1])
        v_copy = copy.copy(v_ends)
        
        if v_copy[0] == -1:
            v_copy[0] = M_ends[0]*numpy.sqrt(1.4*288*Temp[0])
        if v_copy[1] == -1:
            v_copy[1] = M_ends[1]*numpy.sqrt(1.4*288*Temp[self.numElem])

        v = mission.get_v(self.numElem, v_copy)
        #print "current v: ", v
        u(['v',-1])[:] = v
    '''
    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)
    '''
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
        M_ends = p(['M_ends',-1])
        Temp = p(['Temp',-1])
        v_ends = p(['v_ends',-1])

        [dv_dvEnds0, dv_dvEnds1] = mission.get_v_d(self.numElem, v_ends)

        if v_ends[0] == -1:
            dv_dMEnds0 = dv_dvEnds0*numpy.sqrt(1.4*288*Temp[0])
            dv_dvEnds0 = numpy.zeros(numPts)
        if v_ends[1] == -1:
            dv_dMEnds1 = dv_dvEnds1*numpy.sqrt(1.4*288*Temp[numPts-1])
            dv_dvEnds1 = numpy.zeros(numPts)

        if self.mode == 'fwd':
            dg(['v',-1])[:] = 0.0
            if self.get_id('v_ends') in arguments:
                dg(['v',-1])[:] += dv_dvEnds0 * dp(['v_ends',-1])[0] + \
                    dv_dvEnds1 * dp(['v_ends',-1])[1]
            if self.get_id('M_ends') in arguments:
                dg(['v',-1])[:] += dv_dMEnds0 * dp(['M_ends',-1])[0] + \
                    dv_dMEnds1 * dp(['M_ends',-1])[1]
        if self.mode == 'rev':
            if self.get_id('v_ends') in arguments:
                dp(['v_ends',-1])[:] = [dv_dvEnds0, dv_dvEnds1]*dg(['v',-1])
            if self.get_id('M_ends') in arguments:
                dp(['M_ends',-1])[:] = [dv_dMEnds0, dv_dMEnds1]*dg(['v',-1])
            if self.get_id('v') in arguments:
                dp(['v',-1])[:] = 0.0

class Sys_CT(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        t_IC = self.kwargs['t_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['CT',-1], size=numPts, val=t_IC)
        self._declare_argument(['h',-1], indices=iPts)
        self._declare_argument(['tau_init',-1], indices=[0])
        self._declare_argument(['rho',-1], indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['v',-1], indices=iPts)
        self._declare_argument(['cThrustSL',0], indices=[0])

    def apply_G(self):
        climb = self.kwargs['climb']
        if climb != -1:
            numElem = self.numElem
            p = self.vec['p']
            u = self.vec['u']

            tauInit = p(['tau_init',-1])
            rho = p(['rho',-1])
            v = p(['v',-1])
            S = p(['S',0])
            cThrustSL = p(['cThrustSL',0])
            h = p(['h',-1])
            CT = mission.get_ct_init(numElem, tauInit, S, cThrustSL, rho, v, h)

            #print "current CT: ", CT
            u(['CT',-1])[:] = CT

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dGdp(self, arguments):
        p = self.vec['p']
        dp = self.vec['dp']
        dg = self.vec['dg']
        du = self.vec['du']
        
        tau = p(['tau_init',-1])
        S = p(['S',0])
        cThrustSL = p(['cThrustSL',0])
        rho = p(['rho',-1])
        v = p(['v',-1])
        h = p(['h',-1])

        [dCTdcThrustSL, dCTdH, dCTdRho, dCTdV, dCTdS,
         dCTdTau] = mission.get_ct_init_d(self.numElem, tau, S, cThrustSL,
                                          rho, v, h)
        
        if self.mode == 'fwd':
            dg(['CT',-1])[:] = 0.0
            if self.get_id('cThrustSL') in arguments:
                dg(['CT',-1])[:] += dCTdcThrustSL * dp(['cThrustSL',0])
            if self.get_id('h') in arguments:
                dg(['CT',-1])[:] += dCTdH * dp(['h',-1])
            if self.get_id('rho') in arguments:
                dg(['CT',-1])[:] += dCTdRho * dp(['rho',-1])
            if self.get_id('v') in arguments:
                dg(['CT',-1])[:] += dCTdV * dp(['v',-1])
            if self.get_id('S') in arguments:
                dg(['CT',-1])[:] += dCTdS * dp(['S',0])
            if self.get_id('tau_init') in arguments:
                dg(['CT',-1])[:] += dCTdTau * dp(['tau_init',0])
        if self.mode == 'rev':
            if self.get_id('cThrustSL') in arguments:
                dp(['cThrustSL',0])[:] = dCTdcThrustSL * dg(['CT',-1])
            if self.get_id('h') in arguments:
                dp(['h',-1])[:] = dCTdH * dg(['CT',-1])
            if self.get_id('rho') in arguments:
                dp(['rho',-1])[:] = dCTdRho * dg(['CT',-1])
            if self.get_id('v') in arguments:
                dp(['v',-1])[:] = dCTdV * dg(['CT',-1])
            if self.get_id('S') in arguments:
                dp(['S',0])[:] = dCTdS * dg(['CT',-1])
            if self.get_id('tau_init') in arguments:
                dp(['tau_init',0])[:] = dCTdTau * dg(['CT',-1])

class Sys_tau(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable(['tau',-1], size=numPts)
        self._declare_argument(['CT',-1], indices=iPts)
        self._declare_argument(['rho',-1], indices=iPts)
        self._declare_argument(['v',-1], indices=iPts)
        self._declare_argument(['h',-1], indices=iPts)
        self._declare_argument(['cThrustSL',0], indices=[0])
        self._declare_argument(['S',0], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        CT = p(['CT',-1])
        rho = p(['rho',-1])
        v = p(['v',-1])
        h = p(['h',-1])
        cThrustSL = p(['cThrustSL',0])
        S = p(['S',0])

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

        cThrustSL = p(['cThrustSL',0])
        S = p(['S',0])
        h = p(['h',-1])
        CT = p(['CT',-1])
        rho = p(['rho',-1])
        v = p(['v',-1])
        
        [dtdcThrustSL, dtdh, dtdCT, dtdRho, dtdV,
         dtdS] = mission.get_tau_d(self.numElem, cThrustSL, S, h, CT,
                                   rho, v)

        if self.mode == 'fwd':
            dg(['tau',-1])[:] = 0.0
            if self.get_id('cThrustSL') in arguments:
                dg(['tau',-1])[:] += dtdcThrustSL * dp(['cThrustSL',0])
            if self.get_id('h') in arguments:
                dg(['tau',-1])[:] += dtdh * dp(['h',-1])
            if self.get_id('CT') in arguments:
                dg(['tau',-1])[:] += dtdCT * dp(['CT',-1])
            if self.get_id('rho') in arguments:
                dg(['tau',-1])[:] += dtdRho * dp(['rho',-1])
            if self.get_id('v') in arguments:
                dg(['tau',-1])[:] += dtdV * dp(['v',-1])
            if self.get_id('S') in arguments:
                dg(['tau',-1])[:] += dtdS * dp(['S',0])
        if self.mode == 'rev':
            if self.get_id('cThrustSL') in arguments:
                dp(['cThrustSL',0])[:] = dtdcThrustSL * dg(['tau',-1])
            if self.get_id('h') in arguments:
                dp(['h',-1])[:] = dtdh * dg(['tau',-1])
            if self.get_id('CT') in arguments:
                dp(['CT',-1])[:] = dtdCT * dg(['tau',-1])
            if self.get_id('rho') in arguments:
                dp(['rho',-1])[:] = dtdRho * dg(['tau',-1])
            if self.get_id('v') in arguments:
                dp(['v',-1])[:] = dtdV * dg(['tau',-1])
            if self.get_id('S') in arguments:
                dp(['S',0])[:] = dtdS * dg(['tau',-1])
            if self.get_id('tau') in arguments:
                du(['tau',-1])[:] = 0.0

class Sys_CL(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CL_IC = self.kwargs['CL_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['CL',-1], size=numPts, val=CL_IC)
        self._declare_argument(['Wf',-1], indices=iPts)
        self._declare_argument(['gamma',-1], indices=iPts)
        self._declare_argument(['CT',-1], indices=iPts)
        self._declare_argument(['alpha',-1], indices=iPts)
        self._declare_argument(['rho',-1], indices=iPts)
        self._declare_argument(['v',-1], indices=iPts)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['Wac',0], indices=[0])
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        
    def apply_F(self):
        self.numInt = self.kwargs['numInt']

        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        Wf = p(['Wf',-1])
        gamma = p(['gamma',-1])
        CT = p(['CT',-1])
        alpha = p(['alpha',-1])
        rho = p(['rho',-1])
        v = p(['v',-1])
        S = p(['S',0])
        Wac = p(['Wac',0])
        g = 9.81
        CL = u(['CL',-1])
        x_ends = p(['x',0])
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        CLRes = mission.get_cl(self.numElem, self.numInt, Wac, S, g, 
                               x_int, v, rho, CL, Wf, gamma, CT, alpha)
        #print "current CL: ", CL
        #print "current CLRes: ", CLRes
        f(['CL',-1])[:] = CLRes

    def apply_dFdpu0(self, arguments):
        self._apply_dFdpu_FD(arguments)

    def apply_dFdpu(self, arguments):
        p = self.vec['p']
        u = self.vec['u']
        dp = self.vec['dp']
        du = self.vec['du']

        Wac = p(['Wac',0])
        S = p(['S',0])
        x_ends = p(['x_ends',0])
        v = p(['v',-1])
        rho = p(['rho',-1])
        CL = u(['CL',-1])
        Wf = p(['Wf',-1])
        gamma = p(['gamma',-1])
        CT = p(['CT',-1])
        alpha = p(['alpha',-1])

class Sys_alpha(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        a_IC = self.kwargs['a_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['alpha',-1], size=numPts, val=a_IC)
        self._declare_argument(['CL',-1], indices=iPts)
        self._declare_argument(['eta',-1], indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']

        CL = p(['CL',-1])
        eta = p(['eta',-1])
        alpha = mission.get_alpha(self.numElem, CL, eta)
        #print "current alpha: ", alpha
        u(['alpha',-1])[:] = alpha

    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)

class Sys_CD(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CD_IC = self.kwargs['CD_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable(['CD',-1], size=numPts, val=CD_IC)
        self._declare_argument(['CL',-1], indices=iPts)
        self._declare_argument(['AR',0], indices=[0])
        self._declare_argument(['e',0], indices=[0])

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        CL = p(['CL',-1])
        AR = p(['AR',0])
        e = p(['e',0])
        CD = mission.get_cd(self.numElem, AR, e, CL)
        #print "current CD: ", CD
        u(['CD',-1])[:] = CD

    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)

class Sys_Wf(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        Wf_IC = self.kwargs['Wf_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable(['Wf',-1], size=numPts, val=Wf_IC,
                               u_scal=1e7, f_scal=1e7)
        self._declare_argument(['v',-1], indices=iPts)
        self._declare_argument(['gamma',-1], indices=iPts)
        self._declare_argument(['CT',-1], indices=iPts)
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        self._declare_argument(['g',0], indices=[0])
        self._declare_argument(['SFC',-1], indices=iPts)
        self._declare_argument(['rho',-1], indices=iPts)
        self._declare_argument(['WfSeg',0], indices=[self.copy+1])
        
    def apply_F(self):
        self._nln_init()
        self.numInt = self.kwargs['numInt']
        
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        x = p(['x',0])
        g = p(['g',0])
        v = p(['v',-1])
        gamma = p(['gamma',-1])
        CT = p(['CT',-1])
        SFC = p(['SFC',-1])
        rho = p(['rho',-1])
        WfIn = u(['Wf',-1])
        WfSeg = p(['WfSeg',0])
        xInt = numpy.linspace(x[0], x[1], self.numElem+1)*100.0

        Wf = mission.get_wf(self.numInt, self.numElem, xInt, v, gamma,
                            CT, SFC, rho, WfIn, g, WfSeg)
        
        #print "current Wf: ", WfIn
        #print "current WfRes: ", Wf-WfIn
        f(['Wf',-1])[:] = Wf-WfIn
        self._nln_final()

    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)

class Sys_CM(ImplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        CM_IC = self.kwargs['CM_IC']
        numPts = self.numElem+1
        iPts = range(numPts)

        self._declare_variable(['CM',-1], size=numPts, val=CM_IC)
        self._declare_argument(['S',0], indices=[0])
        self._declare_argument(['chord',0], indices=[0])
        self._declare_argument(['x',0], indices=[self.copy,self.copy+1])
        self._declare_argument(['v',-1], indices=iPts)
        self._declare_argument(['rho',-1], indices=iPts)

    def apply_F(self):
        self.numInt = self.kwargs['numInt']
        
        p = self.vec['p']
        u = self.vec['u']
        f = self.vec['f']

        S = p(['S',0])
        chord = p(['chord',0])
        x_ends = p(['x',0])
        v = p(['v',-1])
        rho = p(['rho',-1])
        CM = u(['CM',-1])
        x_int = numpy.linspace(x_ends[0], x_ends[1], self.numElem+1)

        CMRes = mission.get_cm(self.numElem, self.numInt, S, chord,
                               x_int, v, rho, CM)
        #print "current CM: ", CM
        #print "current CMRes: ", CMRes

        f(['CM',-1])[:] = CMRes

    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)

class Sys_eta(ExplicitSystem):
    def _declare(self):
        self.numElem = self.kwargs['numElem']
        e_IC = self.kwargs['e_IC']
        numPts = self.numElem+1
        iPts = range(numPts)
        
        self._declare_variable(['eta',-1], size=numPts, val=e_IC)
        self._declare_argument(['CM',-1], indices=iPts)
        self._declare_argument(['alpha',-1], indices=iPts)

    def apply_G(self):
        p = self.vec['p']
        u = self.vec['u']
        
        CM = p(['CM',-1])
        alpha = p(['alpha',-1])
        eta = mission.get_eta(self.numElem, CM, alpha)
        #print "current eta: ", eta
        u(['eta',-1])[:] = eta

    def apply_dFdpu(self, arguments):
        self._apply_dFdpu_FD(arguments)

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

        self.hPts[self.crtPt] = h
        self.xPts[self.crtPt] = x
        if v == -1:
            self.MPts[self.crtPt] = M
            self.vPts[self.crtPt] = -1
            self.vMType.append('M')
        else:
            self.MPts[self.crtPt] = -1
            self.vPts[self.crtPt] = v
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
            self.h_IC = h_IC

        if a_IC == None:
            self.a_IC = numpy.ones(totalElem+1)*3.0*numpy.pi/180.0
            #self.a_IC = 3.0*numpy.pi/180.0
        else:
            self.a_IC = a_IC       

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
            self.CD_IC = numpy.ones(totalElem+1)*0.01
            #self.CD_IC = 0.01
        else:
            self.CD_IC = CD_IC

        if CM_IC == None:
            self.CM_IC = numpy.zeros(totalElem+1)
            #self.CM_IC = 0.0
        else:
            self.CM_IC = CM_IC

        if Wf_IC == None:
            self.Wf_IC = numpy.linspace(100000.0,self.Wf_final,totalElem+1)
            #self.Wf_IC = 100000.0
        else:
            self.Wf_IC = Wf_IC

        if x_IC == None:
            self.x_IC = numpy.linspace(0.0,self.range,totalElem+1)
            #self.x_IC = 0.0
        else:
            self.x_IC = x_IC
        
        self.SFC_IC = mission.get_sfc(totalElem, self.SFCSL, self.h_IC)
        # potential bug: must allow diff x spacing in diff segments
        self.gamma_IC = mission.get_gamma(totalElem, self.h_IC, self.x_IC)
        self.Temp_IC = mission.get_temp(totalElem, self.h_IC)
        self.rho_IC = mission.get_rho(totalElem, self.g, self.Temp_IC)
        
        v_ends = numpy.zeros(self.numSeg+1)
        self.v_IC = []
        for seg in xrange(self.numSeg+1):
            if self.vPts[seg] == -1:
                v_ends[seg] = self.MPts[seg]*numpy.sqrt(1.4*288*self.Temp_IC[seg])
            else:
                v_ends[seg] = self.vPts[seg]

        for seg in xrange(self.numSeg):
            self.v_IC.extend(mission.get_v(self.numElem[seg], [v_ends[seg], v_ends[seg+1]]))
            self.v_IC.pop()
        self.v_IC.append(v_ends[self.numSeg])

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
        for pt in xrange(self.crtPt-1):
            pt2 += int(self.numElem[pt])
            h_IC = self.h_IC[pt1:pt2+1]
            t_IC = self.t_IC[pt1:pt2+1]
            e_IC = self.e_IC[pt1:pt2+1]
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
            '''
            self.segments.append(SerialSystem('segment',pt,subsystems=[
                        SerialSystem('seg_param',pt,subsystems=[
                                IndVar('h_ends',pt,val=[self.hPts[pt],self.hPts[pt+1]],size=2),
                                IndVar('v_ends',pt,val=[self.vPts[pt],self.vPts[pt+1]],size=2),
                                IndVar('M_ends',pt,val=[self.MPts[pt],self.MPts[pt+1]],size=2),
                                IndVar('h_dot',pt,val=self.hDotPts[pt],size=1),
                                IndVar('tau_init',pt,val=self.tPts[pt],size=1),
                                IndVar('numElem',pt,val=ne[pt],size=1),
                                ], NL='NLN_GS'),
                        #SerialSystem('init',pt,subsystems=[
                                #Sys_h('h',pt,des=self.hDotPts[pt],
                                      #opt=self.opt[pt],h_IC=h_IC,
                                      #numInt=self.numInt,numElem=self.numElem[pt]),
                                #Sys_SFCInit('SFC',pt,numElem=self.numElem[pt]),
                                #Sys_gammaInit('gamma',pt,numElem=self.numElem[pt]),
                                #Sys_TempInit('Temp',pt,numElem=self.numElem[pt]),
                                #Sys_rhoInit('rho',pt,numElem=self.numElem[pt]),
                                #Sys_vInit('v',pt,numElem=self.numElem[pt]),
                                #Sys_hDepInit(pt,numElem=self.numElem[pt]),
                                #Sys_CTInit('CT',pt,numElem=self.numElem[pt],climb=self.tPts[pt],t_IC=t_IC),
                                #], NL='NLN_GS',),
                        SerialSystem('seg_analysis',pt,subsystems=[
                                #IndVar('Wf',pt,val=ones(ne[pt]+1)*10000.0),
                                #IndVar('alpha',pt,val=ones(ne[pt]+1)*3.0*numpy.pi/180.0),
                                #IndVar('CT',pt,val=ones(ne[pt]+1)*0.1),
                                IndVar('CD',pt,val=ones(ne[pt]+1)*0.1),
                                #IndVar('rho',pt,val=ls(1.225,0.5,ne+1)),
                                #IndVar('v',pt,val=ls(150.0,250.0,ne+1)),
                                IndVar('eta',pt,val=numpy.zeros(ne[pt]+1)),

                                Sys_CL('CL',pt,numElem=self.numElem,numInt=self.numInt,CL_IC=self.CL_IC),
                                Sys_alpha('alpha',pt,numElem=self.numElem,a_IC=self.a_IC),
                                #Sys_CD(pt,CD_IC=CD_IC),
                                Sys_h('h',pt,climb=self.tPts[pt],h_IC=self.h_IC,numElem=self.numElem[pt],numInt=self.numInt),
                                #Sys_hDep(pt,numElem=self.numElem[pt]),
                                Sys_SFC('SFC',pt,numElem=self.numElem[pt],SFC_IC=self.SFC_IC),
                                Sys_gamma('gamma',pt,numElem=self.numElem[pt],gamma_IC=self.gamma_IC),
                                Sys_Temp('Temp',pt,numElem=self.numElem[pt],Temp_IC=self.Temp_IC),
                                Sys_rho('rho',pt,numElem=self.numElem[pt],rho_IC=self.rho_IC),
                                Sys_v('v',pt,numElem=self.numElem[pt],v_IC=self.v_IC),
                                #Sys_v(pt,numElem=self.numElem[pt]),
                                Sys_CT('CT',pt,numElem=self.numElem[pt],climb=self.tPts[pt],t_IC=t_IC),
                                #Sys_tau(pt,numElem=self.numElem[pt]),
                                Sys_Wf('Wf',pt,numElem=self.numElem[pt],numInt=self.numInt,Wf_IC=Wf_IC),
                                #Sys_CM(pt,CM_IC=CM_IC),
                                #Sys_eta(pt,e_IC=e_IC),
                                #]),
                                ], NL='NEWTON')
                        ], NL="NLN_GS"))
                        '''
        self.segments.append(IndVar('x',val=numpy.linspace(0,100.0e3,len(self.segments)+1)))
        self.segments.reverse()
        self.segments.append(IndVar('WfSeg',val=numpy.linspace(100000.0,0.0,len(self.segments))))
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
                        ], NL='NLN_GS'),
                SerialSystem('segments',0,subsystems=self.segments,
                             #self.segments.reverse(),
                             #Sys_xRes(x_IC=self.x_IC),
                             NL='NLN_JC')#NL='NLN_JC')
                ], NL='NLN_GS').setup()
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
                        IndVar('x',val=numpy.linspace(0,100.0e3,len(self.segments)+1)),
                        IndVar('h_ends',pt,val=[self.hPts[pt],self.hPts[pt+1]],size=2),
                        IndVar('v_ends',pt,val=[self.vPts[pt],self.vPts[pt+1]],size=2),
                        IndVar('M_ends',pt,val=[self.MPts[pt],self.MPts[pt+1]],size=2),
                        IndVar('h_dot',pt,val=self.hDotPts[pt],size=1),
                        IndVar('tau_init',pt,val=self.tPts[pt],size=1),
                        IndVar('numElem',pt,val=ne[pt],size=1),
                        #IndVar('Wf',pt,val=ones(ne[pt]+1)*10000.0),
                        #IndVar('CD',pt,val=ones(ne[pt]+1)*0.33),
                        #IndVar('eta',pt,val=numpy.zeros(ne[pt]+1)),
                        ], NL="NLN_GS"),
                Sys_CL('CL',pt,numElem=self.numElem,numInt=self.numInt,CL_IC=self.CL_IC),
                Sys_alpha('alpha',pt,numElem=self.numElem,a_IC=self.a_IC),
                Sys_CD('CD',pt,numElem=self.numElem,CD_IC=self.CD_IC),
                Sys_h('h',pt,climb=self.tPts[pt],h_IC=self.h_IC,numElem=self.numElem[pt],numInt=self.numInt),
                Sys_SFC('SFC',pt,numElem=self.numElem[pt],SFC_IC=self.SFC_IC),
                Sys_gamma('gamma',pt,numElem=self.numElem[pt],gamma_IC=self.gamma_IC),
                Sys_Temp('Temp',pt,numElem=self.numElem[pt],Temp_IC=self.Temp_IC),
                Sys_rho('rho',pt,numElem=self.numElem[pt],rho_IC=self.rho_IC),
                Sys_v('v',pt,numElem=self.numElem[pt],v_IC=self.v_IC),
                Sys_CT('CT',pt,numElem=self.numElem[pt],climb=self.tPts[pt],t_IC=t_IC),
                Sys_tau('tau',pt,numElem=self.numElem[pt]),
                Sys_Wf('Wf',pt,numElem=self.numElem[pt],numInt=self.numInt,Wf_IC=Wf_IC),
                Sys_CM('CM',pt,numElem=self.numElem[pt],numInt=self.numInt,CM_IC=self.CM_IC),
                Sys_eta('eta',pt,numElem=self.numElem[pt],e_IC=self.e_IC),
                IndVar('WfSeg',val=numpy.linspace(100000.0,0.0,len(self.segments)))
                ], NL='NEWTON', LN_ilimit=10, NL_ilimit=100, NL_rtol=1e-8).setup()
                

        return self.mainMission

params = {
    'S': 427.8,
    'Wac': 185953*9.81,
    'cThrustSL': 1020000.0,
    'SFCSL': 8.951e-6,
    'chord': 8.15,
    'inertia': 4.1e7,
    'AR': 8.68,
    'e': 0.8,
    'g': 9.81,
    }

n = 1
missionProblem = Trajectory(n)
missionProblem.add_seg_point(0.0,0.0,v=150.0)
for i in xrange(n):
    missionProblem.add_seg_point(8000.0,500.0,M=0.82,tau=1.0,numElem=5)
#missionProblem.add_seg_point(11000.0,1000.0,M=0.82,tau=1.0,numElem=5)
#missionProblem.add_seg_point(11000.0,4000.0,M=0.82,hDot=0.0,numElem=20)
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


print 'compute', 0, problemPtr.vec['u'].array.shape[0]
problemPtr.compute(True).array
problemPtr.vec['u'].array[:] *= problemPtr.vec['u0'].array[:]
#print problemPtr.vec['u']
print problemPtr.vec['u0']

print "---------------------------------------"
for name in ['SFC', 'gamma', 'Temp', 'h', 'rho', 'CT', 'v', 'tau']:
    print 'check ' + name, problemPtr([name,0]).check_derivatives(problemPtr.variables.keys())

#print 'checking h', problemPtr(['h',0]).check_derivatives(([('CD',0)]))

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
