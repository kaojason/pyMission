from __future__ import division
from framework import *
import numpy
import time
import MBI, scipy.sparse

def setup_surrogate(surr_file):

    raw = numpy.loadtxt(surr_file+'_inputs.dat')
    [CL, CD, CM] = numpy.loadtxt(surr_file+'_outputs.dat')

    M_num, a_num, h_num, e_num = raw[:4].astype(int)
    M_surr = raw[4:4 + M_num]
    a_surr = raw[4 + M_num:4 + M_num + a_num]
    h_surr = raw[4 + M_num + a_num:4 + M_num + a_num + h_num]
    e_surr = raw[4 + M_num + a_num + h_num:]

    mbi_CL = numpy.zeros((M_num, a_num, h_num, e_num))
    mbi_CD = numpy.zeros((M_num, a_num, h_num, e_num))
    mbi_CM = numpy.zeros((M_num, a_num, h_num, e_num))
    
    count = 0
    for i in xrange(M_num):
        for j in xrange(a_num):
            for k in xrange(h_num):
                for l in xrange(e_num):
                    mbi_CL[i][j][k][l] = CL[count]
                    mbi_CD[i][j][k][l] = CD[count]
                    mbi_CM[i][j][k][l] = CM[count]
                    count += 1

    CL_arr = MBI.MBI(mbi_CL, [M_surr, a_surr, h_surr, e_surr],
                     [M_num, a_num, h_num, e_num], [4, 4, 4, 4])
    CD_arr = MBI.MBI(mbi_CD, [M_surr, a_surr, h_surr, e_surr],
                     [M_num, a_num, h_num, e_num], [4, 4, 4, 4])
    CM_arr = MBI.MBI(mbi_CM, [M_surr, a_surr, h_surr, e_surr],
                     [M_num, a_num, h_num, e_num], [4, 4, 4, 4])

    nums = {
        'M': M_num,
        'a': a_num,
        'h': h_num,
        'e': e_num,
        }
        
    return [CL_arr, CD_arr, CM_arr, nums]


class SysTripanCLSurrogate(ImplicitSystem):

    def _declare(self):

        self.num_elem = self.kwargs['num_elem']
        self.num = self.kwargs['num']
        self.CL_arr = self.kwargs['CL']
        ind_pts = range(self.num_elem + 1)

        self._declare_variable('alpha', size=self.num_elem+1)
        self._declare_argument('M', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)
        self._declare_argument('eta', indices=ind_pts)
        self._declare_argument('CL_tar', indices=ind_pts)

        self.J_CL = [None for i in range(4)]

    def apply_F(self):

        pvec = self.vec['p']
        uvec = self.vec['u']
        fvec = self.vec['f']

        M_num = self.num['M']
        a_num = self.num['a']
        h_num = self.num['h']
        e_num = self.num['e']

        Mach = pvec('M')
        alpha = uvec('alpha') * 180 / numpy.pi * 1e-1
        alt = pvec('h') * 3.28 * 1e3
        eta = pvec('eta') * 180 / numpy.pi * 1e-1
        CL_tar = pvec('CL_tar')
        alpha_res = fvec('alpha')

        inputs = numpy.zeros((self.num_elem + 1, 4))
        inputs[:, 0] = Mach
        inputs[:, 1] = alpha
        inputs[:, 2] = alt
        inputs[:, 3] = eta

        CL_temp = self.CL_arr.evaluate(inputs)

        CL = numpy.zeros(self.num_elem+1)
        for index in xrange(self.num_elem + 1):
            CL[index] = CL_temp[index, 0]

        flaps = Mach <= 0.4
        flaps = flaps * 5*(0.4-Mach)

        alpha_res[:] = (CL + flaps) - CL_tar

    def linearize(self):

        pvec = self.vec['p']
        uvec = self.vec['u']

        M_num = self.num['M']
        a_num = self.num['a']
        h_num = self.num['h']
        e_num = self.num['e']

        Mach = pvec('M')
        alpha = uvec('alpha') * 180 / numpy.pi * 1e-1
        alt = pvec('h') * 3.28 * 1e3
        eta = pvec('eta') * 180 / numpy.pi * 1e-1

        inputs = numpy.zeros((self.num_elem + 1, 4))
        inputs[:, 0] = Mach
        inputs[:, 1] = alpha
        inputs[:, 2] = alt
        inputs[:, 3] = eta

        for index in xrange(4):
            self.J_CL[index] = self.CL_arr.evaluate(inputs,
                                                    1+index, 0)[:, 0]

    def apply_dFdpu(self, args):

        dpvec = self.vec['dp']
        duvec = self.vec['du']
        dfvec = self.vec['df']
        pvec = self.vec['p']

        Mach = pvec('M')
        flaps = Mach <= 0.4
        flaps = flaps * (-5)

        dMach = dpvec('M')
        dalpha = duvec('alpha')
        dalt = dpvec('h')
        deta = dpvec('eta')
        dCL = dpvec('CL_tar')

        dres = dfvec('alpha')

        if self.mode == 'fwd':
            dres[:] = 0.0
            if self.get_id('M') in args:
                dres[:] += (self.J_CL[0]+flaps) * dMach
            if self.get_id('alpha') in args:
                dres[:] += self.J_CL[1] * dalpha * 180 / numpy.pi * 1e-1
            if self.get_id('h') in args:
                dres[:] += self.J_CL[2] * dalt * 3.28 * 1e3
            if self.get_id('eta') in args:
                dres[:] += self.J_CL[3] * deta * 180 / numpy.pi * 1e-1
            if self.get_id('CL_tar') in args:
                dres[:] -= dCL
        elif self.mode == 'rev':
            dMach[:] = 0.0
            dalpha[:] = 0.0
            dalt[:] = 0.0
            deta[:] = 0.0
            dCL[:] = 0.0
            if self.get_id('M') in args:
                dMach[:] += (self.J_CL[0]+flaps) * dres
            if self.get_id('alpha') in args:
                dalpha[:] += self.J_CL[1] * dres * 180 / numpy.pi * 1e-1
            if self.get_id('h') in args:
                dalt[:] += self.J_CL[2] * dres * 3.28 * 1e3
            if self.get_id('eta') in args:
                deta[:] += self.J_CL[3] * dres * 180 / numpy.pi * 1e-1
            if self.get_id('CL_tar') in args:
                dCL[:] -= dres

class SysTripanCDSurrogate(ExplicitSystem):

    def _declare(self):

        self.num_elem = self.kwargs['num_elem']
        self.num = self.kwargs['num']
        self.CD_arr = self.kwargs['CD']
        ind_pts = range(self.num_elem + 1)

        self._declare_variable('CD', size=self.num_elem+1)
        self._declare_argument('alpha', indices=ind_pts)
        self._declare_argument('M', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)
        self._declare_argument('eta', indices=ind_pts)

        self.J_CD = [None for i in range(4)]

    def apply_G(self):

        pvec = self.vec['p']
        uvec = self.vec['u']

        M_num = self.num['M']
        a_num = self.num['a']
        h_num = self.num['h']
        e_num = self.num['e']

        Mach = pvec('M')
        alpha = pvec('alpha') * 180 / numpy.pi * 1e-1
        alt = pvec('h') * 3.28 * 1e3
        eta = pvec('eta') * 180 / numpy.pi * 1e-1
        CD = uvec('CD')

        inputs = numpy.zeros((self.num_elem + 1, 4))
        inputs[:, 0] = Mach
        inputs[:, 1] = alpha
        inputs[:, 2] = alt
        inputs[:, 3] = eta

        CD_temp = self.CD_arr.evaluate(inputs)

        for index in xrange(self.num_elem + 1):
            CD[index] = CD_temp[index, 0] / 1e-1 + 0.015/1e-1

        flaps = Mach <= 0.4
        flaps = flaps * 0.25*(0.4-Mach)

        CD[:] += flaps/1e-1

    def linearize(self):

        pvec = self.vec['p']

        M_num = self.num['M']
        a_num = self.num['a']
        h_num = self.num['h']
        e_num = self.num['e']

        Mach = pvec('M')
        alpha = pvec('alpha') * 180 / numpy.pi * 1e-1
        alt = pvec('h') * 3.28 * 1e3
        eta = pvec('eta') * 180 / numpy.pi * 1e-1

        inputs = numpy.zeros((self.num_elem + 1, 4))
        inputs[:, 0] = Mach
        inputs[:, 1] = alpha
        inputs[:, 2] = alt
        inputs[:, 3] = eta

        for index in xrange(4):
            self.J_CD[index] = self.CD_arr.evaluate(inputs,
                                                    1+index, 0)[:,0]

    def apply_dGdp(self, args):

        dpvec = self.vec['dp']
        dgvec = self.vec['dg']
        pvec = self.vec['p']

        Mach = pvec('M')

        flaps = Mach <= 0.4
        flaps = flaps * (-0.25)

        dMach = dpvec('M')
        dalpha = dpvec('alpha')
        dalt = dpvec('h')
        deta = dpvec('eta')

        dCD = dgvec('CD')

        if self.mode == 'fwd':
            dCD[:] = 0.0
            if self.get_id('M') in args:
                dCD[:] += (self.J_CD[0]+flaps) * dMach / 1e-1
            if self.get_id('alpha') in args:
                dCD[:] += self.J_CD[1] * dalpha * 180 / numpy.pi
            if self.get_id('h') in args:
                dCD[:] += self.J_CD[2] * dalt * 3.28 * 1e3 / 1e-1
            if self.get_id('eta') in args:
                dCD[:] += self.J_CD[3] * deta * 180 / numpy.pi
        elif self.mode == 'rev':
            dMach[:] = 0.0
            dalpha[:] = 0.0
            dalt[:] = 0.0
            deta[:] = 0.0
            if self.get_id('M') in args:
                dMach[:] += (self.J_CD[0]+flaps) * dCD / 1e-1
            if self.get_id('alpha') in args:
                dalpha[:] += self.J_CD[1] * dCD * 180 / numpy.pi
            if self.get_id('h') in args:
                dalt[:] += self.J_CD[2] * dCD * 3.28 * 1e3 / 1e-1
            if self.get_id('eta') in args:
                deta[:] += self.J_CD[3] * dCD * 180 / numpy.pi

class SysTripanCMSurrogate(ImplicitSystem):

    def _declare(self):

        self.num_elem = self.kwargs['num_elem']
        self.num = self.kwargs['num']
        self.CM_arr = self.kwargs['CM']
        ind_pts = range(self.num_elem + 1)

        self._declare_variable('eta', size=self.num_elem+1)
        self._declare_argument('alpha', indices=ind_pts)
        self._declare_argument('M', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)

        self.J_CM = [None for i in range(4)]

    def apply_F(self):

        pvec = self.vec['p']
        uvec = self.vec['u']
        fvec = self.vec['f']

        M_num = self.num['M']
        a_num = self.num['a']
        h_num = self.num['h']
        e_num = self.num['e']

        Mach = pvec('M')
        alpha = pvec('alpha') * 180 / numpy.pi * 1e-1
        alt = pvec('h') * 3.28 * 1e3
        eta = uvec('eta') * 180 / numpy.pi * 1e-1
        res = fvec('eta')

        inputs = numpy.zeros((self.num_elem + 1, 4))
        inputs[:, 0] = Mach
        inputs[:, 1] = alpha
        inputs[:, 2] = alt
        inputs[:, 3] = eta

        CM_temp = self.CM_arr.evaluate(inputs)

        for index in xrange(self.num_elem + 1):
            res[index] = CM_temp[index, 0]

    def linearize(self):

        pvec = self.vec['p']
        uvec = self.vec['u']

        M_num = self.num['M']
        a_num = self.num['a']
        h_num = self.num['h']
        e_num = self.num['e']

        Mach = pvec('M')
        alpha = pvec('alpha') * 180 / numpy.pi * 1e-1
        alt = pvec('h') * 3.28 * 1e3
        eta = uvec('eta') * 180 / numpy.pi * 1e-1

        inputs = numpy.zeros((self.num_elem + 1, 4))
        inputs[:, 0] = Mach
        inputs[:, 1] = alpha
        inputs[:, 2] = alt
        inputs[:, 3] = eta

        for index in xrange(4):
            self.J_CM[index] = self.CM_arr.evaluate(inputs,
                                                    1+index, 0)[:,0]

    def apply_dFdpu(self, args):

        dpvec = self.vec['dp']
        duvec = self.vec['du']
        dfvec = self.vec['df']

        dMach = dpvec('M')
        dalpha = dpvec('alpha')
        dalt = dpvec('h')
        deta = duvec('eta')
        dCM = dfvec('eta')

        if self.mode == 'fwd':
            dCM[:] = 0.0
            if self.get_id('M') in args:
                dCM[:] += self.J_CM[0] * dMach
            if self.get_id('alpha') in args:
                dCM[:] += self.J_CM[1] * dalpha * 180 / numpy.pi * 1e-1
            if self.get_id('h') in args:
                dCM[:] += self.J_CM[2] * dalt * 3.28 * 1e3
            if self.get_id('eta') in args:
                dCM[:] += self.J_CM[3] * deta * 180 / numpy.pi * 1e-1
        elif self.mode == 'rev':
            dMach[:] = 0.0
            dalpha[:] = 0.0
            dalt[:] = 0.0
            deta[:] = 0.0
            if self.get_id('M') in args:
                dMach[:] += self.J_CM[0] * dCM
            if self.get_id('alpha') in args:
                dalpha[:] += self.J_CM[1] * dCM * 180 / numpy.pi * 1e-1
            if self.get_id('h') in args:
                dalt[:] += self.J_CM[2] * dCM * 3.28 * 1e3
            if self.get_id('eta') in args:
                deta[:] += self.J_CM[3] * dCM * 180 / numpy.pi * 1e-1

'''
class SysTripanSurrogate1(ExplicitSystem):

    def setupGlobalModel(self, Data, theta0, GE=False, bf_flag=False, corrfunction=0, scaling=True):
        # User needs to specify: Data and theta0
        # Data is a dictionary with the following keys:
        # * x = sample locations, size [Ns, Ndv]
        # * yfun = function values, size [Ns]
        # * Ns = number of samples
        # * lb = lower bound, size [Ndv]
        # * ub = upper bound, size [Ndv]
        # * [OPTIONAL] ygrad = derivative information, size [Ns, Ndv]
        
        # theta0: the initial values for the hyperparameters, size [Ndv]
        # just initialize to 10.0*numpy.ones(Ndv) 
        # (you can basically change it to any numbers, but I found that 10.0 is a safe number to start)
        
        # The following variables have default values assigned to them, you might want to pay attention to GE
        # GE: True for Gradient-enhanced kriging, False for ordinary kriging
        # bf_flag: True for universal kriging (with specified basis functions for the global model), False for ordinary kriging model (using a constant global model)
        # corrfunction: 0 for Gaussian correlation function, 1 for cubic spline
        # scaling: True when you scale the input variables to [0,1], False when you don't use scaling (using scaling is better)
        
        x = Data['x'].copy()
        yfun = Data['yfun'].copy()
        Ns = Data['Ns']
        lb = Data['lb']
        ub = Data['ub']

        if 'ygrad' in Data.keys():
            ygrad = Data['ygrad'].copy()
        else:
            ygrad = []

        Ndv = x.shape[1]

        if GE:
            y = numpy.append(yfun, numpy.reshape(ygrad.T, Ndv*Ns, 1))
        else:
            y = yfun

#        model = Kriging(theta0, GE=GE, bf_flag=bf_flag, corrfunction=corrfunction)
        model = RBF()
   
        if scaling:
            model.setSampleData(x, y, Lower=numpy.array(lb), Upper=numpy.array(ub))
        else:
            model.setSampleData(x, y)
        model.build()

        return model

    def _declare(self):
        self.num_elem = self.kwargs['num_elem']
        ind_pts = range(self.num_elem + 1)

        surr_size = 50

        self._declare_variable('CL', size=self.num_elem + 1)
        self._declare_variable('CD', size=self.num_elem + 1)
        self._declare_variable('Cm', size=self.num_elem + 1)
        self._declare_argument('alpha', indices=ind_pts)
        self._declare_argument('M', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)                                                  
        self._declare_argument('eta', indices=ind_pts)

        print 'Loading Tripan data'

        filename = 'AeroData/CRMTransport_4D_CLkriging.pkl'
        cPicklefile = open(filename, 'r')

        CLData = cPickle.load(cPicklefile)
        cPicklefile.close()
        
# Load data to build kriging model for CM
        
        filename = 'AeroData/CRMTransport_4D_CMkriging.pkl'
        cPicklefile = open(filename, 'r')
        
        CMData = cPickle.load(cPicklefile)
        cPicklefile.close()
        
# Load data to build kriging model for CD
        
        filename = 'AeroData/CRMTransport_4D_CDGEK.pkl'
        cPicklefile = open(filename, 'r')
        
        CDData = cPickle.load(cPicklefile)
        
        cPicklefile.close()

        Ndv = 4
        theta0 = 10.0*numpy.ones(Ndv)

        t0 = time.time()
        self.CLmodel = self.setupGlobalModel(CLData, theta0 = theta0, GE=False, bf_flag=False)
        t1 = time.time()
        print t1-t0
        self.CMmodel = self.setupGlobalModel(CMData, theta0 = theta0, GE=False, bf_flag=False)
        t2 = time.time()
        print t2-t1
        self.CDmodel = self.setupGlobalModel(CDData, theta0 = theta0, GE=False, bf_flag=False)
        t3 = time.time()
        print t3-t2

        self.inputs = numpy.zeros((self.num_elem + 1, 4))

        print 'Done loading Tripan data'

        M_surr = numpy.linspace(0, 1.2, surr_num)
        a_surr = numpy.linspace(-10,20, surr_num)
        h_surr = numpy.linspace(0, 50000, surr_num)
        e_surr = numpy.linspace(-20, 20, surr_num)
        surr_inp = [M_surr, a_surr, h_surr, e_surr]

        count = 0
        inputs = numpy.zeros(surr_num^4, 4)
        for entry in surr_inp:
            

    def apply_G(self):
        CL = self.vec['u']('CL')
        CD = self.vec['u']('CD')
        Cm = self.vec['u']('Cm')

        alpha = self.vec['p']('alpha') * 180 / numpy.pi * 1e-1
        M = self.vec['p']('M')
        h = self.vec['p']('h') * 0.3048 * 1e3
        eta = self.vec['p']('eta') * 180 / numpy.pi * 1e-1

        self.inputs[:, 0] = M
        self.inputs[:, 1] = alpha
        self.inputs[:, 2] = h
        self.inputs[:, 3] = eta

        for model in [self.CLmodel, self.CDmodel, self.CMmodel]:
            model.evaluate(self.inputs)

        CL[:] = self.CLmodel.yhat
        CD[:] = self.CDmodel.yhat / 1e-1
        Cm[:] = self.CMmodel.yhat

    def linearize(self):
        for model in [self.CLmodel, self.CDmodel, self.CMmodel]:
            model.computeGradient(self.inputs)

    def apply_dGdp(self, args):
        dCL = self.vec['dg']('CL')
        dCD = self.vec['dg']('CD')
        dCm = self.vec['dg']('Cm')

        dalpha = self.vec['dp']('alpha')
        dM = self.vec['dp']('M')
        dh = self.vec['dp']('h')
        deta = self.vec['dp']('eta')

        CLgrad = self.CLmodel.gradients
        CDgrad = self.CDmodel.gradients
        CMgrad = self.CMmodel.gradients

        if self.mode == 'fwd':
            dCL[:] = 0.0
            dCD[:] = 0.0
            dCm[:] = 0.0
            if self.get_id('alpha') in args:
                dCL[:] += CLgrad[:, 1] * dalpha * 180/numpy.pi * 1e-1
                dCD[:] += CDgrad[:, 1] * dalpha * 180/numpy.pi * 1e-1 / 1e-1
                dCm[:] += CMgrad[:, 1] * dalpha * 180/numpy.pi * 1e-1
            if self.get_id('M') in args:
                dCL[:] += CLgrad[:, 0] * dM
                dCD[:] += CDgrad[:, 0] * dM / 1e-1 
                dCm[:] += CMgrad[:, 0] * dM 
            if self.get_id('h') in args:
                dCL[:] += CLgrad[:, 2] * dh * 0.3028 * 1e3
                dCD[:] += CDgrad[:, 2] * dh * 0.3028 * 1e3 / 1e-1
                dCm[:] += CMgrad[:, 2] * dh * 0.3028 * 1e3
            if self.get_id('eta') in args:
                dCL[:] += CLgrad[:, 3] * deta * 180/numpy.pi * 1e-1
                dCD[:] += CDgrad[:, 3] * deta * 180/numpy.pi * 1e-1 / 1e-1
                dCm[:] += CMgrad[:, 3] * deta * 180/numpy.pi * 1e-1
        elif self.mode == 'rev':
            dalpha[:] = 0.0
            dM[:] = 0.0
            dh[:] = 0.0
            deta[:] = 0.0
            if self.get_id('alpha') in args:
                dalpha[:] += CLgrad[:, 1] * dCL * 180/numpy.pi * 1e-1
                dalpha[:] += CDgrad[:, 1] * dCD * 180/numpy.pi * 1e-1 / 1e-1
                dalpha[:] += CMgrad[:, 1] * dCm * 180/numpy.pi * 1e-1
            if self.get_id('M') in args:
                dM[:] += CLgrad[:, 0] * dCL
                dM[:] += CDgrad[:, 0] * dCD / 1e-1 
                dM[:] += CMgrad[:, 0] * dCm
            if self.get_id('h') in args:
                dh[:] += CLgrad[:, 2] * dCL * 0.3028 * 1e3
                dh[:] += CDgrad[:, 2] * dCD * 0.3028 * 1e3 / 1e-1
                dh[:] += CMgrad[:, 2] * dCm * 0.3028 * 1e3
            if self.get_id('eta') in args:
                deta[:] += CLgrad[:, 3] * dCL * 180/numpy.pi * 1e-1
                deta[:] += CDgrad[:, 3] * dCD * 180/numpy.pi * 1e-1 / 1e-1
                deta[:] += CMgrad[:, 3] * dCm * 180/numpy.pi * 1e-1



class SysTripanSurrogate0(ExplicitSystem):

    def setupGlobalModel(self, Data, theta0, GE=False, bf_flag=False, corrfunction=0, scaling=True):
        # User needs to specify: Data and theta0
        # Data is a dictionary with the following keys:
        # * x = sample locations, size [Ns, Ndv]
        # * yfun = function values, size [Ns]
        # * Ns = number of samples
        # * lb = lower bound, size [Ndv]
        # * ub = upper bound, size [Ndv]
        # * [OPTIONAL] ygrad = derivative information, size [Ns, Ndv]
        
        # theta0: the initial values for the hyperparameters, size [Ndv]
        # just initialize to 10.0*numpy.ones(Ndv) 
        # (you can basically change it to any numbers, but I found that 10.0 is a safe number to start)
        
        # The following variables have default values assigned to them, you might want to pay attention to GE
        # GE: True for Gradient-enhanced kriging, False for ordinary kriging
        # bf_flag: True for universal kriging (with specified basis functions for the global model), False for ordinary kriging model (using a constant global model)
        # corrfunction: 0 for Gaussian correlation function, 1 for cubic spline
        # scaling: True when you scale the input variables to [0,1], False when you don't use scaling (using scaling is better)
        
        x = Data['x'].copy()
        yfun = Data['yfun'].copy()
        Ns = Data['Ns']
        lb = Data['lb']
        ub = Data['ub']

        if 'ygrad' in Data.keys():
            ygrad = Data['ygrad'].copy()
        else:
            ygrad = []

        Ndv = x.shape[1]

        #if GE:
        #    y = numpy.append(yfun, numpy.reshape(ygrad.T, Ndv*Ns, 1))
        #else:
        #    y = yfun

        #model = Kriging(theta0, GE=GE, bf_flag=bf_flag, corrfunction=corrfunction)
        model = RBF()
        model.setSampleData(x, yfun)
    
        #if scaling:
        #    model.setSampleData(x, y, Lower=numpy.array(lb), Upper=numpy.array(ub))
        #else:
        #    model.setSampleData(x, y)
        model.build()

        return model

    def _declare(self):
        self.num_elem = self.kwargs['num_elem']
        ind_pts = range(self.num_elem + 1)

        self._declare_variable('CL', size=self.num_elem + 1)
        self._declare_variable('CD', size=self.num_elem + 1)
        self._declare_variable('Cm', size=self.num_elem + 1)
        self._declare_argument('alpha', indices=ind_pts)
        self._declare_argument('M', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)                                                  
        self._declare_argument('eta', indices=ind_pts)

        print 'Loading Tripan data'

        filename = 'AeroData/CRMTransport_4D_CLkriging.pkl'
        cPicklefile = open(filename, 'r')

        CLData = cPickle.load(cPicklefile)
        cPicklefile.close()
        
# Load data to build kriging model for CM
        
        filename = 'AeroData/CRMTransport_4D_CMkriging.pkl'
        cPicklefile = open(filename, 'r')
        
        CMData = cPickle.load(cPicklefile)
        cPicklefile.close()
        
# Load data to build kriging model for CD
        
        filename = 'AeroData/CRMTransport_4D_CDGEK.pkl'
        cPicklefile = open(filename, 'r')
        
        CDData = cPickle.load(cPicklefile)
        
        cPicklefile.close()

        Ndv = 4
        theta0 = 10.0*numpy.ones(Ndv)

        t0 = time.time()
        self.CLmodel = self.setupGlobalModel(CLData, theta0 = theta0, GE=False, bf_flag=False); print '1'
        t1 = time.time()
        print t1-t0
        self.CMmodel = self.setupGlobalModel(CMData, theta0 = theta0, GE=False, bf_flag=False); print '2'
        t2 = time.time()
        print t2-t1
        self.CDmodel = self.setupGlobalModel(CDData, theta0 = theta0, GE=False, bf_flag=True); print '3'
        t3 = time.time()
        print t3-t2

        self.inputs = numpy.zeros((self.num_elem + 1, 4))

        print 'Done loading Tripan data'

    def apply_G(self):
        CL = self.vec['u']('CL')
        CD = self.vec['u']('CD')
        Cm = self.vec['u']('Cm')

        alpha = self.vec['p']('alpha') * 180 / numpy.pi * 1e-1
        M = self.vec['p']('M')
        h = self.vec['p']('h') * 0.3048 * 1e3
        eta = self.vec['p']('eta') * 180 / numpy.pi * 1e-1

        self.inputs[:, 0] = M
        self.inputs[:, 1] = alpha
        self.inputs[:, 2] = h
        self.inputs[:, 3] = eta

        for model in [self.CLmodel, self.CDmodel, self.CMmodel]:
            model.evaluate(self.inputs)

        CL[:] = self.CLmodel.yhat
        CD[:] = self.CDmodel.yhat / 1e-1
        Cm[:] = self.CMmodel.yhat

    def linearize(self):
        for model in [self.CLmodel, self.CDmodel, self.CMmodel]:
            model.computeGradient(self.inputs)

    def apply_dGdp(self, args):
        dCL = self.vec['dg']('CL')
        dCD = self.vec['dg']('CD')
        dCm = self.vec['dg']('Cm')

        dalpha = self.vec['dp']('alpha')
        dM = self.vec['dp']('M')
        dh = self.vec['dp']('h')
        deta = self.vec['dp']('eta')

        CLgrad = self.CLmodel.gradients
        CDgrad = self.CDmodel.gradients
        CMgrad = self.CMmodel.gradients

        if self.mode == 'fwd':
            dCL[:] = 0.0
            dCD[:] = 0.0
            dCm[:] = 0.0
            if self.get_id('alpha') in args:
                dCL[:] += CLgrad[:, 1] * dalpha * 180/numpy.pi * 1e-1
                dCD[:] += CDgrad[:, 1] * dalpha * 180/numpy.pi * 1e-1 / 1e-1
                dCm[:] += CMgrad[:, 1] * dalpha * 180/numpy.pi * 1e-1
            if self.get_id('M') in args:
                dCL[:] += CLgrad[:, 0] * dM
                dCD[:] += CDgrad[:, 0] * dM / 1e-1 
                dCm[:] += CMgrad[:, 0] * dM 
            if self.get_id('h') in args:
                dCL[:] += CLgrad[:, 2] * dh * 0.3028 * 1e3
                dCD[:] += CDgrad[:, 2] * dh * 0.3028 * 1e3 / 1e-1
                dCm[:] += CMgrad[:, 2] * dh * 0.3028 * 1e3
            if self.get_id('eta') in args:
                dCL[:] += CLgrad[:, 3] * deta * 180/numpy.pi * 1e-1
                dCD[:] += CDgrad[:, 3] * deta * 180/numpy.pi * 1e-1 / 1e-1
                dCm[:] += CMgrad[:, 3] * deta * 180/numpy.pi * 1e-1
        elif self.mode == 'rev':
            dalpha[:] = 0.0
            dM[:] = 0.0
            dh[:] = 0.0
            deta[:] = 0.0
            if self.get_id('alpha') in args:
                dalpha[:] += CLgrad[:, 1] * dCL * 180/numpy.pi * 1e-1
                dalpha[:] += CDgrad[:, 1] * dCD * 180/numpy.pi * 1e-1 / 1e-1
                dalpha[:] += CMgrad[:, 1] * dCm * 180/numpy.pi * 1e-1
            if self.get_id('M') in args:
                dM[:] += CLgrad[:, 0] * dCL
                dM[:] += CDgrad[:, 0] * dCD / 1e-1 
                dM[:] += CMgrad[:, 0] * dCm
            if self.get_id('h') in args:
                dh[:] += CLgrad[:, 2] * dCL * 0.3028 * 1e3
                dh[:] += CDgrad[:, 2] * dCD * 0.3028 * 1e3 / 1e-1
                dh[:] += CMgrad[:, 2] * dCm * 0.3028 * 1e3
            if self.get_id('eta') in args:
                deta[:] += CLgrad[:, 3] * dCL * 180/numpy.pi * 1e-1
                deta[:] += CDgrad[:, 3] * dCD * 180/numpy.pi * 1e-1 / 1e-1
                deta[:] += CMgrad[:, 3] * dCm * 180/numpy.pi * 1e-1
'''
