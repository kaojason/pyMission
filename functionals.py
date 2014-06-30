"""
INTENDED FOR MISSION ANALYSIS USE
This file contains the functional systems used for the optimization
problem. These include objective and constraint functions defined for
the trajectory optimization case
"""

# pylint: disable=E1101
from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
from framework import *
import numpy


class SysHi(ExplicitSystem):
    """ initial altitude point used for constraints """

    def _declare(self):
        """ owned variable: h_i (initial altitude point)
            dependencies: h (altitude points)
        """

        self._declare_variable('h_i')
        self._declare_argument('h', indices=[0])

    def apply_G(self):
        """ assign system to the initial altitude point """

        alt_i = self.vec['u']('h_i')
        alt = self.vec['p']('h')

        alt_i[0] = alt[0]

    def apply_dGdp(self, args):
        """ derivative of this is same as intial altitude point """

        dalt_i = self.vec['dg']('h_i')
        dalt = self.vec['dp']('h')

        if self.mode == 'fwd':
            dalt_i[0] = 0.0
            if self.get_id('h') in args:
                dalt_i[0] += dalt[0]
        if self.mode == 'rev':
            dalt[0] = 0.0
            if self.get_id('h') in args:
                dalt[0] += dalt_i[0]

class SysHf(ExplicitSystem):
    """ final altitude point used for constraints """

    def _declare(self):
        """ owned variable: h_f (final altitude point)
            dependencies: h (altitude points)
        """

        num_elem = self.kwargs['num_elem']
        self._declare_variable('h_f')
        self._declare_argument('h', indices=[num_elem])

    def apply_G(self):
        """ assign system to the final altitude point """

        alt_f = self.vec['u']('h_f')
        alt = self.vec['p']('h')

        alt_f[0] = alt[0]

    def apply_dGdp(self, args):
        """ derivative of this is same as final altitude point """

        dalt_f = self.vec['dg']('h_f')
        dalt = self.vec['dp']('h')

        if self.mode == 'fwd':
            dalt_f[0] = 0.0
            if self.get_id('h') in args:
                dalt_f[0] += dalt[0]
        if self.mode == 'rev':
            dalt[0] = 0.0
            if self.get_id('h') in args:
                dalt[0] += dalt_f[0]

class SysTmin(ExplicitSystem):
    """ KS constraint function for minimum throttle """

    def _declare(self):
        """ owned variable: Tmin (minimum thrust constraint)
            dependencies: tau (throttle setting)
        """

        num_elem = self.kwargs['num_elem']

        self._declare_variable('Tmin')
        self._declare_argument('tau', indices=range(0, num_elem+1))
        self.min = 0.01
        self.rho = 30

    def apply_G(self):
        """ compute the KS function of minimum throttle """

        tmin = self.vec['u']('Tmin')
        tau = self.vec['p']('tau')

        fmax = numpy.max(self.min - tau)
        tmin[0] = fmax + 1/self.rho * \
            numpy.log(numpy.sum(numpy.exp(self.rho*(self.min - tau - fmax))))

    def apply_dGdp(self, args):
        """ compute min throttle KS function derivatives wrt throttle """

        tau = self.vec['p']('tau')

        dtmin = self.vec['dg']('Tmin')
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
            dtmin[0] = 0.0
            if self.get_id('tau') in args:
                dtmin[0] += numpy.sum(deriv * dtau[:])
        if self.mode == 'rev':
            dtau[:] = 0.0
            if self.get_id('tau') in args:
                dtau[:] += deriv * dtmin[0]

class SysTmax(ExplicitSystem):
    """ KS constraint function for maximum throttle """

    def _declare(self):
        """ owned variable: Tmax (maximum thrust constraint)
            dependencies: tau (throttle setting)
        """

        num_elem = self.kwargs['num_elem']
        self._declare_variable('Tmax')
        self._declare_argument('tau', indices=range(0, num_elem+1))
        self.max = 1.0
        self.rho = 30

    def apply_G(self):
        """ compute KS function for max throttle setting """

        tmax = self.vec['u']('Tmax')
        tau = self.vec['p']('tau')

        fmax = numpy.max(tau - self.max)
        tmax[0] = fmax + 1/self.rho * \
            numpy.log(numpy.sum(numpy.exp(self.rho*(tau - self.max - fmax))))

    def apply_dGdp(self, args):
        """ compute max throttle KS function derivatives wrt throttle """

        tau = self.vec['p']('tau')

        dtmax = self.vec['dg']('Tmax')
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
            dtmax[0] = 0.0
            if self.get_id('tau') in args:
                dtmax[0] += numpy.sum(deriv * dtau[:])
        if self.mode == 'rev':
            dtau[:] = 0.0
            if self.get_id('tau') in args:
                dtau[:] += deriv * dtmax[0]

class SysFuelObj(ExplicitSystem):
    """ objective function used for the optimization problem """

    def _declare(self):
        """ owned variable: objective fual burn (initial fuel carried)
            dependencies: weight of fuel (fuel_w)
        """

        self._declare_variable('wf_obj', size=1)
        self._declare_argument('fuel_w', indices=[0])

    def apply_G(self):
        """ set objective fuel weight to initial fuel carried (required for
            mission
        """

        p = self.vec['p']
        u = self.vec['u']
        Wf = p('fuel_w')

        u('wf_obj')[0] = Wf[0]

    def apply_dGdp(self, arguments):
        """ compute objective derivatives (equal to initial fuel weight
            derivative
        """

        dp = self.vec['dp']
        dg = self.vec['dg']
        
        if self.mode == 'fwd':
            dg('wf_obj')[0] = 0.0
            if self.get_id('fuel_w') in arguments:
                dg('wf_obj')[0] += dp('fuel_w')[0]
        if self.mode == 'rev':
            dp('fuel_w')[0] = 0.0
            if self.get_id('fuel_w') in arguments:
                dp('fuel_w')[0] += dg('wf_obj')[0]
