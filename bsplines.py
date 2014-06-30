"""
INTENDED FOR MISSION ANALYSIS USE
b-splines parameterization of altitude, x-distance, and Mach number.
These provide altitude and Mach number profile wrt x-distance.
Gamma (flight path angle) is also computed using the altitude
parameterization
"""

# pylint: disable=E1101
from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
from framework import *
import numpy
import MBI, scipy.sparse

class BSplineSystem(ExplicitSystem):
    """ class used to allow the setup of b-splines """

    def MBI_setup(self):
        """ generate jacobians for b-splines using MBI package """

        num_pts = self.num_pts
        num_pt = self.num_pt

        alt = numpy.linspace(0, 13, num_pts)
        x_dist = numpy.linspace(0, self.range, num_pts)/1e6

        arr = MBI.MBI(alt, [x_dist], [num_pt], [4])
        jac = arr.getJacobian(0, 0)
        jacd = arr.getJacobian(1, 0)

        c_arryx = self.range/1e6 * 0.5*\
            (1-numpy.cos(numpy.pi*numpy.linspace(0, 1, num_pt)))
        d_arryx = jacd.dot(c_arryx)*1e6

        lins = numpy.linspace(0, num_pts-1, num_pts).astype(int)
        diag = scipy.sparse.csc_matrix((1.0/d_arryx,
                                        (lins,lins)))
        jace = diag.dot(jacd)

        self.jac_h = jac
        self.jac_gamma = jace
        self.c_arryx = c_arryx

class SysXBspline(BSplineSystem):
    """ a b-spline parameterization of distance """

    def _declare(self):
        """ owned variable: x (with b-spline parameterization)
            dependencies: x_pt (control points)
        """

        self.num_pts = self.kwargs['num_elem']+1
        self.num_pt = self.kwargs['num_pt']
        self.range = self.kwargs['x_range']*1e6
        x_0 = self.kwargs['x_0']

        self._declare_variable('x', size=self.num_pts, val=x_0)
        self._declare_argument('x_pt', indices=range(self.num_pt))
        self.MBI_setup()

    def apply_G(self):
        """ compute x b-spline values with x control point values
            using pre-calculated MBI jacobian
        """

        x_dist = self.vec['u']('x')
        x_pt = self.vec['p']('x_pt')
        x_dist[:] = self.jac_h.dot(x_pt[:])

    def apply_dGdp(self, args):
        """ compute derivative of x b-spline values wrt x control points
            using pre-calculated MBI jacobian
        """

        dx_dist = self.vec['dg']('x')
        dx_pt = self.vec['dp']('x_pt')

        if self.mode == 'fwd':
            dx_dist[:] = 0.0
            if self.get_id('x_pt') in args:
                dx_dist[:] += self.jac_h.dot(dx_pt[:])
        if self.mode == 'rev':
            dx_pt[:] = 0.0
            if self.get_id('x_pt') in args:
                dx_pt[:] += self.jac_h.T.dot(dx_dist[:])

class SysHBspline(BSplineSystem):
    """ a b-spline parameterization of altitude """

    def _declare(self):
        """ owned variable: h (altitude with b-spline parameterization)
            dependencies: h_pt (control points)
        """

        self.num_pts = self.kwargs['num_elem']+1
        self.num_pt = self.kwargs['num_pt']
        self.range = self.kwargs['x_range']*1e6

        self._declare_variable('h', size=self.num_pts)
        self._declare_argument('h_pt', indices=range(self.num_pt))
        self.MBI_setup()

    def apply_G(self):
        """ compute h b-splines values using h control point values
            using pre-calculated MBI jacobian
        """

        alt = self.vec['u']('h')
        alt_pt = self.vec['p']('h_pt')
        alt[:] = self.jac_h.dot(alt_pt[:])

    def apply_dGdp(self, args):
        """ compute h b-spline derivatives wrt h control points
            using pre-calculated MBI jacobian
        """

        dalt = self.vec['dg']('h')
        dalt_pt = self.vec['dp']('h_pt')

        if self.mode == 'fwd':
            dalt[:] = 0.0
            if self.get_id('h_pt') in args:
                dalt[:] += self.jac_h.dot(dalt_pt[:])
        if self.mode == 'rev':
            dalt_pt[:] = 0.0
            if self.get_id('h_pt') in args:
                dalt_pt[:] += self.jac_h.T.dot(dalt[:])


class SysVBspline(BSplineSystem):
    """ a b-spline parameterization of velocity """

    def _declare(self):
        """ owned variable: v (speed with b-spline parameterization)
            dependencies: v_pt (control points)
        """

        self.num_pts = self.kwargs['num_elem']+1
        self.num_pt = self.kwargs['num_pt']
        self.range = self.kwargs['x_range']*1e6

        self._declare_variable('v', size=self.num_pts)
        self._declare_argument('v_pt', indices=range(self.num_pt))
        self.MBI_setup()

    def apply_G(self):
        """ compute v b-spline values using v control point values
            using pre-calculated MBI jacobian
        """

        speed = self.vec['u']('v')
        speed_pt = self.vec['p']('v_pt')
        speed[:] = self.jac_h.dot(speed_pt[:])

    def apply_dGdp(self, args):
        """ compute v b-spline derivatives wrt v control points
            using pre-calculated MBI jacobian
        """
        dspeed = self.vec['dg']('v')
        dspeed_pt = self.vec['dp']('v_pt')

        if self.mode == 'fwd':
            dspeed[:] = 0.0
            if self.get_id('v_pt') in args:
                dspeed[:] += self.jac_h.dot(dspeed_pt[:])
        if self.mode == 'rev':
            dspeed_pt[:] = 0.0
            if self.get_id('v_pt') in args:
                dspeed_pt[:] += self.jac_h.T.dot(dspeed[:])

class SysMBspline(BSplineSystem):
    """ a b-spline parameterization of Mach number """

    def _declare(self):
        """ owned variable: M (Mach # with b-spline parameterization)
            dependencies: M_pt (control points)
        """

        self.num_pts = self.kwargs['num_elem']+1
        self.num_pt = self.kwargs['num_pt']
        self.range = self.kwargs['x_range']*1e6

        self._declare_variable('M', size=self.num_pts)
        self._declare_argument('M_pt', indices=range(self.num_pt))
        self.MBI_setup()

    def apply_G(self):
        """ compute M b-spline values using M control point values
            using pre-calculated MBI jacobian
        """

        mach = self.vec['u']('M')
        mach_pt = self.vec['p']('M_pt')
        mach[:] = self.jac_h.dot(mach_pt[:])

    def apply_dGdp(self, args):
        """ compute M b-spline derivatives wrt M control points
            using pre-calculated MBI jacobian
        """
        dmach = self.vec['dg']('M')
        dmach_pt = self.vec['dp']('M_pt')

        if self.mode == 'fwd':
            dmach[:] = 0.0
            if self.get_id('M_pt') in args:
                dmach[:] += self.jac_h.dot(dmach_pt[:])
        if self.mode == 'rev':
            dmach_pt[:] = 0.0
            if self.get_id('M_pt') in args:
                dmach_pt[:] += self.jac_h.T.dot(dmach[:])

class SysGammaBspline(BSplineSystem):
    """ dh/dx obtained from b-spline parameterization of altitude """

    def _declare(self):
        """ owned variable: gamma (flight path angle w/ b-spline 
            parameterization
            dependencies: h_pt (altitude control points)
        """

        self.num_pts = self.kwargs['num_elem']+1
        self.num_pt = self.kwargs['num_pt']
        self.range = self.kwargs['x_range']*1e6

        self._declare_variable('gamma', size=self.num_pts)
        self._declare_argument('h_pt', indices=range(self.num_pt))
        self.MBI_setup()

    def apply_G(self):
        """ compute gamma b-spline values using gamma control point values
            using pre-calculated MBI jacobian
        """

        gamma = self.vec['u']('gamma')
        h_pt = self.vec['p']('h_pt')
        gamma[:] = self.jac_gamma.dot(h_pt[:]) * 1e3/1e-1

    def apply_dGdp(self, args):
        """ compute gamma b-spline derivatives wrt gamma control points
            using pre-calculated MBI jacobian
        """

        dgamma = self.vec['dg']('gamma')
        dh_pt = self.vec['dp']('h_pt')

        if self.mode == 'fwd':
            dgamma[:] = 0.0
            if self.get_id('h_pt') in args:
                dgamma[:] += self.jac_gamma.dot(dh_pt[:]) * 1e3/1e-1
        if self.mode == 'rev':
            dh_pt[:] = 0.0
            if self.get_id('h_pt') in args:
                dh_pt[:] += self.jac_gamma.T.dot(dgamma[:]) * 1e3/1e-1
