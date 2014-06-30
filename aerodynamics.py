"""
INTENDED FOR MISSION ANALYSIS USE
This file contains the aerodynamic models used by the mission analysis
code. The present implementation uses linear aerodynamics.
"""

# pylint: disable=E1101
from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
from framework import *
import numpy

class CLLin(object):
    def get_CLLin(alpha, eta):
        lift_c0 = 0.26
        lift_ca = 4.24
        lift_ce = 0.27

        lift_c = lift_c0 + lift_ca * alpha + lift_ce*eta
        return lift_c

    def get_CLLin_d(alpha, eta):
        lift_ca = 4.24
        lift_ce = 0.27

        return lift_ca, lift_ce

class SysCLLin(ExplicitSystem):
    """ linear aerodynamic model for CL """
    
    def _declare(self):
        """ owned variable: CL (lift coefficient)
            dependencies: alpha (angle of attack)
                          eta (tail rotation angle)
        """

        self.num_elem = self.kwargs['num_elem']
        num_pts = self.num_elem+1
        ind_pts = range(num_pts)

        self._declare_variable('CL', size=num_pts)
        self._declare_argument('alpha', indices=ind_pts)
        self._declare_argument('eta', indices=ind_pts)

    def apply_G(self):
        """ compute CL value using linear aerodynamic model and
            values for alpha, and eta
        """

        pvec = self.vec['p']
        uvec = self.vec['u']
        
        alpha = pvec('alpha')
        eta = pvec('eta')
        lift_c = uvec('CL')

        lift_c0 = 0.26
        lift_ca = 4.24
        lift_ce = 0.27

        lift_c[:] = lift_c0 + lift_ca * alpha + lift_ce * eta

    def apply_dGdp(self, args):
        """ compute derivatives of linear CL model wrt alpha and eta """
        
        dpvec = self.vec['dp']
        dgvec = self.vec['dg']

        dalpha = dpvec('alpha')
        deta = dpvec('eta')
        dlift_c = dgvec('CL')

        lift_c0 = 0.26
        lift_ca = 4.24
        lift_ce = 0.27

        if self.mode == 'fwd':
            dlift_c[:] = 0.0
            if self.get_id('alpha') in args:
                dlift_c[:] += lift_ca * dalpha
            if self.get_id('eta') in args:
                dlift_c[:] += lift_ce * deta

        elif self.mode == 'rev':
            dalpha[:] = 0.0
            deta[:] = 0.0
            if self.get_id('alpha') in args:
                dalpha[:] += lift_ca * dlift_c
            if self.get_id('eta') in args:
                deta[:] += lift_ce * dlift_c

class SysCM(ImplicitSystem):
    """ compute the tail rotation angle necessary to maintain pitch moment
        equilibrium
    """

    def _declare(self):
        """ owned variable: eta (tail rotation angle)
            dependencies: alpha (angle of attack)
        """

        self.num_elem = self.kwargs['num_elem']
        ind_pts = range(self.num_elem+1)

        self._declare_variable('eta', size=self.num_elem+1)
        self._declare_argument('alpha', indices=ind_pts)

    def apply_F(self):
        """ compute CM value using alpha and eta, and use the CM value as
            residual for eta
        """

        pvec = self.vec['p']
        uvec = self.vec['u']
        fvec = self.vec['f']

        alpha = pvec('alpha') * 1e-1
        eta = uvec('eta') * 1e-1
        eta_res = fvec('eta')

        mmt_ca = 0.63
        mmt_ce = 1.06

        eta_res[:] = (mmt_ca * alpha + mmt_ce * eta) / 1e-1

    def apply_dFdpu(self, args):
        """ compute the derivatives of tail rotation angle wrt angle of attack
        """

        dpvec = self.vec['dp']
        duvec = self.vec['du']
        dfvec = self.vec['df']

        dalpha = dpvec('alpha')
        deta = duvec('eta')
        deta_res = dfvec('eta')

        mmt_ca = 0.63
        mmt_ce = 1.06

        if self.mode == 'fwd':
            deta_res[:] = 0.0
            if self.get_id('alpha') in args:
                deta_res[:] += mmt_ca * dalpha
            if self.get_id('eta') in args:
                deta_res[:] += mmt_ce * deta

        elif self.mode == 'rev':
            dalpha[:] = 0.0
            deta[:] = 0.0
            if self.get_id('alpha') in args:
                dalpha[:] += mmt_ca * deta_res
            if self.get_id('eta') in args:
                deta[:] += mmt_ce * deta_res

class SysCD(ExplicitSystem):
    """ simple drag model with CD0 and CDi ~ CL^2 """

    def _declare(self):
        """ declare variable and arguments for drag coefficient """

        self.num_elem = self.kwargs['num_elem']
        num_pts = self.num_elem+1
        ind_pts = range(num_pts)

        self._declare_variable('CD', size=num_pts)
        self._declare_argument('CL', indices=ind_pts)
        self._declare_argument(['AR', 0], indices=[0])
        self._declare_argument(['e', 0], indices=[0])

    def apply_G(self):
        """ compute CD using a simple drag polar and CL values """

        pvec = self.vec['p']
        uvec = self.vec['u']

        lift_c = pvec('CL')
        aspect_ratio = pvec(['AR', 0])
        oswald = pvec(['e', 0])
        drag_c = uvec('CD')

        drag_c[:] = (0.018 + lift_c**2 / (numpy.pi*aspect_ratio*oswald)) / 1e-1

    def apply_dGdp(self, arguments):
        """ compute CD derivatives wrt CL, aspect ratio, and oswald's
            efficiency
        """

        pvec = self.vec['p']
        dpvec = self.vec['dp']
        dgvec = self.vec['dg']

        aspect_ratio = pvec(['AR', 0])
        oswald = pvec(['e', 0])
        lift_c = pvec('CL')

        daspect_ratio = dpvec(['AR', 0])
        doswald = dpvec(['e', 0])
        dlift_c = dpvec('CL')
        ddrag_c = dgvec('CD')

        ddrag_c_dar = -lift_c**2 / (numpy.pi*oswald*aspect_ratio**2)
        ddrag_c_doswald = -lift_c**2 / (numpy.pi*oswald**2*aspect_ratio)
        ddrag_c_dlift_c = 2*lift_c / (numpy.pi*oswald*aspect_ratio)

        if self.mode == 'fwd':
            ddrag_c[:] = 0.0
            if self.get_id('AR') in arguments:
                ddrag_c[:] += ddrag_c_dar * daspect_ratio / 1e-1
            if self.get_id('e') in arguments:
                ddrag_c[:] += ddrag_c_doswald * doswald / 1e-1
            if self.get_id('CL') in arguments:
                ddrag_c[:] += ddrag_c_dlift_c * dlift_c / 1e-1

        if self.mode == 'rev':
            daspect_ratio[:] = 0.0
            doswald[:] = 0.0
            dlift_c[:] = 0.0

            if self.get_id('AR') in arguments:
                daspect_ratio[:] += numpy.sum(ddrag_c_dar * ddrag_c) /1e-1
            if self.get_id('e') in arguments:
                doswald[:] += numpy.sum(ddrag_c_doswald * ddrag_c) /1e-1
            if self.get_id('CL') in arguments:
                dlift_c[:] += ddrag_c_dlift_c * ddrag_c /1e-1

class SysTau(ExplicitSystem):
    """ throttle setting determined primarily by thrust coefficient
        A simple linear relationship using the sea-level max thrust
        and a linear dependence on altitude is used
    """

    def _declare(self):
        """ owned variable: tau (throttle setting)
            dependencies: CT (coefficient of thrust)
                          rho (density)
                          v (speed)
                          h (altitude)
                          thrust_sl (maximum sea-level thrust)
                          S (wing area)
        """

        self.num_elem = self.kwargs['num_elem']
        num_pts = self.num_elem+1
        ind_pts = range(num_pts)

        self._declare_variable('tau', size=num_pts)
        self._declare_argument('CT_tar', indices=ind_pts)
        self._declare_argument('rho', indices=ind_pts)
        self._declare_argument('v', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)
        self._declare_argument(['thrust_sl', 0], indices=[0])
        self._declare_argument(['S', 0], indices=[0])

    def apply_G(self):
        """ compute throttle setting primarily using thrust coefficient """

        pvec = self.vec['p']
        uvec = self.vec['u']

        thrust_c = pvec('CT_tar') * 1e-1
        rho = pvec('rho')
        speed = pvec('v') * 1e2
        alt = pvec('h') * 1e3
        thrust_sl = pvec(['thrust_sl', 0]) * 1e6
        wing_area = pvec(['S', 0]) * 1e2
        tau = uvec('tau')

        cThrust = thrust_sl - 0.072 * alt
        Thrust = 0.5*rho*speed**2*wing_area*thrust_c
        tau[:] = Thrust / cThrust

    def linearize(self):
        """ pre-compute the throttle derivatives wrt density, velocity
            wing area, thrust coefficient, sea level thrust, and altitude
        """

        pvec = self.vec['p']

        thrust_sl = pvec(['thrust_sl', 0]) * 1e6
        wing_area = pvec(['S', 0]) * 1e2
        alt = pvec('h') * 1e3
        thrust_c = pvec('CT_tar') * 1e-1
        rho = pvec('rho')
        speed = pvec('v') * 1e2

        self.dt_drho = (0.5*speed**2*wing_area*thrust_c) / (thrust_sl-0.072*alt)
        self.dt_dspeed = (rho*speed*wing_area*thrust_c) / (thrust_sl-0.072*alt)
        self.dt_dS = (0.5*rho*speed**2*thrust_c) / (thrust_sl-0.072*alt)
        self.dt_dthrust_c = (0.5*rho*speed**2*wing_area) / (thrust_sl-0.072*alt)
        self.dt_dthrust_sl = -(0.5*rho*speed**2*wing_area*thrust_c) /\
                             (thrust_sl-0.072*alt)**2
        self.dt_dalt = 0.072 * (0.5*rho*speed**2*wing_area*thrust_c) /\
                       (thrust_sl-0.072*alt)**2

    def apply_dGdp(self, arguments):
        """ assign throttle directional derivatives """

        dpvec = self.vec['dp']
        dgvec = self.vec['dg']

        dthrust_sl = dpvec(['thrust_sl', 0])
        dwing_area = dpvec(['S', 0])
        dalt = dpvec('h')
        dthrust_c = dpvec('CT_tar')
        drho = dpvec('rho')
        dspeed = dpvec('v')
        dtau = dgvec('tau')

        dt_dthrust_sl = self.dt_dthrust_sl
        dt_dalt = self.dt_dalt
        dt_dthrust_c = self.dt_dthrust_c
        dt_drho = self.dt_drho
        dt_dspeed = self.dt_dspeed

        if self.mode == 'fwd':
            dtau[:] = 0.0
            if self.get_id('thrust_sl') in arguments:
                dtau[:] += (dt_dthrust_sl *
                            dthrust_sl) * 1e6
            if self.get_id('h') in arguments:
                dtau[:] += (dt_dalt * dalt) * 1e3
            if self.get_id('CT_tar') in arguments:
                dtau[:] += (dt_dthrust_c * dthrust_c) * 1e-1
            if self.get_id('rho') in arguments:
                dtau[:] += dt_drho * drho
            if self.get_id('v') in arguments:
                dtau[:] += (dt_dspeed * dspeed) * 1e2
        if self.mode == 'rev':
            dthrust_sl[:] = 0.0
            dalt[:] = 0.0
            dthrust_c[:] = 0.0
            drho[:] = 0.0
            dspeed[:] = 0.0
            dwing_area[:] = 0.0

            if self.get_id('thrust_sl') in arguments:
                dthrust_sl[:] += numpy.sum(dt_dthrust_sl * dtau) * 1e6
            if self.get_id('h') in arguments:
                dalt[:] += dt_dalt * dtau * 1e3
            if self.get_id('CT_tar') in arguments:
                dthrust_c[:] += dt_dthrust_c * dtau * 1e-1
            if self.get_id('rho') in arguments:
                drho[:] += dt_drho * dtau
            if self.get_id('v') in arguments:
                dspeed[:] += dt_dspeed * dtau * 1e2

