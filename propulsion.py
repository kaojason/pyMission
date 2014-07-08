"""
comments
"""

# pylint: disable=E1101
from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
from framework import *
import numpy

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
        dt_dS = self.dt_dS

        if self.mode == 'fwd':
            dtau[:] = 0.0
            if self.get_id('S') in arguments:
                dtau[:] += (dt_dS * dwing_area) * 1e2
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
            if self.get_id('S') in arguments:
                dwing_area[:] += numpy.sum(dt_dS * dtau) * 1e2
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

