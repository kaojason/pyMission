"""
INTENDED FOR MISSION ANALYSIS USE
provides propulsion models for the use of mission analysis.
The mission analysis and trajectory optimization tool was developed by:
    Jason Kao*
    John Hwang*

* University of Michigan Department of Aerospace Engineering,
  Multidisciplinary Design Optimization lab
  mdolab.engin.umich.edu

copyright July 2014
"""

# pylint: disable=E1101
from __future__ import division
import sys
from framework import *
import numpy

def setup_prop_surrogate(prop_file):

    alt_num = 11
    mach_num = 10

    tmp = numpy.loadtxt('PAX300.outputFLOPS')

    tmp[:, 3] = (tmp[:, 3] - tmp[:, 4]) # total thrust = thrust - ram drag

    for i in xrange(len(tmp)):
        if tmp[i, 2] == 50:
            tmax = tmp[i, 3]
        tmp[i, 2] = tmp[i, 3] / tmax # change from power code to throttle

    output_array = numpy.zeros((len(tmp), 5))
    output_array[:, 0] = tmp[:, 1] # altitude
    output_array[:, 1] = tmp[:, 0] # mach number
    output_array[:, 2] = tmp[:, 2] # throttle
    output_array[:, 3] = tmp[:, 3] # thrust
    output_array[:, 4] = tmp[:, 6] # TSFC

    mbi_Thrust = numpy.zeros((alt_num, mach_num, throttle_num))
    mbi_TSFC = numpy.zeros((alt_num, mach_num, throttle_num))

    #Thrust_arr = MBI.MBI(mbi_Thrust, [

class SysSFC(ExplicitSystem):
    """ linear SFC model wrt altitude """

    def _declare(self):
        """ owned variable: SFC (specific fuel consumption)
            dependencies: h (altitude)
                          SFCSL (sea-level SFC value)
        """

        self.num_elem = self.kwargs['num_elem']
        self.SFCSL = self.kwargs['SFCSL']
        num_pts = self.num_elem+1
        ind_pts = range(num_pts)

        self._declare_variable('SFC', size=num_pts)
        self._declare_argument('h', indices=ind_pts)

    def apply_G(self):
        """ compute SFC value using sea level SFC and altitude
            the model is a linear correction for altitude changes
        """

        pvec = self.vec['p']
        uvec = self.vec['u']
        alt = pvec('h') * 1e3
        sfcsl = self.SFCSL * 1e-6
        sfc = uvec('SFC')

        sfc_temp = sfcsl + (6.39e-13*9.81) * alt
        sfc[:] = sfc_temp / 1e-6

    def apply_dGdp(self, args):
        """ compute SFC derivatives wrt sea level SFC and altitude """

        dpvec = self.vec['dp']
        dgvec = self.vec['dg']

        dalt = dpvec('h')
        dsfc = dgvec('SFC')

        dsfc_dalt = 6.39e-13 * 9.81

        if self.mode == 'fwd':
            dsfc[:] = 0.0
            if self.get_id('h') in args:
                dsfc[:] += (dsfc_dalt * dalt) * 1e3/1e-6

        if self.mode == 'rev':
            dalt[:] = 0.0

            if self.get_id('h') in args:
                dalt[:] += dsfc_dalt * dsfc * 1e3/1e-6

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
        self.thrust_sl = self.kwargs['thrust_sl']
        self.wing_area = self.kwargs['S']
        num_pts = self.num_elem+1
        ind_pts = range(num_pts)

        self._declare_variable('tau', size=num_pts)
        self._declare_argument('CT_tar', indices=ind_pts)
        self._declare_argument('rho', indices=ind_pts)
        self._declare_argument('v', indices=ind_pts)
        self._declare_argument('h', indices=ind_pts)

    def apply_G(self):
        """ compute throttle setting primarily using thrust coefficient """

        pvec = self.vec['p']
        uvec = self.vec['u']

        thrust_c = pvec('CT_tar') * 1e-1
        rho = pvec('rho')
        speed = pvec('v') * 1e2
        alt = pvec('h') * 1e3
        thrust_sl = self.thrust_sl * 1e6
        wing_area = self.wing_area * 1e2
        tau = uvec('tau')

        cThrust = thrust_sl - 72 * alt
        Thrust = 0.5*rho*speed**2*wing_area*thrust_c
        tau[:] = (Thrust / cThrust)

    def linearize(self):
        """ pre-compute the throttle derivatives wrt density, velocity
            wing area, thrust coefficient, sea level thrust, and altitude
        """

        pvec = self.vec['p']

        thrust_sl = self.thrust_sl * 1e6
        wing_area = self.wing_area * 1e2
        alt = pvec('h') * 1e3
        thrust_c = pvec('CT_tar') * 1e-1
        rho = pvec('rho')
        speed = pvec('v') * 1e2

        self.dt_drho = ((0.5*speed**2*wing_area*thrust_c) / (thrust_sl-72*alt))
        self.dt_dspeed = ((rho*speed*wing_area*thrust_c) / (thrust_sl-72*alt))
        self.dt_dthrust_c = ((0.5*rho*speed**2*wing_area) / (thrust_sl-72*alt))
        self.dt_dalt = 72 * ((0.5*rho*speed**2*wing_area*thrust_c) /\
                       (thrust_sl-72*alt)**2)

    def apply_dGdp(self, arguments):
        """ assign throttle directional derivatives """

        dpvec = self.vec['dp']
        dgvec = self.vec['dg']

        dalt = dpvec('h')
        dthrust_c = dpvec('CT_tar')
        drho = dpvec('rho')
        dspeed = dpvec('v')
        dtau = dgvec('tau')

        dt_dalt = self.dt_dalt
        dt_dthrust_c = self.dt_dthrust_c
        dt_drho = self.dt_drho
        dt_dspeed = self.dt_dspeed

        if self.mode == 'fwd':
            dtau[:] = 0.0
            if self.get_id('h') in arguments:
                dtau[:] += (dt_dalt * dalt) * 1e3
            if self.get_id('CT_tar') in arguments:
                dtau[:] += (dt_dthrust_c * dthrust_c) * 1e-1
            if self.get_id('rho') in arguments:
                dtau[:] += dt_drho * drho
            if self.get_id('v') in arguments:
                dtau[:] += (dt_dspeed * dspeed) * 1e2
        if self.mode == 'rev':
            dalt[:] = 0.0
            dthrust_c[:] = 0.0
            drho[:] = 0.0
            dspeed[:] = 0.0
            if self.get_id('h') in arguments:
                dalt[:] += dt_dalt * dtau * 1e3
            if self.get_id('CT_tar') in arguments:
                dthrust_c[:] += dt_dthrust_c * dtau * 1e-1
            if self.get_id('rho') in arguments:
                drho[:] += dt_drho * dtau
            if self.get_id('v') in arguments:
                dspeed[:] += dt_dspeed * dtau * 1e2

class SysTauSurrogate(ExplicitSystem):
    """ compute the throttle setting from target CT by using existing
        engine data
    """

    def _declare(self):
        """ owned variables: tau (throttle setting)
            dependencies: h (altitude)
                          temp (temperature)
                          v (speed)
                          CT (coefficient of thrust)
                          rho (density of air)
                          S (wing area)
        """

        self.num_elem = self.kwargs['num_elem']
        num_pts = self.num_elem+1
        ind_pts = self.range(num_pts)

        self._declare_variable('tau', size=num_pts)
        self._declare_argument('h', indices=ind_pts)
        self._declare_argument('Temp', indices=ind_pts)
        self._declare_argument('v', indices=ind_pts)
        self._declare_argument('CT_tar', indices=ind_pts)
        self._declare_argument('rho', indices=ind_pts)
        self._declare_argument(['S', 0], indices=[0])

        self.build_surrogate('UHB.outputFLOPS')

    def build_surrogate(self, file_name):
        """ builds the surrogate model from the data stored in the file name
            given in the input arguments
        """

        data_file = open(file_name, 'r')

        for i, l in enumerate(data_file):
            pass

        file_len = i+1
        mach = numpy.zeros(file_len)
        altitude = numpy.zeros(file_len)
        power_code = numpy.zeros(file_len)
        thrust = numpy.zeros(file_len)
        drag = numpy.zeros(file_len)
        TSFC = numpy.zeros(file_len)
        i = 0

        data_file = open(file_name, 'r')

        for line in data_file:
            [mach[i], altitude[i], power_code[i], thrust[i], drag[i], fuel_burn,
             TSFC[i], Nox, area] = line.split()
            i += 1

        mach = [float(i) for i in mach]
        altitude = [float(i) for i in altitude]
        power_code = [float(i) for i in power_code]
        thrust = [float(i) for i in thrust]
        drag = [float(i) for i in drag]
        TSFC = [float(i) for i in TSFC]



#    def apply_G(self):
