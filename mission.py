from __future__ import division
import sys
sys.path.insert(0, '/home/jason/github/CMF')
import numpy
import copy
from framework import *
from optimization import *
from bsplines import *
from atmospherics import *
from coupled_analysis import *
from functionals import *
from aerodynamics import *
from propulsion import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab
import MBI, scipy.sparse

class GlobalizedSystem(SerialSystem):
    ''' doc string '''
    
    def solve_F(self):
        """ Solve f for u, p |-> u """

        kwargs = self.kwargs
        self.solvers['NL']['NLN_GS'](ilimit=kwargs['GL_GS_ilimit'],
                                     atol=kwargs['GL_GS_atol'],
                                     rtol=kwargs['GL_GS_rtol'])
        return self.solvers['NL']['NEWTON'](ilimit=kwargs['GL_NT_ilimit'],
                                            atol=kwargs['GL_NT_atol'],
                                            rtol=kwargs['GL_NT_rtol'])
    
    
class OptTrajectory(object):
    """ class used to define and setup trajectory optimization problem """

    def __init__(self, num_elem, num_cp):
        self.num_elem = num_elem
        self.num_pt = num_cp
        self.h_pts = numpy.zeros(num_elem+1)
        self.M_pts = numpy.zeros(num_elem+1)
        self.x_pts = numpy.zeros(num_elem+1)
        self.wing_area = 0.0
        self.ac_w = 0.0
        self.thrust_sl = 0.0
        self.sfc_sl = 0.0
        self.aspect_ratio = 0.0
        self.oswald = 0.0

    def set_init_h(self, h_init):
        self.h_pts = h_init

    def set_init_M(self, M_init):
        self.M_pts = M_init

    def set_init_x(self, x_init):
        self.x_pts = x_init

    def set_params(self, kw):
        self.wing_area = kw['S']
        self.ac_w = kw['ac_w']
        self.thrust_sl = kw['thrust_sl']
        self.sfc_sl = kw['SFCSL']
        self.aspect_ratio = kw['AR']
        self.oswald = kw['e']

    def initialize(self):

        self.main = SerialSystem('mission',
                                 NL='NLN_GS',
                                 LN='LIN_GS',
                                 LN_ilimit=2,
                                 NL_ilimit=2,
                                 NL_rtol=1e-6,
                                 NL_atol=1e-10,
                                 LN_rtol=1e-6,
                                 LN_atol=1e-10,
                                 output=True,
                                 subsystems=[
                SerialSystem('mission_param',
                             NL='NLN_GS',
                             LN='LIN_GS',
                             LN_ilimit=2,
                             NL_ilimit=2,
                             output=True,
                             subsystems=[
                        IndVar('S', val=self.wing_area, size=1),
                        IndVar('ac_w', val=self.ac_w, size=1),
                        IndVar('thrust_sl', val=self.thrust_sl, size=1),
                        IndVar('SFCSL', val=self.sfc_sl, size=1),
                        IndVar('AR', val=self.aspect_ratio, size=1),
                        IndVar('e', val=self.oswald, size=1),
                        ]),
                SerialSystem('segment',
                             NL='NLN_GS',
                             LN='LIN_GS',
                             LN_ilimit=2,
                             NL_ilimit=2,
                             output=True,
                             subsystems=[
                        SerialSystem('bsplines',
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     LN_ilimit=2,
                                     NL_ilimit=2,
                                     output=True,
                                     subsystems=[
                                IndVar('x_pt', val=self.x_pts, lower=0),
                                IndVar('h_pt', val=self.h_pts, lower=0),
                                IndVar('M_pt', val=self.M_pts, lower=0),
                                SysXBspline('x', num_elem=self.num_elem,
                                            num_pt=self.num_pt,
                                            x_range=self.x_pts[-1],
                                            x_0=numpy.linspace(0.0, self.x_pts[-1], self.num_elem+1)),
                                SysHBspline('h', num_elem=self.num_elem,
                                            num_pt=self.num_pt,
                                            x_range=self.x_pts[-1]),
                                SysMBspline('M', num_elem=self.num_elem,
                                            num_pt=self.num_pt,
                                            x_range=self.x_pts[-1]),
                                SysGammaBspline('gamma', num_elem=self.num_elem,
                                                num_pt=self.num_pt,
                                                x_range=self.x_pts[-1]),
                                ]),
                        SerialSystem('atmosphere',
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     LN_ilimit=2,
                                     NL_ilimit=2,
                                     output=True,
                                     subsystems=[
                                SysSFC('SFC', num_elem=self.num_elem),
                                SysTemp('Temp', num_elem=self.num_elem),
                                SysRho('rho', num_elem=self.num_elem),
                                SysSpeed('v', num_elem=self.num_elem),
                                ]),
                        SerialSystem('coupled_analysis',
                                     NL='NEWTON',   
                                     LN='KSP_PC',
                                     PC='LIN_GS',
                                     LN_ilimit=30,
                                     NL_ilimit=30,
                                     #GL_GS_ilimit=5,
                                     #GL_NT_ilimit=30,
                                     PC_ilimit=3,
                                     NL_rtol=1e-10,
                                     NL_atol=1e-10,
                                     #GL_GS_rtol=1e-6,
                                     #GL_GS_atol=1e-10,
                                     #GL_NT_rtol=1e-10,
                                     #GL_NT_atol=1e-10,
                                     LN_rtol=1e-14,
                                     LN_atol=1e-14,
                                     PC_rtol=1e-6,
                                     PC_atol=1e-10,
                                     output=True,
                                     subsystems=[
                                SerialSystem('vert_eqlm',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=2,
                                             NL_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysCLTar('CL_tar',
                                                 num_elem=self.num_elem),
                                        ]),
                                SerialSystem('drag',
                                             NL='NEWTON',
                                             LN='KSP_PC',
                                             LN_ilimit=30,
                                             NL_ilimit=30,
                                             PC_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysAeroSurrogate('CL', num_elem=self.num_elem),
                                        SysAlpha('alpha',
                                                 num_elem=self.num_elem),
                                        ]),
                                SerialSystem('hor_eqlm',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=2,
                                             NL_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysCTTar('CT_tar',
                                                 num_elem=self.num_elem),
                                        ]),
                                SerialSystem('mmt_eqlm',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=2,
                                             NL_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysCM('eta', num_elem=self.num_elem),
                                        ]),
                                SerialSystem('weight',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=2,
                                             NL_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysFuelWeight('fuel_w',
                                                      num_elem=self.num_elem,
                                                      fuel_w_0=numpy.linspace(1.0, 0.0, self.num_elem+1)),
                                        ]),
                                ]),
                        SerialSystem('functionals',
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     NL_ilimit=2,
                                     LN_ilimit=2,
                                     NL_rtol=1e-10,
                                     NL_atol=1e-10,
                                     LN_rtol=1e-10,
                                     LN_ato1=1e-10,
                                     output=True,
                                     subsystems=[
                                SysTau('tau', num_elem=self.num_elem),
                                SysFuelObj('wf_obj'),
                                SysHi('h_i'),
                                SysHf('h_f', num_elem=self.num_elem),
                                SysTmin('Tmin', num_elem=self.num_elem),
                                SysTmax('Tmax', num_elem=self.num_elem),
                                ]),
                        ]),
                 ]).setup()

        return self.main
