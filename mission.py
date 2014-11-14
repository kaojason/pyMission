from __future__ import division
import numpy
import copy
import os
from framework import *
from optimization import *
from bsplines import *
from atmospherics import *
from coupled_analysis import *
from functionals import *
from aerodynamics_old import *
from propulsion import *
from aeroTripan import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab
import MBI, scipy.sparse
import scipy.sparse.linalg
from history import *

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

    def __init__(self, num_elem, num_cp, first):
        self.num_elem = num_elem
        self.num_cp = num_cp
        self.h_pts = numpy.zeros(num_cp)
        self.M_pts = numpy.zeros(num_cp)
        self.v_pts = numpy.zeros(num_cp)
        self.x_pts = numpy.zeros(num_cp)
        self.wing_area = 0.0
        self.ac_w = 0.0
        self.thrust_sl = 0.0
        self.sfc_sl = 0.0
        self.aspect_ratio = 0.0
        self.oswald = 0.0
        self.v_specified = 0
        self.gamma_lb = numpy.tan(-30.0 * (numpy.pi/180.0))/1e-1
        self.gamma_ub = numpy.tan(30.0 * (numpy.pi/180.0))/1e-1
        self.folder_path = None
        self.name = None
        self.first = first

    def set_init_h(self, h_init):
        self.h_pts = h_init

    def set_init_h_pt(self, h_init_pt):
        A = self.jac_h
        b = h_init_pt
        ATA = A.T.dot(A)
        ATb = A.T.dot(b)
        self.h_pts = scipy.sparse.linalg.gmres(ATA, ATb)[0]

    def set_init_M(self, M_init):
        self.M_pts = M_init

    def set_init_x(self, x_init):
        self.x_pts = x_init

    def set_init_v(self, v_init):
        self.v_pts = v_init
        self.v_specified = 1

    def setup_MBI(self):
        """ generate jacobians for b-splines using MBI package """

        num_pts = self.num_elem + 1
        num_cp = self.num_cp

        alt = numpy.linspace(0, 16, num_pts)
        x_dist = numpy.linspace(0, self.x_pts[-1], num_pts)/1e6

        arr = MBI.MBI(alt, [x_dist], [num_cp], [4])
        jac = arr.getJacobian(0, 0)
        jacd = arr.getJacobian(1, 0)

        c_arryx = self.x_pts
        d_arryx = jacd.dot(c_arryx)*1e6

        lins = numpy.linspace(0, num_pts-1, num_pts).astype(int)
        diag = scipy.sparse.csc_matrix((1.0/d_arryx,
                                        (lins,lins)))
        jace = diag.dot(jacd)

        self.jac_h = jac
        self.jac_gamma = jace

        [self.CL_arr, self.CD_arr, self.CM_arr, self.num] = \
            setup_surrogate(self.surr_file)
        
    def set_params(self, kw):
        self.wing_area = kw['S']
        self.ac_w = kw['ac_w']
        self.thrust_sl = kw['thrust_sl']
        self.sfc_sl = kw['SFCSL']
        self.aspect_ratio = kw['AR']
        self.oswald = kw['e']
        self.surr_file = kw['surr']

    def set_folder(self, folder_path):
        self.folder_path = folder_path

    def set_name(self, name):
        self.name = name

    def initialize_framework(self):

        self.main = SerialSystem('mission',
                                 NL='NLN_GS',
                                 LN='LIN_GS',
                                 LN_ilimit=1,
                                 NL_ilimit=1,
                                 NL_rtol=1e-6,
                                 NL_atol=1e-9,
                                 LN_rtol=1e-6,
                                 LN_atol=1e-10,
                                 output=True,
                                 subsystems=[
                SerialSystem('mission_param',
                             NL='NLN_GS',
                             LN='LIN_GS',
                             LN_ilimit=1,
                             NL_ilimit=1,
                             output=True,
                             subsystems=[
                        IndVar('jason', val=100000, size=1),
                        #IndVar('lambda', 
                        #       val=numpy.ones(2*(self.num_elem+1)), 
                        #       lower=0),
                        #IndVar('slack',
                        #       val=numpy.zeros(2*(self.num_elem+1)),
                        #       lower=0),
                        ]),
                SerialSystem('segment',
                             NL='NLN_GS',
                             LN='LIN_GS',
                             LN_ilimit=1,
                             NL_ilimit=1,
                             NL_atol=1e-9,
                             NL_rtol=1e-9,
                             output=True,
                             subsystems=[
                        SerialSystem('bsplines',
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     LN_ilimit=1,
                                     NL_ilimit=1,
                                     output=True,
                                     subsystems=[
                                IndVar('x_pt', val=self.x_pts, lower=0),
                                IndVar('h_pt', val=self.h_pts, lower=0),
                                IndVar('M_pt', val=self.M_pts, lower=0),
                                IndVar('v_pt', val=self.v_pts, lower=0),
                                SysXBspline('x', num_elem=self.num_elem,
                                            num_cp=self.num_cp,
                                            x_init=self.x_pts,
                                            x_0=numpy.linspace(0.0, self.x_pts[-1], self.num_elem+1),
                                            jac_h=self.jac_h),
                                SysHBspline('h', num_elem=self.num_elem,
                                            num_cp=self.num_cp,
                                            x_init=self.x_pts,
                                            jac_h=self.jac_h),
                                SysMVBspline('MV', num_elem=self.num_elem,
                                             num_cp=self.num_cp,
                                             x_init=self.x_pts,
                                             jac_h=self.jac_h),
                                SysGammaBspline('gamma',
                                                num_elem=self.num_elem,
                                                num_cp=self.num_cp,
                                                x_init=self.x_pts,
                                                jac_gamma=self.jac_gamma),
                                ]),
                        SerialSystem('atmospherics',
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     LN_ilimit=1,
                                     NL_ilimit=1,
                                     output=True,
                                     subsystems=[
                                SysSFC('SFC', num_elem=self.num_elem,
                                       SFCSL=self.sfc_sl),
                                SysTemp('Temp', num_elem=self.num_elem),
                                SysRho('rho', num_elem=self.num_elem),
                                SysSpeed('v', num_elem=self.num_elem,
                                         v_specified=self.v_specified),
                                SysMach('M', num_elem=self.num_elem,
                                        v_specified=self.v_specified),
                                ]),
                        SerialSystem('coupled_analysis',
                                     LN='KSP_PC',
                                     #PC='LIN_GS',
                                     LN_ilimit=50,
                                     GL_GS_ilimit=1,
                                     GL_NT_ilimit=8,
                                     PC_ilimit=2,
                                     GL_GS_rtol=1e-6,
                                     GL_GS_atol=1e-9,#10,
                                     GL_NT_rtol=1e-9,#14,
                                     GL_NT_atol=1e-9,#14,
                                     NL_rtol=1e-9,
                                     NL_atol=1e-9,
                                     LN_rtol=1e-20,#14,
                                     LN_atol=1e-14,#14,
                                     PC_rtol=1e-6,
                                     PC_atol=1e-10,
                                     output=True,
                                     subsystems=[
                                SerialSystem('vert_eqlm',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=1,
                                             NL_ilimit=1,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysCLTar('CL_tar',
                                                 num_elem=self.num_elem,
                                                 S=self.wing_area,
                                                 ac_w = self.ac_w),
                                        ]),
                                SerialSystem('tripan_alpha',
                                             NL='NLN_GS',
                                             LN='LIN_GS',
                                             LN_ilimit=18,
                                             NL_ilimit=1,
                                             PC_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-6,
                                             LN_atol=1e-6,
                                             subsystems=[
                                        SysTripanCLSurrogate('alpha', num_elem=self.num_elem, num=self.num, CL=self.CL_arr),
                                        #SysAeroSurrogate('aero', num_elem=self.num_elem),
                                        #SysAlpha('alpha', num_elem=self.num_elem),
                                        #SysCM('eta', num_elem=self.num_elem),
                                        ]),
                                SerialSystem('tripan_eta',
                                             NL='NLN_GS',
                                             LN='LIN_GS',
                                             LN_ilimit=18,
                                             NL_ilimit=1,
                                             PC_ilimit=2,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-6,
                                             LN_atol=1e-6,
                                             subsystems=[
                                        SysTripanCMSurrogate('eta', num_elem=self.num_elem, num=self.num, CM=self.CM_arr),
                                        ]),
                                SerialSystem('tripan_drag',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=1,
                                             NL_ilimit=1,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysTripanCDSurrogate('drag', num_elem=self.num_elem, num=self.num, CD=self.CD_arr),
                                        ]),
                                SerialSystem('hor_eqlm',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=1,
                                             NL_ilimit=1,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysCTTar('CT_tar',
                                                 num_elem=self.num_elem,
                                                 S=self.wing_area,
                                                 ac_w=self.ac_w),
                                        ]),
                                SerialSystem('weight',
                                             NL='NLN_GS',
                                             LN='KSP_PC',
                                             LN_ilimit=1,
                                             NL_ilimit=1,
                                             NL_rtol=1e-10,
                                             NL_atol=1e-10,
                                             LN_rtol=1e-10,
                                             LN_atol=1e-10,
                                             subsystems=[
                                        SysFuelWeight('fuel_w',
                                                      num_elem=self.num_elem,
                                                      fuel_w_0=numpy.linspace(1.0, 0.0, self.num_elem+1),
                                                      S=self.wing_area),
                                        ]),
                                ]),
                        SerialSystem('functionals',
                                     NL='NLN_GS',
                                     LN='LIN_GS',
                                     NL_ilimit=1,
                                     LN_ilimit=1,
                                     NL_rtol=1e-10,
                                     NL_atol=1e-10,
                                     LN_rtol=1e-10,
                                     LN_ato1=1e-10,
                                     output=True,
                                     subsystems=[
                                SysTau('tau', num_elem=self.num_elem,
                                       thrust_sl=self.thrust_sl,
                                       S=self.wing_area),
                                SysFuelObj('fuelburn', num_elem=self.num_elem),
                                SysHi('h_i'),
                                SysHf('h_f', num_elem=self.num_elem),
                                SysTmin('Tmin', num_elem=self.num_elem),
                                SysTmax('Tmax', num_elem=self.num_elem),
                                SysMi('M_i'),
                                SysMf('M_f', num_elem=self.num_elem),
                                SysBlockTime('time', num_elem=self.num_elem),
                                ]),
                        ]),
                 ]).setup()

        self.history = History(self.num_elem, self.num_cp, self.x_pts[-1], 
                               self.folder_path, self.name, self.first)

        return self.main

    def set_gamma_bound(self, gamma_lb, gamma_ub):
        self.gamma_lb = gamma_lb
        self.gamma_ub = gamma_ub

    def set_takeoff_speed(self, v_to):
        self.v_to = v_to
        self.M_to = v_to / numpy.sqrt(1.4*287*288)

    def set_landing_speed(self, v_ld):
        self.v_ld = v_ld
        self.M_ld = v_ld / numpy.sqrt(1.4*287*288)

    def initialize_opt(self, main):
        gamma_lb = self.gamma_lb
        gamma_ub = self.gamma_ub

        opt = Optimization(main)
        opt.add_design_variable('h_pt', value=self.h_pts, lower=0.0, upper=14.1)
        #opt.add_design_variable('M_pt', value=self.M_pts, lower=0.001, upper=0.949)
        #opt.add_design_variable('lambda',
        #                        value=numpy.zeros(2*(self.num_elem+1)),
        #                        lower=0.0)
        #opt.add_design_variable('slack',
        #                        value=numpy.zeros(2*(self.num_elem+1)),
        #                        lower=0.0)
        #opt.add_design_variable('v_pt', value=self.v_pts, lower=0.0, upper=10.0)
        opt.add_objective('fuelburn')
        opt.add_constraint('h_i', lower=0.0, upper=0.0)
        opt.add_constraint('h_f', lower=0.0, upper=0.0)
        opt.add_constraint('Tmin', upper=0.0)
        opt.add_constraint('Tmax', upper=0.0)
        opt.add_constraint('gamma', lower=gamma_lb, upper=gamma_ub,
                           get_jacs=main('gamma').get_jacs, linear=True)
        #opt.add_constraint('M_i', lower=self.M_to, upper=self.M_to)
        #opt.add_constraint('M_f', lower=self.M_ld, upper=self.M_ld)
        #opt.add_constraint('time', lower=0.0, upper=11*3600.0)
        opt.add_sens_callback(self.callback)
        return opt

    def callback(self):        
        self.history.save_history(self.main.vec['u'], self.wing_area,
                                  self.ac_w)

