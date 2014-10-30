"""
INTENDED FOR MISSION ANALYSIS USE
This module provides functions and classes used to write optimization
history into files, and plot the information from these files.
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
import os
import numpy
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab

class History(object):
    """ class used to write optimization history onto disk """

    def __init__(self, num_elem, num_cp, x_range, folder_path, name, first=False):
        """ initialize variables, set folder name to 
            distxxxxkm-yyyy-zzzz-nnn/
            where xxxx is the distance of the mission
                  yyyy is the number of control points
                  zzzz is the number of elements
              and nnn is the case index to distinguish cases of the same
              parameters
        """

        self.num_elem = num_elem
        self.num_cp = num_cp
        self.x_range = x_range
        self.folder_name = folder_path + name
        self.name = name

        index = 0
        while os.path.exists(self.folder_name+'_%03i' % (index)):
            index += 1

        if first == True:
            self.folder_name = self.folder_name+'_%03i/' % (index)
            os.makedirs(self.folder_name)
            self.index = index
        else:
            index -= 1
            self.folder_name = self.folder_name+'_%03i/' % (index)
            self.index = index

        self.hist_counter = 0

        self.variable_max = numpy.zeros(16)
        self.variable_min = numpy.zeros(16)

    def get_index(self):
        """ returns the case number (self.index), and the iteration number
            (self.hist_counter)
        """

        return self.index, (self.hist_counter - 1)

    def save_history(self, vecu, S, ac_w):
        """ saves all relevant variables within the u vector into the folder
            name mentioned before, and the file name is:
            xxxxkm-yyyy-zzzz-nnnn.dat
            where xxxx is the distance of the mission
                  yyyy is the number of control points
                  zzzz is the number of elements
                  nnnn is the iteration number
        """

        dist = vecu('x') * 1e6
        altitude = vecu('h') * 1e3
        speed = vecu('v') * 1e2
        alpha = vecu('alpha') * 1e-1 * 180/numpy.pi
        throttle = vecu('tau')
        eta = vecu('eta') * 1e-1 * 180/numpy.pi
        fuel = vecu('fuel_w') * 1e5
        rho = vecu('rho')
        thrust = vecu('CT_tar')*0.5*rho*speed**2*S*1e2 * 1e-1
        drag_c = vecu('CD') * 1e-1
        lift_c = vecu('CL_tar')
        gamma = vecu('gamma') * 1e-1 * 180/numpy.pi
        weight = (ac_w*1e6 + vecu('fuel_w')*1e5)
        temp = vecu('Temp') * 1e2
        SFC = vecu('SFC') * 1e-6

        file_name = self.name + '_%04i_%04i.dat' % (self.num_cp,
                                                  self.hist_counter)

        output_file = self.folder_name + file_name

        file_array = [dist, altitude, speed, alpha, throttle, eta, fuel,
                      rho, lift_c, drag_c, thrust, gamma, weight,
                      temp, SFC]
        numpy.savetxt(output_file, file_array)

        self.hist_counter += 1

    def print_max_min(self, vecu, kw):
        """ print the maximum and the minimum of each variable throughout
            the optimization history into the file
            xxxxkm-yyyy-zzzz-maxmin.dat
            where xxxx is the distance of the mission
                  yyyy is the number of control points
                  zzzz is the number of elements
        """

        S = kw['S']
        ac_w = kw['ac_w']

        dist = vecu('x') * 1e6
        altitude = vecu('h') * 1e3
        speed = vecu('v') * 1e2
        alpha = vecu('alpha') * 1e-1 * 180/numpy.pi
        throttle = vecu('tau')
        eta = vecu('eta') * 1e-1 * 180/numpy.pi
        fuel = vecu('fuel_w') * 1e5
        rho = vecu('rho')
        thrust = vecu('CT_tar')*0.5*rho*speed**2*S * 1e-1
        drag_c = vecu('CD') * 1e-1
        lift_c = vecu('CL_tar')
        gamma = vecu('gamma') * 1e-1 * 180/numpy.pi
        weight = (ac_w*1e6 + vecu('fuel_w')*1e5)
        temp = vecu('Temp') * 1e2
        SFC = vecu('SFC') * 1e-6

        array = [dist, altitude, speed, alpha, throttle, eta, fuel,
                 rho, lift_c, drag_c, thrust, gamma, weight,
                 temp, SFC]
        index = 0
        for variable in array:
            self.variable_max[index] = max(variable)
            self.variable_min[index] = min(variable)
            index += 1

        opt_alt = self.compute_est(vecu, kw)
        self.variable_min[index] = opt_alt[0]
        self.variable_max[index] = opt_alt[1]

        file_array = [self.variable_min, self.variable_max]
        file_name = self.name + '_maxmin.dat'
        output_file = self.folder_name+file_name
        numpy.savetxt(output_file, file_array)

    def compute_est(self, vecu, kw):
        """ computes the estimated optimal flight altitude with
            empty weight and takeoff gross weight. this data is
            intended to be stored in the max-min file and be used
            in pltscript.py to provide a sanity check on the
            results on the final plot.
        """

        self.epsilon = 500
        h_lower = 11000 - self.epsilon
        h_upper = 11000 + self.epsilon
        matrix = numpy.array([[h_lower**3, h_lower**2, h_lower, 1],
                              [h_upper**3, h_upper**2, h_upper, 1],
                              [3*h_lower**2, 2*h_lower, 1, 0],
                              [3*h_upper**2, 2*h_upper, 1, 0]])
        rhs = numpy.array([288.16-(6.5e-3)*h_lower, 216.65,
                           -6.5e-3, 0])
        coefs_temp = numpy.linalg.solve(matrix, rhs)

        rhs = numpy.array([101325*(1-0.0065*h_lower/288.16)**5.2561,
                           22632*numpy.exp(-9.81*self.epsilon/(288*216.65)),
                           (-101325*5.2561*(0.0065/288.16)*
                             (1-0.0065*h_lower/288.15)**4.2561),
                           (22632*(-9.81/(288*216.65))*
                            numpy.exp(-9.81*self.epsilon/(288*216.65)))])
        coefs_dens = numpy.linalg.solve(matrix, rhs)

        params = {
            'SFCSL': kw['SFCSL'] * 1e-6,
            'wing_area': kw['S'] * 1e2,
            'weight': (kw['ac_w']*1e6 + vecu('fuel_w')[0]*1e5),
            'aspect_ratio': kw['AR'],
            'oswald': kw['e'],
            'coefs_t': coefs_temp,
            'coefs_d': coefs_dens,
            }

        h_1 = 10000
        h_2 = 15000
        res = 1.0
        M_guess = 0.82
        f_1 = self.alt_est_obj(params, h_1, M_guess)
        f_2 = self.alt_est_obj(params, h_2, M_guess)

        while numpy.abs(res) > 1e-6:
            h_guess = (h_1*f_2 - h_2*f_1)/(f_2 - f_1)
            h_1 = copy.copy(h_2)
            h_2 = copy.copy(h_guess)
            f_1 = copy.copy(f_2)
            f_2 = self.alt_est_obj(params, h_2, M_guess)
            res = f_2
        h_to = copy.copy(h_guess)

        params['weight'] = kw['ac_w'] * 1e6
        h_1 = 10000
        h_2 = 15000
        res = 1.0
        f_1 = self.alt_est_obj(params, h_1, M_guess)
        f_2 = self.alt_est_obj(params, h_2, M_guess)

        while numpy.abs(res) > 1e-6:
            h_guess = (h_1*f_2 - h_2*f_1)/(f_2 - f_1)
            h_1 = copy.copy(h_2)
            h_2 = copy.copy(h_guess)
            f_1 = copy.copy(f_2)
            f_2 = self.alt_est_obj(params, h_2, M_guess)
            res = f_2
        h_empty = copy.copy(h_guess)

        return [h_empty, h_to]

    def alt_est_obj(self, params, h_guess, M_guess):
        """ the expression to be minimized to obtain the
            back-of-the-envelope estimation of the optimal
            cruise altitude
        """

        SFCSL = params['SFCSL']
        wing_area = params['wing_area']
        weight = params['weight']
        aspect_ratio = params['aspect_ratio']
        oswald = params['oswald']
        SFC_slp = 6.39e-13 * 9.81
        cd0 = 0.018
        coefs_temp = params['coefs_t']
        coefs_dens = params['coefs_d']

        if h_guess <= (11000-self.epsilon):
            temp = 288.16 - (6.5e-3) * h_guess
            dtdh = -6.5e-3
            pressure = 101325*(1-0.0065*h_guess/288.16)**5.2561
            dpdh = 101325*5.2561*(-0.0065/288.16)*\
                (1-0.0065*h_guess/288.16)**4.2561
            rho = pressure / (288 * temp)
            drho = dpdh/(288*temp) - (pressure*dtdh)/(288*temp**2)
        elif h_guess >= (11000+self.epsilon):
            temp = 216.65
            dtdh = 0.0
            pressure = 22632*numpy.exp(-9.81*(h_guess-11000)/
                                        (288*216.65))
            dpdh = (22632*(-9.81/(288*216.65))*
                    numpy.exp(9.81*11000/(288*216.65))*
                    numpy.exp(-9.81*h_guess/(288*216.65)))
            rho = pressure / (288 * temp)
            drho = dpdh/(288*temp) - (pressure*dtdh)/(288*temp**2)
        else:
            a = coefs_temp[0]
            b = coefs_temp[1]
            c = coefs_temp[2]
            d = coefs_temp[3]
            temp = a*h_guess**3 + b*h_guess**2 + c*h_guess + d
            dtdh = 3*a*h_guess**2 + 2*b*h_guess + c
            a = coefs_dens[0]
            b = coefs_dens[1]
            c = coefs_dens[2]
            d = coefs_dens[3]
            pressure = a*h_guess**3 + b*h_guess**2 + c*h_guess + d
            dpdh = 3*a*h_guess**2 + 2*b*h_guess + c
            rho = pressure / (288 * temp)
            drho = dpdh/(288*temp) - (pressure*dtdh)/(288*temp**2)

        speed = numpy.sqrt(1.4*288*temp) * M_guess
        dspeed = 0.5*M_guess*numpy.sqrt(1.4*288/temp)*dtdh

        obj = (SFC_slp * ((0.5*rho*speed**2*wing_area*cd0)/weight +
                          weight/(0.5*rho*speed**2*numpy.pi*aspect_ratio*oswald)) +
               (SFCSL + SFC_slp*h_guess)*
               ((0.5*wing_area*cd0*(drho*speed**2 + rho*speed*dspeed))/weight -
                (weight/(0.5*wing_area*numpy.pi*aspect_ratio*oswald*(rho*speed**2)**2)) *
                (drho*speed**2 + rho*speed*dspeed)))

        return obj
