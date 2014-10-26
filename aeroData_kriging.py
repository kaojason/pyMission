# =============================================================================                                                                                                        
# Standard Python modules                                                                                                         
# ============================================================================= 

import os, sys, string, pdb, copy, time, datetime,shutil

# =============================================================================
# External Python modules                                                                                                                  
# =============================================================================

from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import random
import numpy
import matplotlib
import pickle
import argparse

from pykriging import sampling, Kriging

def setupGlobalModel(Data, theta0, GE=False, bf_flag=False, corrfunction=0, scaling=True):
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

    model = Kriging(theta0, GE=GE, bf_flag=bf_flag, corrfunction=corrfunction)
    
    if scaling:
        model.setSampleData(x, y, Lower=numpy.array(lb), Upper=numpy.array(ub))
    else:
        model.setSampleData(x, y)
    model.build()

    return model

# Load data to build kriging model for CL 

filename = 'BWBTransport_4D_CLkriging_coarse.pkl'
picklefile = open(filename, 'r')

CLData = pickle.load(picklefile)
picklefile.close()

# Load data to build kriging model for CM

filename = 'BWBTransport_4D_CMkriging_coarse.pkl'
picklefile = open(filename, 'r')

CMData = pickle.load(picklefile)
picklefile.close()

# Load data to build kriging model for CD

filename = 'BWBTransport_4D_CDGEK_coarse.pkl'
picklefile = open(filename, 'r')

CDData = pickle.load(picklefile)

picklefile.close()

Ndv = 4
theta0 = 10.0*numpy.ones(Ndv)

# build ordinary kriging model for CL 
CLmodel = setupGlobalModel(CLData, theta0 = theta0, GE=False, bf_flag=False)

# build ordinary kriging model for CM
CMmodel = setupGlobalModel(CMData, theta0 = theta0, GE=False, bf_flag=False)

# build GEK model for CD
CDmodel = setupGlobalModel(CDData, theta0 = theta0, GE=True, bf_flag=False)

xevals = numpy.array([[0.75, 0.0, 35000.0, 5.0], [0.82, 4.0, 40000.0, 10.0]])
# xevals = evaluation point locations, size [Nevals, Ndv]
# [Mach, angle of attack, altitude, tail angle]

CLmodel.evaluate(xevals)
CMmodel.evaluate(xevals)
CDmodel.evaluate(xevals)

CLmodel.computeGradient(xevals)
CMmodel.computeGradient(xevals)
CDmodel.computeGradient(xevals)

# size of function approximations [Nevals]
CLfunc = CLmodel.yhat.copy()
CMfunc = CMmodel.yhat.copy()
CDfunc = CDmodel.yhat.copy()

# size of gradient approximation [Nevals, Ndv]
# gradients w.r.t [Mach, angle of attack, altitude, tail angle]
CLgrad = CLmodel.gradients.copy()
CMgrad = CMmodel.gradients.copy()
CDgrad = CDmodel.gradients.copy()

