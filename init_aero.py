from __future__ import division
from pykriging import sampling, Kriging, RBF
import numpy
import cPickle
import matplotlib.pylab as plt
import MBI, scipy.sparse

def setupGlobalModel(Data, theta0, GE=False, bf_flag=False, corrfunction=0, scaling=True):

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
    #model = RBF()

    if scaling:
        model.setSampleData(x, y, Lower=numpy.array(lb), Upper=numpy.array(ub))
    else:
        model.setSampleData(x, y)
    model.build()

    return model

print 'Loading Tripan data'

print_figure = False

filename = 'AeroData/CRMTransport_4D_CLkriging.pkl'
cPicklefile = open(filename, 'r')

CLData = cPickle.load(cPicklefile)
cPicklefile.close()

filename = 'AeroData/CRMTransport_4D_CMkriging.pkl'
cPicklefile = open(filename, 'r')

CMData = cPickle.load(cPicklefile)
cPicklefile.close()

filename = 'AeroData/CRMTransport_4D_CDGEK.pkl'
cPicklefile = open(filename, 'r')

CDData = cPickle.load(cPicklefile)
cPicklefile.close()

Ndv = 4
theta0 = 10.0 * numpy.ones(Ndv)

CLmodel = setupGlobalModel(CLData, theta0=theta0, GE=False, bf_flag=False)
CMmodel = setupGlobalModel(CMData, theta0=theta0, GE=False, bf_flag=False)
CDmodel = setupGlobalModel(CDData, theta0=theta0, GE=False, bf_flag=True)

print 'Done loading Tripan data'

M_num = 11
a_num = 11
h_num = 11
e_num = 11

M_surr = numpy.array([0.0, 0.4, 0.6, 0.7, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95])
#M_surr = numpy.array([0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95])
a_surr = numpy.linspace(-20, 30, a_num)
#a_surr = numpy.linspace(0, 10, a_num)
#numpy.array([-20, -5, -1, 3, 7, 11, 15, 19, 23, 27])
h_surr = numpy.linspace(-50, 60000, h_num)
#numpy.array([-10, 5000, 10000, 15000, 20000, 25000, 30000,
#                      35000, 40000, 45000])
e_surr = numpy.linspace(-30, 30, e_num)
#numpy.array([-30, -12, -8, -4, 0, 4, 8, 12, 16, 30])
surr_inp = [M_surr, a_surr, h_surr, e_surr]

M_plt = 0.8
a_plt = 0
h_plt = 35980
e_plt = 0

M_sweep = []
a_sweep = []
h_sweep = []
e_sweep = []

count = 0
inputs = numpy.zeros((M_num*a_num*h_num*e_num, 4))
for i in xrange(M_num):
    for j in xrange(a_num):
        for k in xrange(h_num):
            for l in xrange(e_num):
                inputs[count][0] = M_surr[i]
                inputs[count][1] = a_surr[j]
                inputs[count][2] = h_surr[k]
                inputs[count][3] = e_surr[l]

                if print_figure:
                    if ((inputs[count][1] == a_plt) and
                        (inputs[count][2] == h_plt) and
                        (inputs[count][3] == e_plt)):
                        M_sweep.append(count)
                    if ((inputs[count][0] == M_plt) and
                        (inputs[count][2] == h_plt) and
                        (inputs[count][3] == e_plt)):
                        a_sweep.append(count)
                    if ((inputs[count][0] == M_plt) and
                        (inputs[count][1] == a_plt) and
                        (inputs[count][3] == e_plt)):
                        h_sweep.append(count)
                    if ((inputs[count][0] == M_plt) and
                        (inputs[count][1] == a_plt) and
                        (inputs[count][2] == h_plt)):
                        e_sweep.append(count)
                count += 1

print 'Inputs generated'

for model in [CLmodel, CDmodel, CMmodel]:
    model.evaluate(inputs)

print 'Models evaluated'

CL = CLmodel.yhat
CD = CDmodel.yhat
Cm = CMmodel.yhat

print 'Saving data'

numpy.savetxt('/home/jason/Documents/surr_inputs.dat', [M_surr, a_surr, h_surr, e_surr])
numpy.savetxt('/home/jason/Documents/surr_outputs.dat', [CL, CD, Cm])

if print_figure:
    CL_M = numpy.zeros(M_num)
    CL_a = numpy.zeros(a_num)
    CL_h = numpy.zeros(h_num)
    CL_e = numpy.zeros(e_num)

    CD_M = numpy.zeros(M_num)
    CD_a = numpy.zeros(a_num)
    CD_h = numpy.zeros(h_num)
    CD_e = numpy.zeros(e_num)

    CM_M = numpy.zeros(M_num)
    CM_a = numpy.zeros(a_num)
    CM_h = numpy.zeros(h_num)
    CM_e = numpy.zeros(e_num)

    Mach = numpy.zeros(M_num)
    alpha = numpy.zeros(a_num)
    alt = numpy.zeros(h_num)
    eta = numpy.zeros(e_num)

    for index in xrange(M_num):
        CL_M[index] = CL[M_sweep[index]]
        CD_M[index] = CD[M_sweep[index]]
        CM_M[index] = Cm[M_sweep[index]]
        Mach[index] = inputs[M_sweep[index]][0]

    for index in xrange(a_num):
        CL_a[index] = CL[a_sweep[index]]
        CD_a[index] = CD[a_sweep[index]]
        CM_a[index] = Cm[a_sweep[index]]
        alpha[index] = inputs[a_sweep[index]][1]

    for index in xrange(h_num):
        CL_h[index] = CL[h_sweep[index]]
        CD_h[index] = CD[h_sweep[index]]
        CM_h[index] = Cm[h_sweep[index]]
        alt[index] = inputs[h_sweep[index]][2]

    for index in xrange(e_num):
        CL_e[index] = CL[e_sweep[index]]
        CD_e[index] = CD[e_sweep[index]]
        CM_e[index] = Cm[e_sweep[index]]
        eta[index] = inputs[e_sweep[index]][3]

    plt.figure(1)
    plt.plot(Mach, CL_M)
    plt.xlabel('Mach Number')
    plt.ylabel('CL')

    plt.figure(2)
    plt.plot(alpha, CL_a)
    plt.xlabel('Alpha')
    plt.ylabel('CL')

    plt.figure(3)
    plt.plot(alt, CL_h)
    plt.xlabel('Altitude')
    plt.ylabel('CL')

    plt.figure(4)
    plt.plot(eta, CL_e)
    plt.xlabel('Eta')
    plt.ylabel('CL')

    plt.figure(5)
    plt.plot(Mach, CD_M)
    plt.xlabel('Mach Number')
    plt.ylabel('CD')

    plt.figure(6)
    plt.plot(alpha, CD_a)
    plt.xlabel('Alpha')
    plt.ylabel('CD')

    plt.figure(7)
    plt.plot(alt, CD_h)
    plt.xlabel('Altitude')
    plt.ylabel('CD')

    plt.figure(8)
    plt.plot(eta, CD_e)
    plt.xlabel('Eta')
    plt.ylabel('CD')

    plt.figure(9)
    plt.plot(Mach, CM_M)
    plt.xlabel('Mach Number')
    plt.ylabel('CM')

    plt.figure(10)
    plt.plot(alpha, CM_a)
    plt.xlabel('Alpha')
    plt.ylabel('CM')

    plt.figure(11)
    plt.plot(alt, CM_h)
    plt.xlabel('Altitude')
    plt.ylabel('CM')

    plt.figure(12)
    plt.plot(eta, CM_e)
    plt.xlabel('Eta')
    plt.ylabel('CM')
    plt.show()

mbi_CL = numpy.zeros((M_num, a_num, h_num, e_num))
mbi_CD = numpy.zeros((M_num, a_num, h_num, e_num))
mbi_CM = numpy.zeros((M_num, a_num, h_num, e_num))

count = 0
for i in xrange(M_num):
    for j in xrange(a_num):
        for k in xrange(h_num):
            for l in xrange(e_num):
                mbi_CL[i, j, k, l] = CL[count]
                mbi_CD[i, j, k, l] = CD[count]
                mbi_CM[i, j, k, l] = Cm[count]
                count += 1

CL_arr = MBI.MBI(mbi_CL, [M_surr, a_surr, h_surr, e_surr],
                 [M_num, a_num, h_num, e_num], [4, 4, 4, 4])
CD_arr = MBI.MBI(mbi_CD, [M_surr, a_surr, h_surr, e_surr],
                 [M_num, a_num, h_num, e_num], [4, 4, 4, 4])
CM_arr = MBI.MBI(mbi_CM, [M_surr, a_surr, h_surr, e_surr],
                 [M_num, a_num, h_num, e_num], [4, 4, 4, 4])

Mach = 0.80
alpha = 0.0
alt = 35980
eta = 0.0

inputs = numpy.zeros((100, 4))
inputs[:, 0] = numpy.linspace(0.0, 0.95, 100)
inputs[:, 1] = alpha #numpy.linspace(-5.0, 5.0, 100)
inputs[:, 2] = alt #numpy.linspace(0, 40000, 100)
inputs[:, 3] = eta

print 'CL evaluate'
CL_results = CL_arr.evaluate(inputs)
print 'CD evaluate'
CD_results = CD_arr.evaluate(inputs)
print 'CM evaluate'
CM_results = CM_arr.evaluate(inputs)

CL_print = numpy.zeros(100)
CD_print = numpy.zeros(100)
CM_print = numpy.zeros(100)
for i in xrange(100):
    CL_print[i] = CL_results[i, 0]
    CD_print[i] = CD_results[i, 0]
    CM_print[i] = CM_results[i, 0]

print CL_print
print CD_print
print CM_print


plt.figure()
for i in xrange(6):
    alpha = i
    inputs[:, 1] = alpha
    CD_results = CD_arr.evaluate(inputs)
    
    for j in xrange(100):
        CD_print[j] = CD_results[j, 0]

    plt.plot(inputs[:, 0], CD_print, label="a=%d"%(i))
plt.xlabel('Mach Number')
plt.ylabel('CD')
plt.axis([0.0, 1.0, 0.0, 0.18])
plt.legend(loc=2)
plt.show()
