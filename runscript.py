from pyMission import *


params = {
    'S': 427.8/1e2,#                  e2
    'Wac': 185953*9.81/1e6,#     e6
    'cThrustSL': 1020000.0/1e6,#       e6
    'SFCSL': 8.951,#           e-6
    'chord': 8.15,
    'inertia': 4.1,#              e7
    'AR': 8.68,
    'e': 0.8,
    'g': 9.81,
    }

n = 1
missionProblem = Trajectory(1)
missionProblem.add_seg_point(0.0,0.0,v=150.0)
#for i in xrange(n):
missionProblem.add_seg_point(8000.0,20.0e3,v=270.0,tau=1.0,numElem=10)
#missionProblem.add_seg_point(11000.0,40.0e3,M=0.82,tau=1.0,numElem=10)
#missionProblem.add_seg_point(15000.0,60.0e3,M=0.82,tau=1.0,numElem=10)
#missionProblem.add_seg_point(8000.0,4500.0,M=0.82,hDot=-15.3,numElem=5)
#missionProblem.add_seg_point(0.0,5000.0,v=150.0,hDot=-15.3,numElem=5)
missionProblem.set_ingn_intl(100)
missionProblem.set_params(params)
missionProblem.set_range(1000.0e3)
missionProblem.set_final_Wf(0.0)
missionProblem.set_IC()

problemPtr = missionProblem.initialize()
#print 'h_IC: ', problemPtr.vec['u']['h',0]
#print 'SFC_IC: ', problemPtr.vec['u']['SFC',0]
#print 'gamma_IC: ', problemPtr.vec['u']['gamma',0]
#print 'Temp_IC: ', problemPtr.vec['u']['Temp',0]
#print 'rho_IC: ', problemPtr.vec['u']['rho',0]

problemPtr.compute(True).array

problemPtr.vec['du'].array[:] = 0.0

print
print 'Computing derivatives'
#print problemPtr.compute_derivatives('fwd', 'h_ends', output=True)#.array
#exit()
'''
problemPtr('segment').kwargs['NL'] = 'NLN_GS'
problemPtr.compute(True).array
'''
problemPtr.check_derivatives_all2()
problemPtr.vec['du'].array[:] = 0.0
exit()
'''
problemPtr('segment').set_mode('fwd', True)
problemPtr('segment').rhs_vec.array[:] = 1.0
problemPtr('segment').solve_dFdu()
'''
#exit()

#problemPtr.compute_derivatives('fwd', 'v_ends', output=True)
#exit()
'''
A_new = problemPtr.compute_derivatives('fwd', 'v_ends', output=True)
print 'ANALYTIC:-------------------------'
print A_new
A = numpy.array(A_new.array)
#print problemPtr.compute_derivatives('fwd', 'S', output=True)
h = 1.0e-3

f0 = numpy.array(problemPtr.vec['u'].array)
print "F0 HERE!"
print problemPtr.vec['u']
problemPtr('v_ends').value[0] += h
#problemPtr('v_ends').value[1] += h
problemPtr.compute(False)
f = numpy.array(problemPtr.vec['u'].array)
print "F HERE!"
print problemPtr.vec['u']
problemPtr.vec['du'].array[:] = (f-f0)/h
#print '--------------------NOW COMPUTING DERIVATIVES'
B = numpy.array(problemPtr.vec['du'].array)
print 'FD:-------------------------------'
#print B
print problemPtr.vec['du']
print numpy.vstack([A,B]).T
#print A
exit()
'''
'''
analytic_array = problemPtr.compute_derivatives('fwd','SFCSL').array[:]
problemPtr.vec['u'](['SFCSL',0])[:] += 1e-5
new_array = problemPtr.compute(True).array
diff_array = (new_array-old_array)/1e-10
print analytic_array
print diff_array
print analytic_array-diff_array
'''
#problemPtr.vec['u'].array[:] *= problemPtr.vec['u0'].array[:]
#print problemPtr.vec['u']
#print problemPtr.vec['u0']
'''
print "---------------------------------------"
for name in ['SFC', 'gamma', 'Temp', 'h', 'rho', 'CT', 'v', 'tau', 'CL',
             'alpha', 'CD', 'Wf', 'CM', 'eta']:
    result = problemPtr([name,0]).check_derivatives(problemPtr.variables.keys())
    print 'check ' + name
    print result[0]
    print result[1]
'''

totalElem = 0
total_x = numpy.array([])
total_h = numpy.array([])
total_v = numpy.array([])
total_a = numpy.array([])
total_t = numpy.array([])
total_e = numpy.array([])
total_w = numpy.array([])
total_p = numpy.array([])
for i in xrange(missionProblem.numSeg):
    totalElem += missionProblem.numElem[i]
    x_ends = problemPtr.vec['u'](['x', 0]) #* 1e6
    x_int = numpy.linspace(x_ends[i], x_ends[i+1], missionProblem.numElem[i]+1)
    total_x = numpy.append(total_x, x_int)
    total_h = numpy.append(total_h, problemPtr.vec['u'](['h', i]))
    total_v = numpy.append(total_v, problemPtr.vec['u'](['v', i]))
    total_a = numpy.append(total_a, problemPtr.vec['u'](['alpha', i]))
    total_t = numpy.append(total_t, problemPtr.vec['u'](['tau', i]))
    total_e = numpy.append(total_e, problemPtr.vec['u'](['eta', i]))
    total_w = numpy.append(total_w, problemPtr.vec['u'](['Wf', i]))
    total_p = numpy.append(total_p, problemPtr.vec['u'](['rho', i]))
    numpy.delete(total_x, totalElem)
    numpy.delete(total_h, totalElem)
    numpy.delete(total_v, totalElem)
    numpy.delete(total_a, totalElem)
    numpy.delete(total_t, totalElem)
    numpy.delete(total_e, totalElem)
    numpy.delete(total_w, totalElem)
    numpy.delete(total_p, totalElem)
numSeg = missionProblem.numSeg
numpy.append(total_x, numpy.array([x_ends[numSeg]]))
numpy.append(total_h, numpy.array([missionProblem.hPts[numSeg]]))
numpy.append(total_v, numpy.array([problemPtr.vec['u'](['v', numSeg-1])[-1]]))
numpy.append(total_a, numpy.array([problemPtr.vec['u'](['alpha', numSeg-1])[-1]]))
numpy.append(total_t, numpy.array([problemPtr.vec['u'](['tau', numSeg-1])[-1]]))
numpy.append(total_e, numpy.array([problemPtr.vec['u'](['eta', numSeg-1])[-1]]))
numpy.append(total_w, numpy.array([problemPtr.vec['u'](['Wf', numSeg-1])[-1]]))
numpy.append(total_p, numpy.array([problemPtr.vec['u'](['rho', numSeg-1])[-1]]))
    
pylab.figure
pylab.subplot(711)
pylab.plot(total_x/1000.0, total_h)
pylab.ylabel('Altitude (km)')
pylab.subplot(712)
pylab.plot(total_x/1000.0, total_v*1e2)
pylab.ylabel('Velocity (m/s)')
pylab.subplot(713)
pylab.plot(total_x/1000.0, total_a*1e-1*180.0/numpy.pi)
pylab.ylabel('AoA (deg)')
pylab.subplot(714)
pylab.plot(total_x/1000.0, total_t)
pylab.ylabel('Throttle')
pylab.subplot(715)
pylab.plot(total_x/1000.0, total_e*1e-1*180.0/numpy.pi)
pylab.ylabel('Trim Angle (deg)')
pylab.subplot(716)
pylab.plot(total_x/1000.0, total_w*1e6/(9.81*0.804))
pylab.ylabel('Fuel (L)')
pylab.subplot(717)
pylab.plot(total_x/1000.0, total_p)
pylab.ylabel('rho')
pylab.xlabel('Distance (km)')
pylab.show()

'''
name = 'v'
result = problemPtr(['CL',0]).check_derivatives()
print 'checking CL-' + name
print result[0]
print result[1]
'''
"""
print problemPtr.vec['u']['h',0]
print problemPtr.vec['p']['SFC',0]['h',0]

main = problemPtr
if main(['h2',0]).comm is not None:
    print 'h2-derv:', main(['h2',0]).check_derivatives('fwd')
if main(['SFC',0]).comm is not None:
    print 'SFC-derv:'
    print main(['SFC',0]).check_derivatives('fwd')
"""
