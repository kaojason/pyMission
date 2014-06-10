from pyMission import *

params = {
    'S': 427.8/1e2,#                  e2
    'Wac': 185953*9.81/1e6*0.75,#     e6
    'cThrustSL': 1020000.0/1e6,#       e6
    'SFCSL': 8.951,#           e-6
    'chord': 8.15,
    'inertia': 4.1,#              e7
    'AR': 8.68,
    'e': 0.8,
    'g': 9.81,
    }

n = 1
numElem = 30
numCP = 10

v_init = 280.0 * numpy.ones(numElem+1)
v_init[0] = 150.0
v_init[-1] = 150.0

missionProblem = Trajectory(1, opt=True)
missionProblem.set_opt(700.0e3, numElem, numCP)
missionProblem.set_ingn_intl(100)
missionProblem.set_params(params)
missionProblem.set_final_Wf(0.0)
missionProblem.set_IC(v_IC=v_init)

problemPtr = missionProblem.initialize()
problemPtr.compute(True)

problemPtr.vec['du'].array[:] = 0.0

h_init =  numpy.ones(numCP)
h_init[0] = 0
h_init[-1] = 0

if 0:
    problemPtr.check_derivatives_all2()
    problemPtr.vec['du'].array[:] = 0.0
    exit()

opt = Optimization(problemPtr)
opt.add_design_variable('h_CP', value=h_init, lower=0.0)
opt.add_objective('Wf_obj')
opt.add_constraint('h_i', lower=0.0, upper=0.0)
opt.add_constraint('h_f', lower=0.0, upper=0.0)
opt.add_constraint('Tmin', upper=0.0)
opt.add_constraint('Tmax', upper=0.0)
opt('SNOPT')

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
    x_ends = problemPtr.vec['u'](['x', 0]) * 1e6
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
    
fig = matplotlib.pylab.figure(figsize=(12.0,10.0))
fig.add_subplot(711).plot(total_x/1000.0, total_h)
fig.add_subplot(711).set_ylabel('Altitude (km)')
fig.add_subplot(712).plot(total_x/1000.0, total_v*1e2)
fig.add_subplot(712).set_ylabel('Velocity (m/s)')
fig.add_subplot(713).plot(total_x/1000.0, total_a*1e-1*180.0/numpy.pi)
fig.add_subplot(713).set_ylabel('AoA (deg)')
fig.add_subplot(714).plot(total_x/1000.0, total_t)
fig.add_subplot(714).set_ylabel('Throttle')
fig.add_subplot(715).plot(total_x/1000.0, total_e*1e-1*180.0/numpy.pi)
fig.add_subplot(715).set_ylabel('Trim Angle (deg)')
fig.add_subplot(716).plot(total_x/1000.0, total_w*1e6/(9.81*0.804))
fig.add_subplot(716).set_ylabel('Fuel (L)')
fig.add_subplot(717).plot(total_x/1000.0, total_p)
fig.add_subplot(717).set_ylabel('rho')
fig.add_subplot(717).set_xlabel('Distance (km)')
fig.savefig("OptFig.pdf")

