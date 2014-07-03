from mission import *

params = {
    'S': 427.8/1e2,
    'ac_w': 210000*9.81/1e6,
    'thrust_sl': 1020000.0/1e6,
    'SFCSL': 8.951,
    'AR': 8.68,
    'e': 0.8,
    }

num_elem = 100
num_cp = 30
x_range = 700.0e3

h_init = numpy.ones(num_cp)*8
h_init[0] = 0.0
h_init[-1] = 0.0

M_init = numpy.ones(num_cp)*0.8
x_init = numpy.linspace(0.0, x_range, num_cp)/1e6

traj = OptTrajectory(num_elem, num_cp)
traj.set_init_h(h_init)
traj.set_init_M(M_init)
traj.set_init_x(x_init)
traj.set_params(params)
main = traj.initialize()

main.compute(True)
'''
print
print
main.vec['du'].array[:] = 0.0
main.compute(True)
'''
print
print
#main.compute_derivatives('fwd', 'h_pt', output=True)
#main.compute_derivatives('rev', 'wf_obj', output=True)
#exit()

'''
print
print
main.compute_derivatives('rev', 'wf_obj', output=True)
'''

if 0:
    # print figure option #
    fig = matplotlib.pylab.figure(figsize=(12.0,12.0))
    v = main.vec['u']
    nr, nc = 6, 2
    fig.add_subplot(nr,nc,1).plot(v('x')*1000.0, v('h'))
    fig.add_subplot(nr,nc,1).set_ylabel('Altitude (km)')
    fig.add_subplot(nr,nc,2).plot(v('x')*1000.0, v('v')*1e2)
    fig.add_subplot(nr,nc,2).set_ylabel('Velocity (m/s)')
    fig.add_subplot(nr,nc,3).plot(v('x')*1000.0, v('alpha')*1e-1*180.0/numpy.pi)
    fig.add_subplot(nr,nc,3).set_ylabel('AoA (deg)')
    fig.add_subplot(nr,nc,4).plot(v('x')*1000.0, v('tau'))
    fig.add_subplot(nr,nc,4).set_ylabel('Throttle')
    fig.add_subplot(nr,nc,5).plot(v('x')*1000.0, v('eta')*1e-1*180.0/numpy.pi)
    fig.add_subplot(nr,nc,5).set_ylabel('Trim Angle (deg)')
    fig.add_subplot(nr,nc,6).plot(v('x')*1000.0, v('fuel_w')*1e6/(9.81*0.804))
    fig.add_subplot(nr,nc,6).set_ylabel('Fuel (L)')
    fig.add_subplot(nr,nc,7).plot(v('x')*1000.0, v('rho'))
    fig.add_subplot(nr,nc,7).set_ylabel('rho')
    fig.add_subplot(nr,nc,8).plot(v('x')*1000.0, v('CL_tar'))
    fig.add_subplot(nr,nc,8).set_ylabel('CL_tar')
    fig.add_subplot(nr,nc,9).plot(v('x')*1000.0, v('CD')*0.1)
    fig.add_subplot(nr,nc,9).set_ylabel('CD')
    fig.add_subplot(nr,nc,10).plot(v('x')*1000.0, v('CT_tar')*0.1)
    fig.add_subplot(nr,nc,10).set_ylabel('CT_tar')
    fig.add_subplot(nr,nc,11).plot(v('x')*1000.0, v('gamma')*0.1)
    fig.add_subplot(nr,nc,11).set_ylabel('gamma')
    fig.add_subplot(nr,nc,12).plot(v('x')*1000.0, (v('fuel_w')+v('ac_w'))*1e6/9.81*2.2)
    fig.add_subplot(nr,nc,12).set_ylabel('W (lb)')
    fig.add_subplot(nr,nc,6).set_xlabel('Distance (km)')
    fig.add_subplot(nr,nc,12).set_xlabel('Distance (km)')
    fig.savefig("OptFig.pdf")
    fig.savefig("OptFig.png")

    exit()

if 0:
    # derivatives check #
    main.check_derivatives_all2()
    exit()

opt = Optimization(main)
opt.add_design_variable('h_pt', value=h_init, lower=0.0)
opt.add_objective('wf_obj')
opt.add_constraint('h_i', lower=0.0, upper=0.0)
opt.add_constraint('h_f', lower=0.0, upper=0.0)
opt.add_constraint('Tmin', upper=0.0)
opt.add_constraint('Tmax', upper=0.0)
opt('SNOPT')

# PRINT FIGURE #

fig = matplotlib.pylab.figure(figsize=(12.0,12.0))
v = main.vec['u']
nr, nc = 6, 2
fig.add_subplot(nr,nc,1).plot(v('x')*1000.0, v('h'))
fig.add_subplot(nr,nc,1).set_ylabel('Altitude (km)')
fig.add_subplot(nr,nc,2).plot(v('x')*1000.0, v('v')*1e2)
fig.add_subplot(nr,nc,2).set_ylabel('Velocity (m/s)')
fig.add_subplot(nr,nc,3).plot(v('x')*1000.0, v('alpha')*1e-1*180.0/numpy.pi)
fig.add_subplot(nr,nc,3).set_ylabel('AoA (deg)')
fig.add_subplot(nr,nc,4).plot(v('x')*1000.0, v('tau'))
fig.add_subplot(nr,nc,4).set_ylabel('Throttle')
fig.add_subplot(nr,nc,5).plot(v('x')*1000.0, v('eta')*1e-1*180.0/numpy.pi)
fig.add_subplot(nr,nc,5).set_ylabel('Trim Angle (deg)')
fig.add_subplot(nr,nc,6).plot(v('x')*1000.0, v('fuel_w')*1e6/(9.81*0.804))
fig.add_subplot(nr,nc,6).set_ylabel('Fuel (L)')
fig.add_subplot(nr,nc,7).plot(v('x')*1000.0, v('rho'))
fig.add_subplot(nr,nc,7).set_ylabel('rho')
fig.add_subplot(nr,nc,8).plot(v('x')*1000.0, v('CL_tar'))
fig.add_subplot(nr,nc,8).set_ylabel('CL')
fig.add_subplot(nr,nc,9).plot(v('x')*1000.0, v('CD')*0.1)
fig.add_subplot(nr,nc,9).set_ylabel('CD')
fig.add_subplot(nr,nc,10).plot(v('x')*1000.0, v('CT_tar')*0.1)
fig.add_subplot(nr,nc,10).set_ylabel('CT_tar')
fig.add_subplot(nr,nc,11).plot(v('x')*1000.0, v('gamma')*0.1)
fig.add_subplot(nr,nc,11).set_ylabel('gamma')
fig.add_subplot(nr,nc,12).plot(v('x')*1000.0, (v('fuel_w')+v('ac_w'))*1e6/9.81*2.2)
fig.add_subplot(nr,nc,12).set_ylabel('W (lb)')
fig.add_subplot(nr,nc,6).set_xlabel('Distance (km)')
fig.add_subplot(nr,nc,12).set_xlabel('Distance (km)')
fig.savefig("OptFig.pdf")
fig.savefig("OptFig.png")
