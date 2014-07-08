from mission import *
import time

params = {
    'S': 427.8/1e2,
    'ac_w': 210000*9.81/1e6,
    'thrust_sl': 1020000.0/1e6/3,
    'SFCSL': 8.951,
    'AR': 8.68,
    'e': 0.8,
    }

num_elem = 1000
num_cp = 300
x_range = 7000.0e3

h_init = numpy.ones(num_cp)*0.5
h_init[0] = 0.0
h_init[-1] = 0.0

v_init = numpy.ones(num_cp)*2.3
x_init = numpy.linspace(0.0, x_range, num_cp)/1e6

traj = OptTrajectory(num_elem, num_cp)
traj.set_init_h(h_init)
traj.set_init_v(v_init)
traj.set_init_x(x_init)
traj.set_params(params)
main = traj.initialize()

main.compute(True)

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
start = time.time()
opt('SNOPT')
print 'OPTIMIZATION TIME', time.time() - start

# PRINT FIGURE #
if 0:
    total_x = main.vec['u']('x')
    total_h = main.vec['u']('h')
    total_v = main.vec['u']('v')
    total_a = main.vec['u']('alpha')
    total_t = main.vec['u']('tau')
    total_e = main.vec['u']('eta')
    total_w = main.vec['u']('fuel_w')
    dist = len(main.vec['u']('x'))
    temp_arr = numpy.zeros((dist,2))
    for i in xrange(dist):
        temp_arr[i] = [total_x[i]/1e3, total_h[i]]
    numpy.savetxt('figure_1_h.dat', temp_arr)
    for i in xrange(dist):
        temp_arr[i] = [total_x[i]/1e3, total_v[i]*1e2]
    numpy.savetxt('figure_2_v.dat', temp_arr)
    for i in xrange(dist):
        temp_arr[i] = [total_x[i]/1e3, total_a[i]*1e-1*180.0/numpy.pi]
    numpy.savetxt('figure_3_a.dat', temp_arr)
    for i in xrange(dist):
        temp_arr[i] = [total_x[i]/1e3, total_t[i]]
    numpy.savetxt('figure_4_t.dat', temp_arr)
    for i in xrange(dist):
        temp_arr[i] = [total_x[i]/1e3, total_e[i]*1e-1*180.0/numpy.pi]
    numpy.savetxt('figure_5_e.dat', temp_arr)
    for i in xrange(dist):
        temp_arr[i] = [total_x[i]/1e3, total_w[i]*1e6/(9.81*0.804)]
    numpy.savetxt('figure_6_f.dat', temp_arr)

fig = matplotlib.pylab.figure(figsize=(12.0,7.0))
v = main.vec['u']
nr, nc = 3, 2
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
fig.add_subplot(nr,nc,5).set_xlabel('Distance (km)')
fig.add_subplot(nr,nc,6).set_xlabel('Distance (km)')
fig.savefig("./OptFig.pdf")
fig.savefig("./OptFig.png")
