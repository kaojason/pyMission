from mission import *

params = {
    'S': 427.8/1e2,
    'ac_w': 210000*9.81/1e6,
    'thrust_sl': 1020000.0/1e6,
    'SFCSL': 8.951,
    'AR': 8.68,
    'e': 0.8,
    }

num_elem = 30
num_cp = 10
x_range = 700.0e3

h_init = numpy.ones(num_cp)*8
h_init[0] = 0.0
h_init[-1] = 0.0

M_init = numpy.ones(num_cp)*0.8
x_init = numpy.linspace(0.0, x_range, num_cp)

traj = OptTrajectory(num_elem, num_cp)
traj.set_init_h(h_init)
traj.set_init_M(M_init)
traj.set_init_x(x_init)
traj.set_params(params)
main = traj.initialize()

main.compute(True)

main.vec['du'].array[:] = 0.0

if 1:
    exit()

if 0:
    main.check_derivatives_all2()
    exit()

opt = Optimization(main)
opt.add_design_variable('h_CP', value=h_init, lower=0.0)
opt.add_objective('wf_obj')
opt.add_constraint('h_i', lower=0.0, upper=0.0)
opt.add_constraint('h_f', lower=0.0, upper=0.0)
opt.add_constraint('Tmin', upper=0.0)
opt.add_constraint('Tmax', upper=0.0)
opt('SNOPT')
