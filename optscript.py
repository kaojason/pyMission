"""
MISSION ANALYSIS/TRAJECTORY OPTIMIZATION
This is the runscript used for the trajectory optimization problem.
For details regarding the setup of the analysis problem, see mission.py
The mission analysis and trajectory optimization tool was developed by:
    Jason Kao*
    John Hwang*

* University of Michigan Department of Aerospace Engineering,
  Multidisciplinary Design Optimization lab
  mdolab.engin.umich.edu
"""


from mission import *
from history import *
import time
from subprocess import call

# USER SPECIFIED DATA

params = {
    'S': 427.8/1e2,
    'ac_w': 210000*9.81/1e6,
    'thrust_sl': 1020000.0/1e6,
    'SFCSL': 40,#8.951,
    'AR': 8.68,
    'e': 0.8,
    }

num_elem = 2000
num_cp = 200
x_range = 150.0
folder_name = '/home/jason/Documents/Results/test-'

# END USER SPECIFIED DATA

v_init = numpy.ones(num_cp)*2.3
#x_init = numpy.linspace(0.0, x_range, num_cp)/1e6
x_init = x_range * 1e3 * (1-numpy.cos(numpy.linspace(0, 1, num_cp)*numpy.pi))/2/1e6

h_init = 1 * numpy.sin(numpy.pi * x_init / (x_range/1e3))

gamma_lb = numpy.tan(-10.0 * (numpy.pi/180.0))/1e-1
gamma_ub = numpy.tan(10.0 * (numpy.pi/180.0))/1e-1

traj = OptTrajectory(num_elem, num_cp)
traj.set_init_h(h_init)
traj.set_init_v(v_init)
traj.set_init_x(x_init)
traj.set_params(params)
traj.set_folder_name(folder_name)
main = traj.initialize()

main.compute(True)

opt = Optimization(main)
opt.add_design_variable('h_pt', value=h_init, lower=0.0, upper=20.0)
opt.add_objective('wf_obj')
opt.add_constraint('h_i', lower=0.0, upper=0.0)
opt.add_constraint('h_f', lower=0.0, upper=0.0)
opt.add_constraint('Tmin', upper=0.0)
opt.add_constraint('Tmax', upper=0.0)
opt.add_constraint('gamma', lower=gamma_lb, upper=gamma_ub, 
                   get_jacs=main('gamma').get_jacs, linear=True)

start = time.time()
opt('SNOPT')
print 'OPTIMIZATION TIME', time.time() - start
main.history.print_max_min(main.vec['u'])
run_case, last_itr = main.history.get_index()

folder_name = folder_name + 'dist'+str(int(x_range))\
    +'km-'+str(num_cp)+'-'+str(num_elem)+'-'+str(run_case)+'/.'
call (["mv", "./SNOPT_print.out", folder_name])

