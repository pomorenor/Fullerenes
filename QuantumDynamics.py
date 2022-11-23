import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import complex_ode
from matplotlib import animation as anim
import imageio



def basis_set(x,n,l):
    number_sets = [(i+1) for i in range(0,n)]
    basis_sets = [np.sqrt(1/l)*np.sin((i*np.pi*x)/l) for i in number_sets]
    return basis_sets

## space wavefunction

# Esto lo vemos mañana
#def sp_wave_function(basis_set, time_coefficients):
#    Psi = 0
#    for i,j in zip(basis_set, time_coefficients):
#        Psi += basis_set[i]*time_coefficients
#    return Psi


class Particle_in_a_box:

    def __init__(self, mass, box_length, num_basis_elements):
        self.mass = mass
        self.box_length = box_length
        self.num_basis_elements = num_basis_elements

    def basis_functions(self):
        n = [i for i in range(0, self.num_basis_elements)]
        return n

    def set_initial_conditions(self, initial_conditions):
        self.initial_conditions = initial_conditions
        return initial_conditions

    def Hamiltonian(self, hamilton):
        self.hamilton = hamilton
        return hamilton

    def system_of_equations(self, t, y):
        return [-1j*y[0] -1j*y[1], -1j*y[0] -1j*y[1]]

    def solve_system(self, time_interval, dt, coeff_equations, ini_condi):
        self.time_interval = time_interval
        self.dt = dt
        times_array = []
        first_sol_array = []
        second_sol_array = []
        #t = np.linspace(0,time_interval,101)
        #sol = odeint(coeff_equations, ini_condi, t)
        sol = complex_ode(coeff_equations, jac = None)
        sol.set_initial_value(ini_condi, 0)

        tf = time_interval
        delta_t = dt
        while sol.t < tf:
            sol.integrate(sol.t+delta_t)
            times_array.append(sol.t)
            first_sol_array.append(sol.y[0])
            second_sol_array.append(sol.y[1])

        real_values_fist_sol = [np.abs(i)**2 for i in first_sol_array]
        real_values_second_sol = [np.abs(i)**2 for i in second_sol_array]
        crossed_terms = [np.abs(i*j)**2 for i,j in zip(first_sol_array, second_sol_array)]

            #print(sol.t, sol.y)
        return times_array, real_values_fist_sol, real_values_second_sol










electron_in_a_box = Particle_in_a_box(1,1,4)
print(electron_in_a_box.basis_functions())
inicon = electron_in_a_box.set_initial_conditions([0.7,0.3])
solutions = electron_in_a_box.solve_system(10,0.01 ,electron_in_a_box.system_of_equations,inicon)

time_steps = solutions[0]
first_solution = solutions[1]
second_solution = solutions[2]

L = 6
first_space_wave_function = basis_set(np.arange(0,6,0.1),2,L)[0]
second_space_wave_function = basis_set(np.arange(0,6,0.1),2,L)[1]


Probability_density_1 = []
Probability_density_2 = []
PD = []

for i in first_solution:
    Probability_density_1.append([i*j for j in first_space_wave_function ])

for i in second_solution:
    Probability_density_2.append([i*j*j for j in second_space_wave_function ])

for i,j in  zip(Probability_density_1, Probability_density_2):
    PD.append(i+j)
#print(first_part_of_pd)


#plt.figure()
#plt.plot(time_steps, first_solution)
#plt.plot(time_steps, second_solution)
#plt.grid()
#plt.show()

#plt.figure()
#plt.plot(np.arange(0,6,0.1), basis_set(np.arange(0,6,0.1),2,6)[0])
#plt.plot(np.arange(0,6,0.1), basis_set(np.arange(0,6,0.1),2,6)[1])

#plt.show()

#plt.figure()
#plt.grid()
#plt.plot(np.arange(0,6,0.1),first_part_of_pd)
#plt.show()


fig, axs = plt.subplots()

def animate(frame):

    plt.cla()
    plt.grid()
    plt.xlim(0, 6)
    plt.ylim(0.0,0.4)
    plt.plot(np.arange(0,6,0.1), Probability_density_2[frame])


    plt.tight_layout()

movie = anim.FuncAnimation(fig, animate, frames = 999, interval = 25 ,repeat = False)
plt.show()
