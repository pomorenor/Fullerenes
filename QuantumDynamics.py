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

def expected_values(c_1,c_2, c_mixed,wave_functions, interval):
    first_integrand = []
    second_integrand = []
    interference_integrand = []
    position_expected_value = []

    for i,j in zip(interval,wave_functions[0]):
        first_integrand.append(j*i*j)

    for i,j in zip(interval,wave_functions[1]):
        second_integrand.append(j*i*j)

    for i,j,k in zip(interval, wave_functions[0], wave_functions[1]):
        interference_integrand.append(j*i*k)

    #for i,j,k in zip(c_1,c_2, c_mixed):
    #    position_expected_value.append([i*np.trapz(l, interval) + j*np.trapz(m,interval) + np.re(k)*np.trapz(n,interval) for l,m,n in zip(first_integrand, second_integrand, interference_integrand)])

    for i,j,k in zip(c_1,c_2,c_mixed):
        position_expected_value.append([i*np.trapz(first_integrand, interval)+j*np.trapz(first_integrand, interval)+np.real(k)*np.trapz(interference_integrand,interval)])


    return position_expected_value

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
        crossed_terms = [np.conjugate(i)*j for i,j in zip(first_sol_array, second_sol_array)]

            #print(sol.t, sol.y)
        return times_array, real_values_fist_sol, real_values_second_sol, crossed_terms










electron_in_a_box = Particle_in_a_box(1,1,4)
#print(electron_in_a_box.basis_functions())
inicon = electron_in_a_box.set_initial_conditions([0.7,0.3])
solutions = electron_in_a_box.solve_system(10,0.01 ,electron_in_a_box.system_of_equations,inicon)

time_steps = solutions[0]
first_solution = solutions[1]
second_solution = solutions[2]
crossed_solution = solutions[3]


L = 6
first_space_wave_function = basis_set(np.arange(0,6,0.1),2,L)[0]
second_space_wave_function = basis_set(np.arange(0,6,0.1),2,L)[1]


Probability_density_1 = []
Probability_density_2 = []
Mixed_Probability_density = []
PD = []

for i in first_solution:
    Probability_density_1.append([i*j*j for j in first_space_wave_function ])

for i in second_solution:
    Probability_density_2.append([i*j*j for j in second_space_wave_function ])

for i in crossed_solution:
    Mixed_Probability_density.append([2*np.real(i*j*j*k*k) for j,k in zip(first_space_wave_function, second_space_wave_function)])

Expected_pos = expected_values(first_solution, second_solution, crossed_solution,basis_set(np.arange(0,6,0.1),2,L),np.arange(0,6,0.1) )


fig, axs = plt.subplots()

def animate(frame):

    plt.cla()
    plt.grid()
    plt.xlabel("Distancia (u.a.)")
    plt.ylabel("Densidad de probabilidad")
    plt.xlim(0, 6)
    plt.ylim(0.0,0.4)

    plt.fill(np.arange(0,6,0.1), Probability_density_1[frame])
    plt.fill(np.arange(0,6,0.1), Probability_density_2[frame])
    plt.fill(np.arange(0,6,0.1), Mixed_Probability_density[frame])


    plt.tight_layout()

movie = anim.FuncAnimation(fig, animate, frames = 999, interval = 25 ,repeat = True)
movie.save("interference.gif")
#plt.show()


plt.figure()
plt.plot(time_steps, Expected_pos)
plt.grid()
plt.xlabel("Time")
plt.ylabel("Expected_pos")
plt.savefig("position_expected_value.png")
