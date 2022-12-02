import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import complex_ode
from matplotlib import animation as anim
from numpy import diff
import imageio



def basis_set(x,n,l):
    number_sets = [(i+1) for i in range(0,n)]
    basis_sets = [np.sqrt(2/l)*np.sin((i*np.pi*x)/l) for i in number_sets]
    return basis_sets

def Ex_pos(c_1,c_2, c_mixed,wave_functions, interval):
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

    for i,j,k in zip(c_1,c_2, c_mixed):
        position_expected_value.append([i*np.trapz(first_integrand, interval)+j*np.trapz(first_integrand, interval)+2*np.real(k*np.trapz(interference_integrand,interval))])




    return position_expected_value


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
        return [(-1j/2)*(36/np.pi**2)*y[0] -1j*y[1], -1j*y[0] -(1j/2)*(9/np.pi**2)*y[1]]

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
inicon = electron_in_a_box.set_initial_conditions([1/np.sqrt(6),np.sqrt(5)/np.sqrt(6)])
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

for i,j,k in zip(Probability_density_1, Probability_density_2, Mixed_Probability_density):
       PD.append([sum(x) for x in zip(i,j,k)])

normalized_solution1 = [i/np.sqrt(i) for i in first_solution]
normalized_solution2 = [i/np.sqrt(i) for i in second_solution]
normalized_crossed_solution = [np.real(i)*np.real(i) for i in crossed_solution]

expected_x = Ex_pos(first_solution, second_solution, crossed_solution,basis_set(np.arange(0,6,0.1),2,L),np.arange(0,6,0.1))

fig, axs = plt.subplots()

def animate(frame):

    plt.cla()
    plt.grid()
    plt.xlabel("$\\frac{x}{L}$")
    plt.ylabel("$|\Psi(x,t)|^2 $ ")
    plt.xlim(0, 6)
    plt.ylim(0.0,1.5)

    #plt.fill(np.arange(0,6,0.1), Probability_density_1[frame])
    #plt.fill(np.arange(0,6,0.1), Probability_density_2[frame])
    #plt.fill(np.arange(0,6,0.1), Mixed_Probability_density[frame])

    plt.fill(np.arange(0,6,0.1), PD[frame], color = 'purple')


    plt.tight_layout()

movie = anim.FuncAnimation(fig, animate, frames = 999, interval = 25 ,repeat = True)
movie.save("Dynamics_Particle_in_a_box.gif")
#plt.show()


plt.figure()
plt.plot(time_steps, normalized_solution1, linestyle = 'dashed',label = "$\mathbb{P}_{\psi_1}$")
plt.plot(time_steps, normalized_solution2, linestyle = 'dashed',label = "$\mathbb{P}_{\psi_2}$")
plt.plot(time_steps, normalized_crossed_solution,linestyle = 'dashed', label = "$\mathbb{P}_I$" )
plt.axvline(x = 1.269, linestyle = 'dashed', color = 'black')
plt.axvline(x = 3.82, linestyle = 'dashed' ,color = 'black')
plt.grid()
plt.xlabel("$t(u.a)$")
plt.ylabel("$\mathbb{P}$")
plt.legend(loc = 'upper left')


plt.xlim(0,4)
plt.tight_layout()

plt.savefig("P_vs_t.png")



plt.figure()
plt.plot(time_steps, expected_x, color = 'black')
plt.grid()
plt.xlabel("$t(u.a)$")
plt.ylabel("$\langle x \\rangle_t$")
plt.savefig("position_expected_value.png")
