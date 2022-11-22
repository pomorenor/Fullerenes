import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import complex_ode


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
        return [1j*y[0] +1j*y[1], 1j*y[0] +1j*y[1]]

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

        real_values_fist_sol = [np.abs(i) for i in first_sol_array]
        real_values_second_sol = [np.abs(i) for i in second_sol_array]

            #print(sol.t, sol.y)
        return times_array, real_values_fist_sol, real_values_second_sol


    ### Now we try another approach








electron_in_a_box = Particle_in_a_box(1,1,4)
print(electron_in_a_box.basis_functions())
inicon = electron_in_a_box.set_initial_conditions([0.7,0.3])
solutions = electron_in_a_box.solve_system(10,0.01 ,electron_in_a_box.system_of_equations,inicon)

time_steps = solutions[0]
first_solution = solutions[1]
second_solution = solutions[2]

plt.figure()
plt.plot(time_steps, first_solution)
plt.plot(time_steps, second_solution)
plt.grid()
plt.show()
