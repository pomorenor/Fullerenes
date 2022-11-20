import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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

    def system_of_equations(self, y):
        for i in range(0, self.num_basis_elements):
            






electron_in_a_box = Particle_in_a_box(1,1,4)
print(electron_in_a_box.basis_functions())
