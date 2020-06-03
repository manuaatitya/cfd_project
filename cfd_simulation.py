#include the necessary headers
import numpy as np
import matplotlib.pyplot as plt

# The simulation class
class simulation:

    # Class initialization
    def __init__(self,n,t,dt,tolerance,D,initial_concentration):
        self.elements = n
        # Number of nodes = no. of elements + 1
        self.nodes = n + 1
        self.l_e = 1 / self.elements
        self.t = t
        self.dt = dt
        self.tolerance = tolerance
        self.D = D
        self.initial_concentration = initial_concentration

        # The initial guess for the solution is taken as all zeros
        self.X = np.zeros(self.nodes)

        # Initialize the solution matrix with zeros
        self.solution_matrix = np.zeros((200,self.nodes))

        # Based on the initial condition
        self.solution_matrix[0] = self.initial_concentration * np.ones(self.nodes)

    def generate_stiffness_matrix(self):
        K = np.zeros((self.nodes,self.nodes))
        for i in range (0,self.nodes):
            if(i == 0):
                K[i][i] = 1
                K[i][i + 1] = -1
            elif(i == self.nodes - 1):
                K[i][i] = 1
                K[i][i - 1] = -1
            else:
                K[i][i] = 2
                K[i][i + 1] = -1
                K[i][i - 1] = -1
        K = (self.D / self.l_e) * K
        self.K = K

    def generate_mass_matrix(self):
        M = np.zeros((self.nodes,self.nodes))
        for i in range (0,self.nodes):
            if(i == 0):
                M[i][i] = 2
                M[i][i + 1] = 1
            elif(i == self.nodes - 1):
                M[i][i] = 2
                M[i][i - 1] = 1
            else:
                M[i][i] = 4
                M[i][i + 1] = 1
                M[i][i - 1] = 1
        M = (self.l_e / 6) * M
        self.M = M

    # Function to check if the solution is within the specified limits
    def check_x_solution_within_limit(self,temp_x,prev_x):
        for i in range (0,self.nodes):
            if(not(abs(temp_x[i] - prev_x[i]) <= self.tolerance)):
                return False
        return True
    
    # Gauss Seidal Method to solve the system AX = B
    # X is the initial trial solution
    def gauss_seidal(self,A,B):
        temp_x_solution = np.zeros(self.nodes)
        for i in range(0,self.nodes):
            temp = B[i]
            for j in range(0,self.nodes):
                if(i != j):
                    temp -= A[i][j] * temp_x_solution[j]
            # Update value for the next iteration
            temp_x_solution[i] = temp/A[i][i]

        if(not(self.check_x_solution_within_limit(temp_x_solution,self.prev_X))):
            self.prev_X = temp_x_solution
            self.gauss_seidal(A,B)
        return temp_x_solution

    # Find the concentration at each grid point
    def generate_solution_of_equations(self):
        for i in range(1,self.nodes):
            self.prev_X = np.zeros(self.nodes)
            temp = self.solution_matrix[i - 1].transpose()
            B = self.M @ temp - self.dt * (self.K @ temp)
            self.solution_matrix[i] = self.gauss_seidal(self.M,B)
    # Print concentration values
    def printsolnmatrix(self):       
        for row in self.solution_matrix:
              for val in row:
                  print (val)
              print
