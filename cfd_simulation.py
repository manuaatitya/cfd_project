#include the necessary headers
import numpy as np
import matplotlib.pyplot as plt

# The simulation class
class simulation:

    # Class initialization
    def __init__(self,n,t,dt,tolerance):
        self.elements = n
        # Number of nodes = no. of elements + 1
        self.nodes = n + 1
        self.t = t
        self.dt = dt
        self.tolerance = tolerance

        # The initial guess for the solution is taken as all zeros
        self.X = np.zeros(self.nodes)

        # Initialize the solution matrix with zeros
        self.solution_matrix = np.zeros((self.nodes,self.t / self.dt))

    def check_x_solution_within_limit(self,temp_x):
        for i in range (0,self.nodes):
            if(not(abs(temp_x[i] - self.X[i]) <= self.tolerance)):
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
            temp_x_solution[i] = temp/A[i][i]
        if(not(check_x_solution_within_limit(temp_x_solution))):
            self.X = temp_x_solution
            gauss_seidal(A,B)
        self.X = temp_x_solution
         