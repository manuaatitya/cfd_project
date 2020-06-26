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

        # Initialize the solution matrix with zeros
        self.solution_matrix = np.zeros((200,self.nodes))

        # Based on the start condition
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
        for i in range (0,len(prev_x)):
            if(not(abs(temp_x[i] - prev_x[i]) <= self.tolerance)):
                return False
        return True
    
    # Gauss Seidal Method to solve the system AX = B
    # X is the initial trial solution
    def gauss_seidal(self,A,B,X):
        for i in range(0,len(A)):
            # Temporary variable to store intermediate variables
            temp = B[i]
            for j in range(0,len(A)):
                if(i != j):
                    temp -= A[i][j] * X[j] # temp_x_solution[j]
            # Update value for the next iteration
            X[i] = temp/A[i][i]

        if(not(self.check_x_solution_within_limit(X,self.prev_X))):
            self.prev_X = X
            self.gauss_seidal(A,B,X)
        return X

    # Find the concentration at each grid point
    def generate_solution_of_equations(self):
        self.M = self.M[1 : self.nodes - 1, 1 : self.nodes - 1]
        self.K = self.K[1 : self.nodes - 1, 1 : self.nodes - 1]
        print("Printing K ",self.K.shape)
        for i in range(1,200):
            # temp = self.solution_matrix[i - 1].transpose()
            temp = self.solution_matrix[i-1][1 : self.nodes - 1].transpose()
            B = self.M @ temp - self.dt * (self.K @ temp)
            self.prev_X = np.zeros(len(self.M))
            self.solution_matrix[i][1 : self.nodes - 1] = self.gauss_seidal(self.M,B,np.zeros(len(self.M)))
    
    # Plot the values in a graph
    def plot_concentration(self):
        x_axis = np.linspace(0,1,self.nodes)
        plt.plot(x_axis,self.solution_matrix[5],label = "t =5")

        plt.plot(x_axis,self.solution_matrix[50], label = "t = 50")

        plt.plot(x_axis,self.solution_matrix[80], label = "t = 80")
        plt.legend()
        plt.show()
    
