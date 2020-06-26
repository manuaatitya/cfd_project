##################################################
#                CFD Assignment                  #
##################################################

# Import necessary header files
import numpy as np
import matplotlib.pyplot as plt
from cfd_simulation import simulation


# Main Driver Function for the application
def main():
    # Get the number of elements input from the user
    number_elements = int(input('Enter the number of nodes needed in the domain \t'))
    
    # Get the total time (seconds) input from the user
    evolution_time = int(input('Enter the time (s) required for the evolution simulation \t'))

    timesteps = evolution_time/200

    tolerance = float(input('Enter the tolerance required for the Gauss Siedal Convergence \t'))

    D = float(input("Enter the D value in the unsteady diffusion equation \t"))

    initial_concentration = float(input("Enter the initial concentration \t"))

    print('Starting Simulation ... \n')

    # Simulation start and functions
    simulation1 = simulation(number_elements,evolution_time,timesteps,tolerance,D,initial_concentration)
    print("Simulation parameters initialised")
    simulation1.generate_mass_matrix()
    print("Mass Matrix M generated")
    simulation1.generate_stiffness_matrix()
    print("Stiffness Matrix K generated")
    simulation1.generate_solution_of_equations()
    print("Concentration values at each grid point calculated")

    print(simulation1.solution_matrix)

    simulation1.plot_concentration()

if __name__ == "__main__":
    main()