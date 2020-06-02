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
    number_elements = int(input('Enter the number of nodes needed in the domain \n'))
    
    # Get the total time (seconds) input from the user
    evolution_time = int(input('Enter the time (s) required for the evolution simulation \n'))

    timesteps = evolution_time/200

    tolerance = float(input('Enter the tolerance required for the Gauss Siedal Convergence'))

    simulation1 = simulation(number_elements,evolution_time,timesteps,tolerance)


if __name__ == "__main__":
    main()