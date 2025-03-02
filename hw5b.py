# region imports
import hw5a as pta
import random as rnd
import numpy as np
from matplotlib import pyplot as plt
# endregion

# region functions
def ffPoint(Re, rr):
    """
    This function takes Re and rr as parameters and outputs a friction factor according to the following:
    1.  if Re>4000 use Colebrook Equation
    2.  if Re<2000 use f=64/Re
    3.  else calculate a probabilistic friction factor where the distribution has a mean midway between the prediction
        of the f=64/Re and Colebrook Equations and a standard deviation of 20% of this mean
    :param Re:  the Reynolds number
    :param rr:  the relative roughness
    :return:  the friction factor
    """
    if Re >= 4000:
        return pta.ff(Re, rr, CBEQN=True)  # Use Colebrook equation for turbulent flow
    elif Re <= 2000:
        return pta.ff(Re, rr)  # Use laminar equation for laminar flow
    else:
        # Transitional flow: interpolate between laminar and turbulent predictions
        CBff = pta.ff(Re, rr, CBEQN=True)  # Prediction of Colebrook Equation in Transition region
        Lamff = pta.ff(Re, rr)  # Prediction of Laminar Equation in Transition region
        mean = Lamff + (CBff - Lamff) * (Re - 2000) / 2000  # Linear interpolation
        sig = 0.2 * mean  # Standard deviation is 20% of the mean
        return rnd.normalvariate(mean, sig)  # Randomly select a value from a normal distribution

def PlotPoint(Re, f):
    """
    This function plots a point on the Moody diagram with an upward triangle if the flow is transitional
    or a circle otherwise.
    :param Re: Reynolds number
    :param f: Friction factor
    """
    if 2000 < Re < 4000:
        marker = '^'  # Upward triangle for transitional flow
    else:
        marker = 'o'  # Circle for laminar or turbulent flow
    pta.plotMoody(plotPoint=True, pt=(Re, f), marker=marker)

def main():
    """
    Main function to solicit user input, calculate head loss, and plot the Moody diagram.
    """
    while True:
        # Get user input
        D = float(input("Enter the pipe diameter in inches: "))  # Pipe diameter in inches
        epsilon = float(input("Enter the pipe roughness in micro-inches: "))  # Pipe roughness in micro-inches
        Q = float(input("Enter the flow rate in gallons per minute: "))  # Flow rate in gallons per minute

        # Convert inputs to consistent units
        D_ft = D / 12  # Convert diameter to feet
        epsilon_ft = epsilon * 1e-6 / 12  # Convert roughness to feet
        Q_ft3s = Q * 0.002228  # Convert flow rate to cubic feet per second

        # Calculate Reynolds number (Re) and relative roughness (rr)
        nu = 1.08e-5  # Kinematic viscosity of water at 60°F in ft^2/s
        V = Q_ft3s / (np.pi * (D_ft / 2) ** 2)  # Velocity in ft/s
        Re = V * D_ft / nu  # Reynolds number
        rr = epsilon_ft / D_ft  # Relative roughness

        # Calculate friction factor (f)
        f = ffPoint(Re, rr)

        # Calculate head loss per foot (hf/L)
        g = 32.2  # Acceleration due to gravity in ft/s^2
        hf_L = f * (V ** 2) / (2 * g * D_ft)  # Head loss per foot

        # Display results
        print(f"Reynolds number (Re): {Re:.2f}")
        print(f"Relative roughness (ϵ/D): {rr:.6f}")
        print(f"Friction factor (f): {f:.4f}")
        print(f"Head loss per foot (hf/L): {hf_L:.6f} ft/ft")

        # Plot the point on the Moody diagram
        PlotPoint(Re, f)

        # Ask the user if they want to re-specify parameters
        repeat = input("Do you want to re-specify parameters? (yes/no): ").strip().lower()
        if repeat != 'yes':
            break

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion