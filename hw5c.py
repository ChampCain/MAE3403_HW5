#Champ Cain
#MAE 3403
#HW5c

# I used AI to help with implementing the ODE's. I was initially pushing out negative
# pressures and they never did converge. It helped me identify what I could do
# to help the code converge.


# region imports
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# endregion

# region functions
def ode_system(t, X, *params):
    '''
    The ode system is defined in terms of state variables.
    I have as unknowns:
    x: position of the piston (This is not strictly needed unless I want to know x(t))
    xdot: velocity of the piston
    p1: pressure on right of piston
    p2: pressure on left of the piston
    For initial conditions, we see: x=x0=0, xdot=0, p1=p1_0=p_a, p2=p2_0=p_a
    :param X: The list of state variables.
    :param t: The time for this instance of the function.
    :param params: the list of physical constants for the system.
    :return: The list of derivatives of the state variables.
    '''
    #unpack the parameters
    A, Cd, ps, pa, V, beta, rho, Kvalve, m, y = params

    #state variables
    #conveniently rename the state variables
    x = X[0]
    xdot = X[1]
    p1 = X[2]
    p2 = X[3]

    #calculate derivatives
    xddot = (p1 - p2) * A / m
    p1dot = (y * Kvalve * (ps - p1) - rho * A * xdot) * (beta / (V * rho))
    p2dot = (-y * Kvalve * (p2 - pa) + rho * A * xdot) * (beta / (rho * V))

    # Ensure pressures do not go negative
    if p1 < 0:
        p1 = 0
        p1dot = 0
    if p2 < 0:
        p2 = 0
        p2dot = 0

    #return the list of derivatives of the state variables
    return [xdot, xddot, p1dot, p2dot]

def main():
    #After some trial and error, I found all the action seems to happen in the first 0.02 seconds
    t = np.linspace(0, 0.02, 200)
    #myargs=(A, Cd, Ps, Pa, V, beta, rho, Kvalve, m, y)
    myargs = (4.909E-4, 0.6, 1.4E7, 1.0E5, 1.473E-4, 2.0E9, 850.0, 2.0E-5, 30, 0.002)
    #because the solution calls for x, xdot, p1 and p2, I make these the state variables X[0], X[1], X[2], X[3]
    #ic=[x=0, xdot=0, p1=pa, p2=pa]
    pa = 1.0E5
    ic = [0, 0, pa, pa]
    #call solve_ivp with ode_system as callback
    sln = solve_ivp(ode_system, [0, 0.02], ic, args=myargs, t_eval=t)

    #unpack result into meaningful names
    xvals = sln.y[0]
    xdot = sln.y[1]
    p1 = sln.y[2]
    p2 = sln.y[3]

    #plot the result
    plt.figure(figsize=(10, 8))

    plt.subplot(2, 1, 1)
    plt.plot(t, xdot, 'b-', label='$\dot{x}$')
    plt.title('Velocity of the Piston vs Time')
    plt.xlabel('Time, s')
    plt.ylabel('$\dot{x}$ (m/s)')
    plt.legend(loc='upper right')

    plt.subplot(2, 1, 2)
    plt.plot(t, p1, 'b-', label='$P_1$')
    plt.plot(t, p2, 'r-', label='$P_2$')
    plt.title('Pressures $P_1$ and $P_2$ vs Time')
    plt.xlabel('Time, s')
    plt.ylabel('Pressure (Pa)')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()
# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion