#Champ Cain
#MAE 3403
#HW5

# I used AI to help me identify which iterative process would me the graph.
# I initially used fsolve, but it was getting messed up as it was being divided by
# large numbers when the flow had a very high Re Number. So I ended up using the
# Swamee-Jain Approximation method.

# region imports
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
# endregion

# region functions
def ff(Re, rr, CBEQN=False):
    """
    This function calculates the friction factor for a pipe based on the
    notion of laminar, turbulent, and transitional flow.
    :param Re: the Reynolds number under question.
    :param rr: the relative pipe roughness (expect between 0 and 0.05)
    :param CBEQN: boolean to indicate if I should use Colebrook (True) or laminar equation
    :return: the (Darcy) friction factor
    """
    if CBEQN:
        # Swamee-Jain approximation for the Colebrook equation
        f = 0.25 / (np.log10((rr / 3.7) + (5.74 / Re**0.9)))**2
        return f
    else:
        return 64 / Re

def plotMoody(plotPoint=False, pt=(0, 0), marker='o'):
    """
    This function produces the Moody diagram for a Re range from 1 to 10^8 and
    for relative roughness from 0 to 0.05 (20 steps). The laminar region is described
    by the simple relationship of f=64/Re whereas the turbulent region is described by
    the Colebrook equation.
    :param plotPoint: Whether to plot a point on the Moody diagram.
    :param pt: The (Re, f) coordinates of the point to plot.
    :param marker: The marker style for the point (e.g., 'o' for circle, '^' for triangle).
    :return: just shows the plot, nothing returned
    """
    # Step 1: Create logspace arrays for ranges of Re
    ReValsCB = np.logspace(np.log10(4000), np.log10(1e8), 100)  # Turbulent range
    ReValsL = np.logspace(np.log10(600.0), np.log10(2000.0), 20)  # Laminar range
    ReValsTrans = np.logspace(np.log10(2000.0), np.log10(4000.0), 20)  # Transition range

    # Step 2: Create array for range of relative roughnesses
    rrVals = np.array([0, 1E-6, 5E-6, 1E-5, 5E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4, 1E-3, 2E-3, 4E-3, 6E-3, 8E-3, 1.5E-2, 2E-2, 3E-2, 4E-2, 5E-2])

    # Step 3: Calculate the friction factor in the laminar range
    ffLam = np.array([ff(Re, 0, CBEQN=False) for Re in ReValsL])
    ffTrans = np.array([ff(Re, 0, CBEQN=False) for Re in ReValsTrans])

    # Step 4: Calculate friction factor values for each rr at each Re for turbulent range
    ffCB = np.array([[ff(Re, rr, CBEQN=True) for Re in ReValsCB] for rr in rrVals])

    # Step 5: Construct the plot
    plt.loglog(ReValsL, ffLam, 'b-', label='Laminar')  # Laminar part as a solid line
    plt.loglog(ReValsTrans, ffTrans, 'b--', label='Transition')  # Transition part as a dashed line
    for nRelR in range(len(ffCB)):
        plt.loglog(ReValsCB, ffCB[nRelR], color='k', label=None)  # Turbulent region without label
        plt.annotate(xy=(1e8, ffCB[nRelR][-1]), text=f'{rrVals[nRelR]}', fontsize=8)  # Label at the end of each curve

    plt.xlim(600, 1e8)
    plt.ylim(0.008, 0.10)
    plt.xlabel(r"Reynolds number $Re$", fontsize=16)
    plt.ylabel(r"Friction factor $f$", fontsize=16)
    plt.text(2.5e8, 0.02, r"Relative roughness $\frac{\epsilon}{d}$", rotation=90, fontsize=16)
    ax = plt.gca()  # Capture the current axes for use in modifying ticks, grids, etc.
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=12)  # Format tick marks
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)
    ax.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
    plt.grid(which='both')
    if plotPoint:
        plt.plot(pt[0], pt[1], marker=marker, markersize=12, markeredgecolor='red', markerfacecolor='none')

    plt.legend()  # Show legend with only "Laminar" and "Transition"
    plt.show()

def main():
    plotMoody()
# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion