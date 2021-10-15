"""Plot particles movement, and error where we have an analytical solution.
"""

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def get_solution(filename: str) -> np.ndarray:
    """Reads data from data/<filename>.txt.

    Arguments:
        filename: Name of file in `data`-folder, without the .txt file ending. 
                  Data must be on a format that np.read_csv can read, with the 
                  first column being time, and the next three columns being x, y
                  and z positions.

    Return:
        Array with all values.
    """
    return np.read_csv(f"data/{filename}.txt")

def plot_3d_solution(ax: matplotlib.axes, filename: str, label: str):
    """Plot the curve described in a file on `ax`, in 3d.

    Arguments:
        ax: Axes instance to plot on. Must be Axes3D.
        filename: Name of file to read. See the function `get_solution` for 
                  details on format.
        label: Label to put on plot.
    """
    _, x, y, z = get_solution(filename)
    ax.plot(x, y, z, label=label)

def plot_relative_error(ax: matplotlib.axes, analytical_filename: str, 
                        numerical_filename: str, label: str):
    """Plot the total relative error of numerical solution.

    Total relative error is the sum of relative errors for each axis.

    Arguments:
        ax: Axes instance to plot on.
        analytical_filename: Name of file to read, with analytical solution. See
                             the function `get_solution` for details on format.
        numerical_filename: Name of file to read, with numerical solution. See
                            the function `get_solution` for details on format.
        label: Label to put on plot.
    """
    t, analytical_x, analytical_y, analytical_z = get_solution(analytical_filename)
    _, numerical_x, numerical_y, numerical_z = get_solution(numerical_filename)

    total_error = np.sum(
        relative_error(analytical_x, numerical_x)
      + relative_error(analytical_y, numerical_y)
      + relative_error(analytical_z, numerical_z), axis=0)

    plt.plot(t, total_error, label=label)

def relative_error(exact: np.ndarray, approximated: np.ndarray) -> np.ndarray:
    """Computes the relative error.

    Arguments:
        exact: Exact computed.
        approximated: An approximate array.

    Return:
        Relative error.
    """
    return np.abs((exact - approximated)/exact)

def main():
    pass

if __main__:
    main()