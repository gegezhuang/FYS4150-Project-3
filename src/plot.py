import numpy as np
import pandas as pd
import seaborn as sns
from pandas.core.frame import DataFrame
import argparse

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sns.set()


def plot_3d_solution(ax, filename: str, label: str):
    """Plot the curve described in a file on `ax`, in 3d.

    Arguments:
        ax: Axes instance to plot on. Must be Axes3D.
        filename: Name of file to read. See the function `get_solution` for
                  details on format.
        label: Label to put on plot.
    """
    df = pd.read_csv(f"data/{filename}.csv")
    ax.plot(df.x, df.y, df.z, label=label)


def plot_2d_solution_two_particles(filename: str, label: str):
    """Plot the curve described in a file on `ax`, in 3d.

    Arguments:
        filename: Name of file to read. See the function `get_solution` for
                  details on format.
        label: Label to put on plot.
    """
    df = pd.read_csv(f"data/{filename}.csv")
    plt.plot(df.x1, df.y1, label=label + " particle 1")
    plt.plot(df.x2, df.y2, label=label + " particle 2")


def plot_3d_solution_two_particles(ax, filename: str, label: str):
    """Plot the curve described in a file on `ax`, in 3d.

    Arguments:
        ax: Axes instance to plot on. Must be Axes3D.
        filename: Name of file to read. See the function `get_solution` for
                  details on format.
        label: Label to put on plot.
    """
    df = pd.read_csv(f"data/{filename}.csv")
    ax.plot(df.x1, df.y1, df.z1, label=label + " particle 1")
    ax.plot(df.x2, df.y2, df.z2, label=label + " particle 2")


def relative_err(computed: DataFrame, expected: DataFrame):
    return np.sqrt(
        (computed["x"] - expected["x"]) ** 2
        + (computed["y"] - expected["y"]) ** 2
        + (computed["z"] - expected["z"]) ** 2
        / (computed["x"] * 2 + computed["y"] * 2 + computed["z"] ** 2)
    )


def absolute_error(computed: DataFrame, expected: DataFrame):
    return np.sqrt(
        (computed["x"] - expected["x"]) ** 2
        + (computed["y"] - expected["y"]) ** 2
        + (computed["z"] - expected["z"]) ** 2
    )


def plot_error():
    for filename, method_name in zip(
        ["data/rk4h=", "data/feh="], ["Runge Kutta 4", "Forward Euler"]
    ):
        for h in ["1e-4", "5e-2", "1e-1", "5e-1", "1"]:
            df_a = pd.read_csv(f"data/analyticalh={h}.csv")
            df_num = pd.read_csv(f"{filename}{h}.csv")
            epsilon = relative_err(df_num, df_a)

            t = np.linspace(0, 100, len(df_a))
            h = h.replace("e", "\\cdot 10^{") + "}"
            plt.plot(t, epsilon, label=fr"$h = {h}$")
        #  plt.title(f"Relative error of {method_name}")
        plt.ylabel("Error (log)")
        plt.xlabel("Time ($\\mu s$)")
        plt.yscale("log")
        plt.legend()
        plt.savefig(f"data/{method_name.lower().replace(' ', '_')}_relativeError.pdf")
        plt.close()


def plot_convergence_rate():
    for filename, method_name in zip(
        ["data/rk4h=", "data/feh="], ["Runge Kutta 4", "Forward Euler"]
    ):

        # Find max error for each h
        max_delta = []
        h_list = ["1", "5e-1", "1e-1", "5e-2", "1e-4"]
        for h in h_list:
            df_a = pd.read_csv(f"data/analyticalh={h}.csv")
            df_num = pd.read_csv(f"{filename}{h}.csv")
            delta = absolute_error(df_num, df_a)
            max_delta.append(max(delta))

        fh_list = [float(h) for h in h_list]

        r_err = (
            1
            / 4
            * sum(
                [
                    np.log(max_delta[i] / max_delta[i - 1])
                    / np.log(fh_list[i] / fh_list[i - 1])
                    for i in range(1, 5)
                ]
            )
        )

        plt.plot(fh_list, max_delta)
        plt.xlabel("Step size (log, $\\mu s$)")
        plt.gca().invert_xaxis()
        plt.ylabel("Max absolute error (log)")
        plt.xscale("log")
        plt.yscale("log")
        #  plt.title(f"Convergence rate for {method_name} ($r_{'{err}'} = {r_err:.3g}$)")
        plt.savefig(f"data/{method_name.lower().replace(' ', '_')}_convergence.pdf")
        plt.close()


def plot_all_solutions_one_particle():
    fig = plt.figure()
    ax = Axes3D(fig)
    plot_3d_solution(ax, "forward_euler_one_particle", "Forward Euler solution")
    plot_3d_solution(ax, "runge_kutta_positons_one_particle", "Runge-Kutta 4 solution")
    plot_3d_solution(ax, "analytical_positons", "Analytical solution")
    plt.legend()
    #  ax.set_title(
    #      r"Simulation of particle position for $p_0=(100, 0, 100)$, $v_0=(0, 1, 0)$, $\Delta t=10^{-3}$, $N=10^5$"
    #  )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.savefig(f"data/position_estimates.pdf")
    plt.show()


def plot_two_particles_2d():
    plot_2d_solution_two_particles(
        "runge_kutta_positons_two_particles_with_interactions",
        "Runge-Kutta 4 solution with interactions",
    )
    plot_2d_solution_two_particles(
        "runge_kutta_positons_two_particles_without_interactions",
        "Runge-Kutta 4 solution without interactions",
    )
    plt.legend()
    #  plt.title(r"Motion of two particles in the xy-plane with and without iteraction")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(f"data/two_particles_2d.pdf")
    plt.show()


def plot_two_particles_3d():
    fig = plt.figure()
    ax = Axes3D(fig)
    plot_3d_solution_two_particles(
        ax,
        "runge_kutta_positons_two_particles_with_interactions",
        "Runge-Kutta 4 solution with interactions",
    )
    plot_3d_solution_two_particles(
        ax,
        "runge_kutta_positons_two_particles_without_interactions",
        "Runge-Kutta 4 solution without interactions",
    )
    plt.legend()
    #  plt.title(r"Motion of two particles in the xy-plane with and without iteraction")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.savefig(f"data/two_particles_3d.pdf")
    plt.show()


def plot_all_solutions():
    plot_all_solutions_one_particle()
    plot_two_particles_2d()
    plot_two_particles_3d()


def plot_frequencies_rough():
    fig, axs = plt.subplots(3, sharex=True, sharey=True)
    #  fig.suptitle(r"Particles left after $500 \mu s$")
    fig.text(0.5, 0.03, r"$\omega_V \, [MHz]$", ha="center", fontsize="small")
    fig.text(
        0.04, 0.5, "#particles", va="center", rotation="vertical", fontsize="small"
    )
    filenames = [
        "particles_left_f=0.100000.csv",
        "particles_left_f=0.400000.csv",
        "particles_left_f=0.700000.csv",
    ]
    # fig.xlabel(r"$\omega_V \, [MHz]$")
    fs = ["0,1", "0,4", "0,7"]
    for i in range(3):
        df = pd.read_csv(f"data/{filenames[i]}")
        df.columns = df.columns.str.replace(" ", "_")
        axs[i].set_title(f"Amplitude = {fs[i]}")
        axs[i].plot(df.omega_V, df.particles_left, "o", markersize=2)
    plt.savefig(f"data/particles_left_rough_grained.pdf")
    plt.show()


def plot_frequencies_fine():
    filenames = [
        "fine_grained_no_coulomb_interactions.csv",
        "fine_grained_with_coulomb_interactions.csv",
    ]
    legends = [
        "No coulomb interactions",
        "With coulomb interactions",
    ]
    for filename, legend in zip(filenames, legends):
        df = pd.read_csv(f"data/{filename}")
        df.columns = df.columns.str.replace(" ", "_")
        plt.plot(df.omega_V, df.particles_left, "o", markersize=8, label=legend)
    #  plt.title("Fine grained simulation")
    plt.xlabel("frequency")
    plt.ylabel("# of particles left in trap")
    plt.legend()
    plt.savefig(f"data/particles_left_fine_grained.pdf")
    plt.show()


def plot_phase():
    #  _, ax = plt.subplots(nrows=3, ncols=2)

    filenames = [
        "data/runge_kutta_positons_two_particles_without_interactions.csv",
        "data/runge_kutta_positons_two_particles_with_interactions.csv",
    ]
    axises = ["x", "y", "z"]

    for y, axis in enumerate(axises):
        _, ax = plt.subplots(nrows=1, ncols=2)
        for x, filename in enumerate(filenames):
            df = pd.read_csv(filename)
            h = 1e-3
            v1_axis = np.diff(df[f"{axis}1"]) / h
            v2_axis = np.diff(df[f"{axis}2"]) / h
            ax[x].plot(df[f"{axis}1"][1:], v1_axis, label="Particle 1")
            ax[x].plot(df[f"{axis}2"][1:], v2_axis, label="Particle 2")
            ax[x].legend()
            ax[x].set_xlabel("x")
            ax[x].title.set_text(
                #  f"Phase portrait ({axis}, v{axis}) {'without' if 'without' in filename else 'with'} coulomb interactions"
                f"{'Without' if 'without' in filename else 'With'} coulomb interactions"
            )
            ax[x].set_ylabel("y")
            #  plt.savefig("data/phase_test.pdf")
        plt.tight_layout()
        plt.savefig(f"data/phase_plots_two_particles_{axis}.pdf")
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get plots for the simulation for the penning trap problem"
    )
    parser.add_argument(
        "-s",
        "--solution",
        help="To plot RK4, FE and analytical solution.",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--rough",
        help="To plot rough-grained scan of frequencies for particles left",
        action="store_true",
    )
    parser.add_argument(
        "-f",
        "--fine",
        help="To plot fine-grained scan of frequencies for particles left",
        action="store_true",
    )
    parser.add_argument(
        "-e",
        "--error",
        help="To plot relative error of the Runge-Kutta and Forward-Euler",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--phase",
        help="To plot phase portraits",
        action="store_true",
    )
    parser.add_argument(
        "-c",
        "--converge",
        help="To plot convergence rate of the Runge-Kutta and Forward-Euler",
        action="store_true",
    )
    parser.add_argument(
        "-a",
        "--all",
        help="To plot all",
        action="store_true",
    )
    args = parser.parse_args()
    if not any(vars(args).values()):
        parser.print_help()
    if args.solution or args.all:
        plot_all_solutions()
    if args.rough or args.all:
        plot_frequencies_rough()
    if args.fine or args.all:
        plot_frequencies_fine()
    if args.error or args.all:
        plot_error()
    if args.phase or args.all:
        plot_phase()
    if args.converge or args.all:
        plot_convergence_rate()
