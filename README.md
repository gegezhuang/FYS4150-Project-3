# FYS4150 Project 3

# Usage

## C++

To compile the c++ file `make`. This will compile `src/main.cpp` and all other `cpp`-files and `hpp`-files and run the test.

```
Usage
        ./runner [flags]

Options:
        -h      Show this help message
        -t      Run all tests
        -p10    Solves problem 9 in the problems set
        -p10    Solves problem 10 in the problems set
        -fg     Simulates the fine grained resonance with and without coulomb interactions
```

## Python

```
usage: plot.py [-h] [-s] [-r] [-f] [-e] [-p] [-c] [-a]

Get plots for the simulation for the penning trap problem

optional arguments:
  -h, --help      show this help message and exit
  -s, --solution  To plot RK4, FE and analytical solution.
  -r, --rough     To plot rough-grained scan of frequencies for particles left
  -f, --fine      To plot fine-grained scan of frequencies for particles left
  -e, --error     To plot relative error of the Runge-Kutta and Forward-Euler
  -p, --phase     To plot phase portraits
  -c, --converge  To plot convergence rate of the Runge-Kutta and Forward-Euler
  -a, --all       To plot all
```
