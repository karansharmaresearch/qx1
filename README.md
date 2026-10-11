# qx+1 for q = 2^k - 1

code and data for the map n -> (q*n + 1) / 2^v on odd n, with q = 7, 15, 31 up to 1023

code/cycles.c searches for cycles (needs GMP and OpenMP)
code/tree.c finds numbers that reach 1
code/bounds.py gives the bounds on the smallest number in a cycle
the check_ files compare the C programs with python
data has the output of all runs

to build

    gcc -O2 -fopenmp -o code/cycles.exe code/cycles.c -lgmp
    gcc -O2 -fopenmp -o code/tree.exe code/tree.c

then

    python code/check_cycles.py
    python code/run_cycles.py 7 91312
    python code/run_tree.py 7 32

no cycle other than 1 was found with up to this many odd steps

    q      steps
    7      91312
    15     44001
    31     31765
    63     12631
    127    67960
    255    94924
    511    83317
    1023   94354


