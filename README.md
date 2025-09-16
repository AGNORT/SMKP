[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An efficient exact algorithm for the submodular knapsack problem

This repository accompanies the paper submitted to the [INFORMS Journal on Computing (IJOC)](https://pubsonline.informs.org/journal/ijoc). All source code is released under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on the paper *paper DOI link* by Feilong Wang, Filippo Ranza, Alice Raffaele, Roberto Roberti, and Renata Mansini.

## 1. Citation

If you use this repository or the accompanying paper in your work, please cite both the paper and this repository using their respective DOIs.

<!-- TODO: Add DOI of the paper -->
<!-- TODO: Add DOI of this repository -->

BibTeX for citing this snapshot of the repository:

<!-- TODO: Insert BibTeX here -->

## 2. Requirements

The algorithms are compiled and executed in a Linux (Ubuntu) environment with **GCC 13.3**, **CMake 3.28**, and **Gurobi 12.0.1**. Please ensure that the following dependencies are installed before compilation:

```bash
GCC     # version 12 or higher
CMake   # version 3.21 or higher
Gurobi  # version 12.0 or higher
```

## 3. Description

This repository provides the implementation of our algorithms to reproduce the computational results reported in the paper. The software allows users to solve both the **Submodular Knapsack Problem (SMKP)** and the **Submodular Bin Packing Problem (SMBP)**

## 3.1 Building

The code is designed for **Linux** environments. 

1. Clone the repository:

```bash
git clone <repo-link>
cd <repo-folder>
```

2. Build the code using **CMake**:

```bash
mkdir build && cd build
cmake ..
make
```

This will generate the executable `SMKMP` in the `build` directory. Since we use Gurobi as the commercial solver, the users need to change the path related to the Gurobi solver in the CMakeLists.txt.

```makefile
# Set Gurobi paths (modify according to your installation)
set(GUROBI_HOME /path/to/your/gurobi/library)
set(GRB_LICENSE_FILE /path/to/your/gurobi/license)
target_link_libraries(SMKMP PRIVATE gurobi_c++ gurobi120 pthread) # Replace gurobi120 with your installed version (e.g., gurobi110 )version
```



## 3.2 Running the Software

 The executable requires command-line arguments. There are **8 parameters in total**, three mandatory and five optional:

```makefile
-i   <file>   # input .txt file (mandatory)
-o   <file>   # output .txt file (mandatory)
-m   <Gurobi>   # solution method (optional, default: DP for SMKP)
```

The "-i" and "-o" are mandatory parameters, "-m" is a optional parameter. The users can also use the "runAll.sh" shell file in the *build* folder to run all instances.

## 3.3 Example Commands

Solve SMKP (Ryn Park. (2021) [[1]](#ref1)  instances):

```bash
cd build
./SMKP -i path/to/input.txt -o path/to/output.txt
```

Solve SMKP with Gurobi:

```bash
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -m Gurobi
```



## 3.4 Example outputs

**SMKP with DP**

```
Solve SMKP use DP!

heuristic labeling get lower bound: 22256
heuristic labeling use time: 0.65s
The best solution obtained by my labelsetting algorithm is: 22256
The solution obtained by TS heuristic is: 22256
The solution obtained by heuristic labeling algorithm is: 22256
The time for my labelsetting algorithm is: 1.453s
The time for heuristic primal solution is: 0.018s
The time for linear relaxation (dual bound) is: 0.121s
The time for heuristic DP is: 0.65s
The time for exact dual bound is: 0s
The time for final exact DP is: 0.646s
The total number of generated labels is: 30982270
The number of non-dominated labels is: 4
The number of CB fathomed labels is: 1254
The number of dominated labels is: 15488946
```
**SMKP with Gurobi**

```
Solve SMKP use Gurobi!
Optimal solution found (tolerance 1.00e-06)
Best objective 2.225600000000e+04, best bound 2.225600000000e+04, gap 0.0000%
***************The results obtained by Gurobi**************
The optimal objective value is :22256
The solution time is :2.51661
The total number of chosen items is: 451
```


## Reference

<a id="ref1">[1]</a> Ryu, J., & Park, S. (2021). Robust solutions for stochastic and distributionally robust chance-constrained binary knapsack problems. *arXiv preprint arXiv:2105.11875*. ([ https://doi.org/10.48550/arXiv.2105.11875](https://doi.org/10.48550/arXiv.2105.11875)).  
