This project centers around the Complete Neglect of Differential Overlap/2 (CNDO/2) method, a semi-empirical quantum mechanical used to estimate molecular electronic structure and potential energy surface. 

This program features a large, complex, and organized class to fascillate the subroutines associated with the CNDO/2 algorithm. In order to improve the overall performance, the Fock Matrix and Core Hamiltonian calculations were parallelized using OpenMP. This allows my program to handle relatively large molecules.

# Build the code
```
make
```

# Run the code 
```
./main <filename>
```

To remove the executables and object files from the directory, run
```
make clean
```

Implementation of the CNDO/2 procedure can be found in `source/CNDO2.cpp` (below the overlap integral calculation functions) and `lib/CNDO2.h`.