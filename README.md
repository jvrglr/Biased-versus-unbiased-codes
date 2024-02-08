# Biased versus unbiased numerical methods for stochastic simulations: application to contagion processes

Python and FORTRAN codes used to generate plots in J. Aguilar, J.J. Ramasco, and R. Toral  "Biased versus unbiased numerical methods for stochastic simulations: application to contagion processes".

This collection of code is provided under an open-source license for free use and modification. If you find these codes useful in your research or project, we kindly request that you acknowledge the source by citing the following article:


Aguilar, J., Ramasco, J. J., & Toral, R.]. "Numerical methods for stochastic simulations: application to contagion processes. ."  2003, arXiv:2305.02902.

If you have any questions or improvements, feel free to reach out (jvrglrschz@gmail.com).

## Code list FORTRAN
FORTRAN 2008 is used in the examples dealing with meta-population systems (Section VI.B in our article). 

Compilation is achieved using the gfortran compiler:

```
gfortran  -o exe.x Random2.f dranxor.f90 ignbin.f  M_declarations.f08 M_functions.f08 M_subroutines.f08 main_compare_B_U.f08 
```

Once compiled, you can execute the program with:
```
time ./exe.x
```

1. **main_compare_B_U.f08** Calls main functions and subroutines to do the computation.


