# Beliefs, Aggregate Risk, and the U.S. Housing Boom 
### By Margaret Jacobson, Federal Reserve Board

These are the views of the author and not those of the Federal Reserve Board, System or Staff.

The fortran code to solve the quantitative model is designed to run on a super computer on intel fortran with the mkl library. The same code can call two versions of the model

1) Standard Krussell Smith

2) Krusell Smith with learning

Model 1 can be called from `run_KS.sh` and uses `coefficients2.csv`

Modeld 2 can be called from `run_KS_learning.sh` and uses `coefficients3.csv`

These scripts run the makefile. Adding `$(debug)` to line 17 of  `makefile` will include debugger flags.

All versions can be called from various configurations of `parameters.f90`, which then calls in the requisite versions of coefficients from the various `coefficients*.xlsx` files



## Tables 1 and 3 (Table 2 is just parameters)
- `parameters.f90`: solve_for_coeffs=0, high_lamb=0, type_sim=1, no_MIT=0, no_pref=1, low_rate=0, high_rate=0
- `coefficients2.csv`: paste in column A from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. These are converged coefficients, to solve for coefficients as a fixed point set solve_for_coeffs=1.
- use `run_KS.sh`,
### output
`Text_tables.m`
 
## Figure 3b

###  Red dotted lines (credit conditions + income): 
- `parameters.f90`: solve_for_coeffs=0, **high_lamb=1**, type_sim=1, no_MIT=0, no_pref=1, low_rate=0, high_rate=0
- `coefficients2.csv`: paste in column C from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. These are converged coefficients, to solve for coefficients as a fixed point set solve_for_coeffs=1.
- use `run_KS.sh`

### Solid blue lines (beliefs + credit conditions + income): 
- `parameters.f90`: solve_for_coeffs=0, high_lamb=1, **type_sim=2**, no_MIT=0, no_pref=1, low_rate=0, high_rate=0, **gain=0.503**
- `coefficients2.csv`: paste in **column A** from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. 
- use `run_KS_learning.sh`
- Note: if you set gain=0, you should get the same results from `run_KS_learning.sh` as `run_KS.sh` provided that `coefficients2.csv` is the same.

### output
`Graphs_simulation.m`

## Figure 4
- `parameters.f90`: solve_for_coeffs=0, high_lamb=1, type_sim=2, no_MIT=0, no_pref=1, low_rate=0, high_rate=0, **gain=0.503, 0.59, 0.55, 0.4725, 0.4225**
- `coefficients2.csv`: paste in column A from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. 
- use `run_KS_learning.sh`

### output
 `Graphs_simulation_gain_sensitivity.m`

   
## Figures 5 to 7 and table 4
### Solid blue lines (same as Figure 3b)
- `parameters.f90`: solve_for_coeffs=0, high_lamb=1, type_sim=2, no_MIT=0, no_pref=1, low_rate=0, high_rate=0, gain=0.503
- `coefficients2.csv`: paste in column A from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. 
- use `run_KS_learning.sh`

### Dashed purple lines
- `parameters.f90`: solve_for_coeffs=0, high_lamb=1, **type_sim=1**, **no_MIT=1**, **no_pref=0**, low_rate=0, high_rate=0, gain=0.503
- **`coefficients3.csv`**: paste in column A from **`coefficients2_save.xlsx`** into column A of `coefficients2.csv`. 
- use `run_KS.sh'

### output
`Graphs_simulation.m`

## Figure 8

## Figure 9 

## Figure 10

## Figure 11

## Table 5 (Table 6 is updated by hand)
- `parameters.f90`: solve_for_coeffs=0, high_lamb=0, type_sim=1, no_MIT=0, no_pref=1, low_rate=0, high_rate=0
- `coefficients2.csv`: paste in column A from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. These are converged coefficients, to solve for coefficients as a fixed point set solve_for_coeffs=1.
- use `run_KS.sh`
- 
### output
`Text_tables.m`

## Tables 7 and 8
- `parameters.f90`: solve_for_coeffs=0, high_lamb=0, type_sim=1, no_MIT=0, no_pref=1, low_rate=0, high_rate=0
- `coefficients2.csv`: paste in column A from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. These are converged coefficients, to solve for coefficients as a fixed point set solve_for_coeffs=1.
- use `run_KS.sh`

### output 
`Tables_segmentation.m`

## Table 9 (To-do)

## Tables 10 and 11
- `parameters.f90`: solve_for_coeffs=0, **high_lamb=1**, type_sim=1, **no_MIT=1**, no_pref=1, low_rate=0, high_rate=0
- `coefficients2.csv`: paste in **column D** from `coefficients2_save.xlsx` into column A of `coefficients2.csv`. These are converged coefficients, to solve for coefficients as a fixed point set solve_for_coeffs=1.
- use `run_KS.sh`

### output
`Text_Tables_Appendix_Markov.m`

## Figure 16

## Figure 17

## Figure 18

## Figure 19

## Figure 20

## Figure 21

## Figure 22

## Figure 23



