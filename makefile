#Setting compiler variables
compiler = ifort
# -save 
debug=-w -g -Og -check all -check bounds -warn all -debug all 
# -fp-stack-check  -ffpe-trap  -w-fpe0 -save -check all -check bounds -check noarg_temp_created  -debug all -fdump-tree-original -diag-enable warn -fp-model  -fcheck=all -fbounds-check -traceback -ftrap=common  -fp-stack-check -Wall -Wextra -fno-builtin
# -w-fpe0 -save -check all -check bounds -check noarg_temp_created  -debug all -fdump-tree-original -diag-enable warn -fp-model  -fcheck=all -fbounds-check -traceback -ftrap=common  -fp-stack-check -Wall -Wextra -fno-builtin
# --leak-check=full --show-leak-kinds=all 
#  -g -Og 
#  -pedantic-errors -ffpe-trap=zero,invalid,overflow,underflo -init=snan,arrays
#-check all -check bounds -check noarg_temp_created  -fbounds-check
#O1 is 3x faster in VFI than Ofast or O3 (7.3 minutes vs. 20.5 minutes)
#O1 small simulation is 3.72 minutes vs. 4.2 minutes with Ofast
#O1 large simulation is X minutes
#-fp-model=precise


options= -O1  -no-ipo  -qopenmp -fp-model=source -parallel -mkl -mtune=native  -march=native -shared-intel -mcmodel=large -traceback -no-wrap-margin 
#-fno-ipa-cp
#  -qoverride-limits 
#-heap-arrays
objects=parameters.o  mod_globals.o mod_interp.o mod_functions.o golden.o mod_matlab.o mod_rouwenhorst.o mod_discreteAR1.o orderpack.o mod_grids_fixed.o mod_plot.o mod_root1dim.o mod_initial_Y.o mod_initial_dist_full.o mod_initial_dist_full2.o mod_coeffs.o mod_bellman_bequest.o mod_continuation_value.o mod_bellman_new.o mod_EEE.o mod_VFI.o mod_function_interp.o mod_simulate_steady.o mod_simulate_bor_full.o mod_KS_regression2.o  mod_simulate_learning.o mod_csv.o main_full3.o

#Fixing m2c not found error (not sure why this works)
%.o : %.mod

#Making the executable by linking the object files
main.out : $(objects)
	$(compiler) -o main.o $(options) $(objects) 

#Making parameters.o
parameters.o: parameters.f90
	$(compiler) $(options) -c  parameters.f90

#Making mod_globals.o
mod_globals.o: mod_globals.f90 parameters.f90
	$(compiler) $(options) -c  mod_globals.f90

#Making mod_interp.o
mod_interp.o : mod_interp.f90 
	$(compiler) $(options) -c mod_interp.f90 
	
#Making mod_functions.o
mod_functions.o : mod_functions.f90 parameters.f90 mod_globals.f90 mod_interp.f90
	$(compiler) $(options) -c mod_functions.f90 
	
#Making golden.o
golden.o : golden.f90 
	$(compiler) $(options) -c golden.f90 
	
#Making mod_matlab.o 
mod_matlab.o: mod_matlab.f90
	$(compiler) $(options) -c mod_matlab.f90 

#Making mod_rouwenhorst.o
mod_rouwenhorst.o: mod_rouwenhorst.f90 parameters.f90 mod_matlab.f90 
	$(compiler) $(options) -c  mod_rouwenhorst.f90

#Making mod_discreteAR1.f90
mod_discreteAR1.o: mod_discreteAR1.f90 parameters.f90 mod_matlab.f90 
	$(compiler) $(options) -c mod_discreteAR1.f90

#Making orderpack.o
orderpack.o: orderpack.f90 
	$(compiler) $(options) -c orderpack.f90

#Making mod_grids_fixed.o
mod_grids_fixed.o: mod_grids_fixed.f90 parameters.f90 mod_functions.f90 mod_globals.f90 mod_matlab.f90 mod_rouwenhorst.f90 mod_discreteAR1.f90 orderpack.f90
	$(compiler) $(options) -c mod_grids_fixed.f90

#Making mod_plot.o
mod_plot.o: mod_plot.f90 
	$(compiler) $(options) -c mod_plot.f90

#Making mod_root1dim.o
mod_root1dim.o: mod_root1dim.f90 mod_plot.f90
	$(compiler) $(options) -c mod_root1dim.f90

#Making mod_initial_Y.o
mod_initial_Y.o: mod_initial_Y.f90 parameters.f90 mod_matlab.f90 mod_interp.f90
	$(compiler) $(options) -c mod_initial_Y.f90

#Making mod_initial_dist_full.o
mod_initial_dist_full.o: mod_initial_dist_full.f90 parameters.f90 mod_matlab.f90 mod_interp.f90
	$(compiler) $(options) -c mod_initial_dist_full.f90

#Making mod_initial_dist_full2.o
mod_initial_dist_full2.o: mod_initial_dist_full2.f90 parameters.f90 mod_matlab.f90 mod_interp.f90
	$(compiler) $(options) -c mod_initial_dist_full2.f90

#Making mod_coeffs.o
mod_coeffs.o: mod_coeffs.f90 parameters.f90 mod_matlab.f90 mod_discreteAR1.f90
	$(compiler) $(options) -c mod_coeffs.f90

#Making mod_bellman_bequest.o
mod_bellman_bequest.o: parameters.f90 mod_matlab.f90 mod_globals.f90 mod_functions.f90
	$(compiler) $(options) -c mod_bellman_bequest.f90

#Making mod_continuation_value.o
mod_continuation_value.o: mod_continuation_value.f90 parameters.f90 mod_matlab.f90 mod_globals.f90 mod_interp.f90 mod_functions.f90
	$(compiler) $(options) -c mod_continuation_value.f90

#Making mod_bellman_new.o
mod_bellman_new.o: mod_bellman_new.f90 parameters.f90 mod_matlab.f90 mod_globals.f90 mod_functions.f90 mod_interp.f90
	$(compiler) $(options) -c mod_bellman_new.f90

#Making mod_EEE.o
mod_EEE.o: mod_EEE.f90 parameters.f90 mod_interp.f90 mod_matlab.f90
	$(compiler) $(options) -c mod_EEE.f90

#Making mod_VFI.o
mod_VFI.o: mod_VFI.f90 parameters.f90 mod_globals.f90 mod_matlab.f90 mod_grids.f90 mod_interp.f90 mod_bellman_bequest.f90 mod_continuation_value.f90 mod_bellman_new.f90 mod_EEE.f90
	$(compiler) $(options) -c mod_VFI.f90
	
#Making mod_funtion_interp.o
mod_function_interp.o: parameters.f90 mod_functions.f90 mod_matlab.f90 mod_interp.f90 
	$(compiler) $(options) -c mod_function_interp.f90

#Making mod_simulate_steady.o
mod_simulate_steady.o: parameters.f90 mod_functions.f90 mod_matlab.f90 mod_interp.f90 mod_root1dim.f90
	$(compiler) $(options) -c mod_simulate_steady.f90

#Making mod_simulate_bor_full.o
mod_simulate_bor_full.o: mod_simulate_bor_full.f90 parameters.f90 mod_functions.f90 mod_matlab.f90 mod_interp.f90 mod_root1dim.f90
	$(compiler) $(options) -c mod_simulate_bor_full.f90

#Making mod_KS_regression.o
mod_KS_regression2.o: mod_KS_regression2.f90 parameters.f90 mod_matlab.f90
	$(compiler) $(options) -c mod_KS_regression2.f90
	

#Making mod_simulate_learning.o
mod_simulate_learning.o: mod_simulate_learning.f90 parameters.f90 mod_functions.f90 mod_matlab.f90 mod_interp.f90 mod_root1dim.f90
	$(compiler) $(options) -c mod_simulate_learning.f90

#Making mod_csv.o
mod_csv.o: mod_csv.f90 
	$(compiler) $(options) -c mod_csv.f90

#Making main_full3.o
main_full3.o: main_full3.f90 parameters.f90 mod_globals.f90 mod_interp.f90 mod_functions.f90 golden.f90 mod_matlab.f90 mod_rouwenhorst.f90 mod_discreteAR1.f90 orderpack.f90 mod_grids_fixed.f90 mod_plot.f90 mod_root1dim.f90 mod_initial_Y.f90 mod_initial_dist_full.f90 mod_initial_dist_full2.f90 mod_coeffs.f90 mod_bellman_bequest.f90 mod_continuation_value.f90 mod_bellman_new.f90 mod_EEE.f90 mod_VFI.f90 mod_function_interp.f90 mod_simulate_steady.f90 mod_simulate_bor_full.f90 mod_KS_regression2.f90 mod_simulate_learning.f90  mod_csv.f90
	$(compiler) $(options) -c main_full3.f90

#Making Clean Command
clean:
	rm -f *.o *.mod *.out
