################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../Main.f95 

OBJS += \
./Main.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

Main.o: ../Main.f95 source/logger/Logger_module.o source/solvers/MethodOfEstablishing_Plate.o source/solvers/NavieStocksGas.o source/solvers/PrandtlSolverGas.o source/solvers/PrandtlSolver_Plate.o test/TestAllModules.o


