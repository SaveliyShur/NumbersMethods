################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../source/solvers/MethodOfEstablishing_Plate.f95 \
../source/solvers/NavieStocksGas.f95 \
../source/solvers/PrandtlSolverGas.f95 \
../source/solvers/PrandtlSolver_Plate.f95 

OBJS += \
./source/solvers/MethodOfEstablishing_Plate.o \
./source/solvers/NavieStocksGas.o \
./source/solvers/PrandtlSolverGas.o \
./source/solvers/PrandtlSolver_Plate.o 


# Each subdirectory must supply rules for building sources it contributes
source/solvers/%.o: ../source/solvers/%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

source/solvers/MethodOfEstablishing_Plate.o: ../source/solvers/MethodOfEstablishing_Plate.f95 source/listener/Listener.o source/logger/Logger_module.o

source/solvers/NavieStocksGas.o: ../source/solvers/NavieStocksGas.f95 source/listener/Listener.o source/logger/Logger_module.o source/solvers/writers/Writer.o

source/solvers/PrandtlSolverGas.o: ../source/solvers/PrandtlSolverGas.f95 source/listener/Listener.o source/logger/Logger_module.o source/methods/SweepMethods.o

source/solvers/PrandtlSolver_Plate.o: ../source/solvers/PrandtlSolver_Plate.f95 source/listener/Listener.o source/logger/Logger_module.o source/methods/SweepMethods.o


