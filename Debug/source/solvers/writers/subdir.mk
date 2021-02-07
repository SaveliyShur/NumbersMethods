################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../source/solvers/writers/Writer.f95 

OBJS += \
./source/solvers/writers/Writer.o 


# Each subdirectory must supply rules for building sources it contributes
source/solvers/writers/%.o: ../source/solvers/writers/%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

source/solvers/writers/Writer.o: ../source/solvers/writers/Writer.f95 source/logger/Logger_module.o


