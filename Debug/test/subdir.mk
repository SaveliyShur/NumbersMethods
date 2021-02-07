################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../test/TestAllModules.f90 

OBJS += \
./test/TestAllModules.o 


# Each subdirectory must supply rules for building sources it contributes
test/%.o: ../test/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

test/TestAllModules.o: ../test/TestAllModules.f90 test/testclasses/Logger_Test.o test/testclasses/SweepMethods_Test.o


