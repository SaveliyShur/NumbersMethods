################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../test/testclasses/Logger_Test.f95 

F90_SRCS += \
../test/testclasses/SweepMethods_Test.f90 

OBJS += \
./test/testclasses/Logger_Test.o \
./test/testclasses/SweepMethods_Test.o 


# Each subdirectory must supply rules for building sources it contributes
test/testclasses/%.o: ../test/testclasses/%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

test/testclasses/Logger_Test.o: ../test/testclasses/Logger_Test.f95 source/logger/Logger_module.o

test/testclasses/%.o: ../test/testclasses/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

test/testclasses/SweepMethods_Test.o: ../test/testclasses/SweepMethods_Test.f90 source/methods/SweepMethods.o


