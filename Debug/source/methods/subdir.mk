################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../source/methods/SweepMethods.f95 

OBJS += \
./source/methods/SweepMethods.o 


# Each subdirectory must supply rules for building sources it contributes
source/methods/%.o: ../source/methods/%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

source/methods/SweepMethods.o: ../source/methods/SweepMethods.f95


