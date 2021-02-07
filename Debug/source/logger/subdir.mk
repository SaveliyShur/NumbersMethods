################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../source/logger/Logger_module.f95 

OBJS += \
./source/logger/Logger_module.o 


# Each subdirectory must supply rules for building sources it contributes
source/logger/%.o: ../source/logger/%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

source/logger/Logger_module.o: ../source/logger/Logger_module.f95


