################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../source/listener/Listener.f95 

OBJS += \
./source/listener/Listener.o 


# Each subdirectory must supply rules for building sources it contributes
source/listener/%.o: ../source/listener/%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

source/listener/Listener.o: ../source/listener/Listener.f95 source/logger/Logger_module.o


