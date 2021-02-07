################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include test/testclasses/subdir.mk
-include test/subdir.mk
-include source/solvers/writers/subdir.mk
-include source/solvers/subdir.mk
-include source/methods/subdir.mk
-include source/logger/subdir.mk
-include source/listener/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: NumbersMethods

# Tool invocations
NumbersMethods: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Fortran Linker'
	gfortran -fopenmp -o "NumbersMethods" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) NumbersMethods
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
