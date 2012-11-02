################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BlasWrappers.cpp \
../Cholesky.cpp \
../MVEE.cpp \
../Utils.cpp \
../main.cpp \
../minVolume.cpp 

OBJS += \
./BlasWrappers.o \
./Cholesky.o \
./MVEE.o \
./Utils.o \
./main.o \
./minVolume.o 

CPP_DEPS += \
./BlasWrappers.d \
./Cholesky.d \
./MVEE.d \
./Utils.d \
./main.d \
./minVolume.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/spyros/gsl-1.15 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


