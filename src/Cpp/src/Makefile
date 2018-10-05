CC = g++

LIBS  = gsl gslcblas m 
CCFLAGS += -g3 -Wall -std=c++11 -fmax-errors=5 
CCFLAGS += -O2 -ffast-math -ftree-vectorize 
CCFLAGS += -fmessage-length=0 -MMD -MP #-Werror 


OBJS = LISA.o MorletGabor.o TDI.o Wavelet.o G_B_MCMC.o IO.o 

all : $(OBJS) g_b

LISA.o : LISA.cpp LISA.h Constants.h
	$(CC) $(CCFLAGS) -c LISA.cpp
	
MorletGabor.o : MorletGabor.cpp MorletGabor.h
	$(CC) $(CCFLAGS) -c MorletGabor.cpp
	
TDI.o : TDI.cpp TDI.h LISA.h Constants.h
	$(CC) $(CCFLAGS) -c TDI.cpp	
	
Wavelet.o : Wavelet.cpp Wavelet.h LISA.h Constants.h
	$(CC) $(CCFLAGS) -c Wavelet.cpp

G_B_MCMC.o : G_B_MCMC.cpp G_B_MCMC.h TDI.h Constants.h LISA.h Wavelet.h
	$(CC) $(CCFLAGS) -c G_B_MCMC.cpp

IO.o : IO.cpp IO.h LISA.h Wavelet.h G_B_MCMC.h
	$(CC) $(CCFLAGS) -c IO.cpp

g_b : $(OBJS) Glitch_Burst.cpp ProgressBar.hpp cxxopts.hpp LISA.h TDI.h IO.h Wavelet.h G_B_MCMC.h Constants.h MorletGabor.h
	$(CC) $(CCFLAGS) -o g_b Glitch_Burst.cpp $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

clean: 
	rm *.o g_b
