compiler = mpiifort

flagsGNU = -Ofast -march=native -mtune=native -pipe -fopenmp
flagsIntel = -xHOST -O3 -no-prec-div -fp-model fast=2 -mkl -parallel -align -inline-level=2 -march=native -mtune=native -mcpu=native -qopenmp -mcmodel=medium -shared-intel

# -autodouble -shared-intel -mcmodel=medium 

libsGNU = -llapack -lblas
libsIntel =  

flags = $(flagsIntel)
libs = $(libsIntel)

.PHONY: all clean

all:	src/jmp.f90 jmpMod.o
	$(compiler) $(flags) -o bin/jmp *.o src/jmp.f90 ${libs}

MPIMod.o:	src/MPIMod.f90
	$(compiler) $(flags) -c src/MPIMod.f90

sim.o:	src/sim.f90
	$(compiler) $(flags) -c src/sim.f90

Param.o:	src/Param.f90
	$(compiler) $(flags) -c src/Param.f90

NL.o:	src/NL.f90 Param.o
	$(compiler) $(flags) -c src/NL.f90

jmpMod.o:	src/jmpMod.f90 NL.o MPIMod.o sim.o Param.o
	$(compiler) $(flags) -c src/jmpMod.f90

clean:
	rm -f bin/jmp ./*.mod ./*.o

cleanresults:
	rm -f results/*.txt

