#!/bin/make

#Compilerbefehl im Terminal (Compilerexecutable)
F = gfortran 
#Compilerflags
CFLAGS = -O -Wall #-pedantic -w -fbounds-check -g#-fopenmp
#eventueller Pfad
S =./src_new/

OBJ = $(S)global.o $(S)readeph2.o $(S)getephem.o  $(S)input.o $(S)transform.o $(S)output.o $(S)closeenc.o $(S)radau.o $(S)symp.o  $(S)main.o


all: $(OBJ)
	mv *.mod $(S)
#	rm *.mod
	$(F) $(CFLAGS) $(OBJ) -o nine.exe


$(S)global.o: $(S)global_m.f90 
	$(F) $(CFLAGS) -c $(S)global_m.f90 -o $(S)global.o

$(S)readeph2.o: $(S)readeph2.f
	$(F) $(CFLAGS) -c $(S)readeph2.f -o $(S)readeph2.o

$(S)getephem.o: $(S)global_m.f90  $(S)getephem_m.f90 $(S)readeph2.f 
	$(F) $(CFLAGS) -c $(S)getephem_m.f90 -o $(S)getephem.o

$(S)input.o:$(S)input_m.f90 $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90
	$(F) $(CFLAGS) -c $(S)input_m.f90 -o $(S)input.o

$(S)transform.o:$(S)transform_m.f90 $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90
	$(F) $(CFLAGS) -c $(S)transform_m.f90 -o $(S)transform.o

$(S)output.o:$(S)output_m.f90 $(S)transform_m.f90 $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90
	$(F) $(CFLAGS) -c $(S)output_m.f90 -o $(S)output.o

$(S)closeenc.o:$(S)closeenc_m.f90 $(S)output_m.f90 $(S)transform_m.f90 $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90
	$(F) $(CFLAGS) -c $(S)closeenc_m.f90 -o $(S)closeenc.o

	
$(S)radau.o: $(S)radau_m.f90 $(S)output_m.f90 $(S)transform_m.f90 $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90 $(S)closeenc_m.f90 
	$(F) $(CFLAGS) -c $(S)radau_m.f90 -o $(S)radau.o

$(S)symp.o: $(S)symp_m.f90 $(S)output_m.f90 $(S)transform_m.f90 $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90 $(S)closeenc_m.f90 
	$(F) $(CFLAGS) -c $(S)symp_m.f90 -o $(S)symp.o


$(S)main.o:$(S)main.f90  $(S)output_m.f90  $(S)transform_m.f90  $(S)input_m.f90  $(S)global_m.f90 $(S)readeph2.f $(S)getephem_m.f90 $(S)radau_m.f90 $(S)symp_m.f90 $(S)closeenc_m.f90 
	$(F) $(CFLAGS) -c $(S)main.f90 -o $(S)main.o


clean:	
	rm $(S)*.mod  $(S)*.o 
