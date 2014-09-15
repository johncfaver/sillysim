
compiler=gfortran
flags=-std=f2003 -fopenmp  
debug=-g
warn=-Wall
opt=-Ofast 

target=main.exe
sourcefiles=main.f90 atom_class.f90 molecule_class.f90 file_handler_class.f90 \
            geometry.f90 potential_class.f90 simulator_class.f90

######################################
all: $(sourcefiles) $(obj); $(compiler) $(flags) $(opt) $(sourcefiles) -o $(target)

debug: $(sourcefiles) $(obj); $(compiler) $(flags) $(opt) $(debug) $(warn) $(sourcefiles) -o $(target)

obj: $(sourcefiles); $(compiler) $(flags) $(opt) -c $(sourcefiles)

clean: ; rm *.mod *.o
