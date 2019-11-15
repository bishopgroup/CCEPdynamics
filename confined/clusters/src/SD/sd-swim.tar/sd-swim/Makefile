# The location of the compiler executable
COMP =  /opt/intel/Compiler/11.0/083/bin/intel64/./ifort

# Any optimation or debugging options for the compiler
FLAGS = -g -traceback
#FLAGS = -O3 -msse4.1 -parallel -i-static

FLINK = -Vaxlib

sd.x: glob.o blas.o chol.o init.o lub.o dist.o mob.o vel.o traj.o out.o sd.o 
	$(COMP) $(FLAGS) glob.o blas.o chol.o init.o lub.o dist.o mob.o vel.o traj.o out.o sd.o -o sd.x $(FLINK)

glob.o: glob.f90
	$(COMP) $(FLAGS) glob.f90 -c

blas.o: blas.f
	$(COMP) $(FLAGS) blas.f -c

chol.o: chol.f90
	$(COMP) $(FLAGS) chol.f90 -c

init.o: init.f90
	$(COMP) $(FLAGS) init.f90 -c

lub.o: lub.f90
	$(COMP) $(FLAGS) lub.f90 -c

dist.o: dist.f90
	$(COMP) $(FLAGS) dist.f90 -c

mob.o: mob.f90
	$(COMP) $(FLAGS) mob.f90 -c

vel.o: vel.f90
	$(COMP) $(FLAGS) vel.f90 -c

traj.o: traj.f90
	$(COMP) $(FLAGS) traj.f90 -c

out.o: out.f90
	$(COMP) $(FLAGS) out.f90 -c

sd.o: sd.f90
	$(COMP) $(FLAGS) sd.f90 -c


clean:
	rm *.mod
	rm *.o
