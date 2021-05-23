APP?=PhaseOscilators

all: ${APP} run

${APP}: $(APP).cu
	nvcc -O2 -o $@.exe $< -w -Xcompiler -fopenmp 

run: $(APP)	
	./$(APP).exe