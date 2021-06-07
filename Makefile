MAIN?=PhaseOscilators
INCLUDE_ENSEMBLE?=Ensemble_param
INCLUDE_CONNECTIONS?=Ensemble_connections
INCLUDE_DINAMIC?=Ensemble_Dinamic
INCLUDE_EVOLVE?=Ensemble_evolve

all: compile clean_o run

compile: ${MAIN} ${INCLUDE_ENSEMBLE} ${INCLUDE_CONNECTIONS} ${INCLUDE_DINAMIC} ${INCLUDE_EVOLVE}
	nvcc ${MAIN}.o ${INCLUDE_ENSEMBLE}.o ${INCLUDE_CONNECTIONS}.o ${INCLUDE_DINAMIC}.o ${INCLUDE_EVOLVE}.o -o ${MAIN}.exe -w -Xcompiler -fopenmp 

${INCLUDE_EVOLVE}: $(INCLUDE_EVOLVE).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${INCLUDE_DINAMIC}: $(INCLUDE_DINAMIC).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${INCLUDE_CONNECTIONS}: $(INCLUDE_CONNECTIONS).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${INCLUDE_ENSEMBLE}: $(INCLUDE_ENSEMBLE).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${MAIN}: $(MAIN).cu 
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

run:
	./$(MAIN).exe

clean_o:
	rm *.o

clean:
	rm ${MAIN}.exe