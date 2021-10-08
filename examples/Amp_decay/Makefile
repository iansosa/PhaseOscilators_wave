MAIN?=PhaseOscilators
INCLUDE_ENSEMBLE_PATH?=../../Ensemble_param
INCLUDE_ENSEMBLE?=Ensemble_param
INCLUDE_CONNECTIONS_PATH?=../../Ensemble_connections
INCLUDE_CONNECTIONS?=Ensemble_connections
INCLUDE_DINAMIC_PATH?=../../Ensemble_Dinamic
INCLUDE_DINAMIC?=Ensemble_Dinamic
INCLUDE_EVOLVE_PATH?=../../Ensemble_evolve
INCLUDE_EVOLVE?=Ensemble_evolve

all: compile clean_o run

compile: ${MAIN} ${INCLUDE_ENSEMBLE} ${INCLUDE_CONNECTIONS} ${INCLUDE_DINAMIC} ${INCLUDE_EVOLVE}
	nvcc ${MAIN}.o ${INCLUDE_ENSEMBLE}.o ${INCLUDE_CONNECTIONS}.o ${INCLUDE_DINAMIC}.o ${INCLUDE_EVOLVE}.o -o ${MAIN}.exe -w -Xcompiler -fopenmp 

${INCLUDE_EVOLVE}: $(INCLUDE_EVOLVE_PATH).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${INCLUDE_DINAMIC}: $(INCLUDE_DINAMIC_PATH).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${INCLUDE_CONNECTIONS}: $(INCLUDE_CONNECTIONS_PATH).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${INCLUDE_ENSEMBLE}: $(INCLUDE_ENSEMBLE_PATH).cu
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

${MAIN}: $(MAIN).cu 
	nvcc -c $< -o $@.o -w -Xcompiler -fopenmp

run:
	./$(MAIN).exe

clean_o:
	rm *.o

clean:
	rm ${MAIN}.exe