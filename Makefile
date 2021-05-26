MAIN?=PhaseOscilators
INCLUDE_ENSEMBLE?=Ensemble_param
INCLUDE_CONNECTIONS?=Ensemble_connections
INCLUDE_DINAMIC?=Ensemble_Dinamic

all: compile clean_o run

compile: ${MAIN} ${INCLUDE_ENSEMBLE} ${INCLUDE_CONNECTIONS} ${INCLUDE_DINAMIC}
	nvcc ${MAIN}.o ${INCLUDE_ENSEMBLE}.o ${INCLUDE_CONNECTIONS}.o ${INCLUDE_DINAMIC}.o -o ${MAIN}.exe -w

${INCLUDE_DINAMIC}: $(INCLUDE_DINAMIC).cu
	nvcc -c $< -o $@.o -w

${INCLUDE_CONNECTIONS}: $(INCLUDE_CONNECTIONS).cu
	nvcc -c $< -o $@.o -w

${INCLUDE_ENSEMBLE}: $(INCLUDE_ENSEMBLE).cu
	nvcc -c $< -o $@.o -w

${MAIN}: $(MAIN).cu
	nvcc -c $< -o $@.o -w

run:
	./$(MAIN).exe

clean_o:
	rm *.o

clean:
	rm ${MAIN}.exe