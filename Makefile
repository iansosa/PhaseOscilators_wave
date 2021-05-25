MAIN?=PhaseOscilators
INCLUDE_ENSEMBLE?=Ensemble_param
INCLUDE_CONNECTIONS?=Ensemble_connections

all: compile clean_o run

compile: ${MAIN} ${INCLUDE_ENSEMBLE} ${INCLUDE_CONNECTIONS}
	nvcc ${MAIN}.o ${INCLUDE_ENSEMBLE}.o ${INCLUDE_CONNECTIONS}.o -o ${MAIN}.exe -w

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