MAIN?=PhaseOscilators
INCLUDE_ENSEMBLE?=Ensemble_param

all: compile clean_o run

compile: ${MAIN} ${INCLUDE_ENSEMBLE}
	nvcc ${MAIN}.o ${INCLUDE_ENSEMBLE}.o -o ${MAIN}.exe -w 

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