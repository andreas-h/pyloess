CFLAGS = -cckr -O
FFLAGS = -O
OBJ = loessc.o loess.o predict.o misc.o loessf.o

gas: gas.x
	gas.x
gas.x: gas.o $(OBJ)
	cc -o gas.x gas.o $(OBJ) -llinpack -lblas -lm -lF77

madeup: madeup.x
	madeup.x
madeup.x: madeup.o $(OBJ)
	cc -o madeup.x madeup.o $(OBJ) -llinpack -lblas -lm -lF77

ethanol: ethanol.x
	ethanol.x
ethanol.x: ethanol.o $(OBJ)
	cc -o ethanol.x ethanol.o $(OBJ) -llinpack -lblas -lm -lF77

air: air.x
	air.x
air.x: air.o $(OBJ)
	cc -o air.x air.o $(OBJ) -llinpack -lblas -lm -lF77

galaxy: galaxy.x
	galaxy.x
galaxy.x: galaxy.o $(OBJ)
	cc -o galaxy.x galaxy.o $(OBJ) -llinpack -lblas -lm -lF77

clean:
	rm -f *.o *.x core
