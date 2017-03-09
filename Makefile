###################################################################
#
#
##################################################################

CC=	gcc

CSRC=	./src/

CFLAGS=	-O3 -g -Wall `/usr/bin/xml2-config --cflags`
#CFLAGS=	-O3 -g -fPIC -Wall

INC=	-I./include

LFLAGS= -lnetcdf `/usr/bin/xml2-config --libs` -lm

COBJ=	$(CSRC)main.o \
	$(CSRC)jutil.o \
	$(CSRC)readXML.o \
	$(CSRC)vortex.o \
	$(CSRC)mrvsmooth.o


OBJ=	$(COBJ) 

EXEC=	./bin/synth

$(EXEC):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LFLAGS)

$(COBJ) : %.o : %.c
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

clean:
	rm $(COBJ)
	rm $(EXEC)
