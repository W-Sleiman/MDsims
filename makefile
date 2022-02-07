NAME = vv

CC = icc
#CFLAGS = -g 
CFLAGS = -O3 -pg
LFLAGS = -lm

SOURCE = vv_md.c io.c integrate.c measure.c gasdev.c

OBJECT = vv_md.o io.o integrate.o measure.o gasdev.o


$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)

io.o integrate.o gasdev.o: defs.h

clean:
		rm -f $(OBJECT)
