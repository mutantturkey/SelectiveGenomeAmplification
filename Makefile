VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra

all: strstream melting_range strstreamone

strstream: strstream.c
	$(CC) strstream.c -o strstream $(CLIBS) $(CFLAGS)
strstreamone: strstreamone.c
	$(CC) strstreamone.c -o strstreamone $(CLIBS) $(CFLAGS)
melting_range: melting_range.c
	$(CC) melting_range.c -o melting_range $(CLIBS) $(CFLAGS)

clean:
	rm -vf strstream melting_range 
