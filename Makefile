VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra

all: strstream

strstream: strstream.c
	$(CC) strstream.c -o strstream $(CLIBS) $(CFLAGS)

clean:
	rm -vf kmer_total_count kmer_counts_per_sequence libkmer.so libkmer.a libkmer.o
