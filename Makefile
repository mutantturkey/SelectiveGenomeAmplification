VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra
SRC = src
BIN = bin
DEST = /usr/local/bin/

all: output_dir $(BIN)/strstream $(BIN)/melting_range $(BIN)/strstreamone

output_dir: 
	mkdir -p $(BIN)
$(BIN)/strstream: $(SRC)/strstream.c
	$(CC) $(SRC)/strstream.c -o $(BIN)/strstream $(CLIBS) $(CFLAGS)
$(BIN)/strstreamone: $(SRC)/strstreamone.c
	$(CC) $(SRC)/strstreamone.c -o $(BIN)/strstreamone $(CLIBS) $(CFLAGS)
$(BIN)/melting_range: $(SRC)/melting_range.c
	$(CC) $(SRC)/melting_range.c -o $(BIN)/melting_range $(CLIBS) $(CFLAGS)

clean:
	rm -vf $(BIN)/strstream $(BIN)/melting_range $(BIN)/strstreamone

install: all
	install -c $(BIN)/strstream $(BIN)/melting_range $(BIN)/strstreamone SelectiveGenomeAmplification select_mers.py score_mers.py $(DEST)

