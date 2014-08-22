VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra
DEST = /usr/local/bin/

all: output_dir bin/strstream bin/filter_melting_range bin/strstreamone bin/sequence_end_points

output_dir: 
	mkdir -p bin

bin/strstream: src/strstream.c
	$(CC) src/strstream.c -o bin/strstream $(CLIBS) $(CFLAGS)
bin/strstreamone: src/strstreamone.c
	$(CC) src/strstreamone.c -o bin/strstreamone $(CLIBS) $(CFLAGS)
bin/sequence_end_points: src/sequence_end_points.c 
	$(CC) src/sequence_end_points.c -o bin/sequence_end_points $(CLIBS) $(CFLAGS)
bin/filter_melting_range: src/filter_melting_range.c
	$(CC) src/filter_melting_range.c -o bin/filter_melting_range $(CLIBS) $(CFLAGS)

clean:
	rm -vf bin/* -Rv

install: all
	# c tools
	install -c bin/strstream $(DEST)
	install -c bin/filter_melting_range $(DEST)
	install -c bin/strstreamone $(DEST)
	install -c bin/sequence_end_points $(DEST)
	# bash scripts
	install -c SelectiveWholeGenomeAmplification $(DEST)
	install -c SelectiveWholeGenomeAmplificationUI $(DEST)
	install -c src/lock $(DEST)
	# python scripts
	install -c src/select_mers.py $(DEST)
	install -c src/score_mers.py $(DEST)
	install -c src/score_wrapper.sh $(DEST)
	install -c src/filter_melting_temperature.py $(DEST)
	install -c src/filter_max_consecutive_binding.py $(DEST)
	install -c src/filter_max_bg_mers.py $(DEST)
	install -c src/filter_average_binding.py $(DEST)
	install -c src/remove_mers.py $(DEST)
	install -c src/remove_mers_from_file.py $(DEST)
