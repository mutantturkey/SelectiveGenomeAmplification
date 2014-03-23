// find string in 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

char **load_mers_from_file(FILE *fh, ssize_t *len) {
	char line[4096];

	size_t realloc_size = 1;
	*len = 0;

	char **mers = NULL;

	while((fgets(line, 4096, fh)) != NULL) { 
		size_t line_len = strlen(line);
		if(line_len == 0) 
			continue;

		line[line_len - 1] = '\0';

		mers = realloc(mers, sizeof(char *) * realloc_size);
		if(mers == NULL) {
			fprintf(stderr, "could not realloc mers\n");
			exit(EXIT_FAILURE);
		}

		char *cpy = malloc(line_len + 1);
		if(cpy == NULL) {
			fprintf(stderr, "could not alloc mers\n");
			exit(EXIT_FAILURE);
		}
		strncpy(cpy, line, line_len);
		mers[*len] = cpy;
		*len += 1;

		realloc_size++;
	}

	if(*len != 0)
		return mers;

	*len = -1;
	return NULL;
}

int main(int argc, char **argv){

	char buffer[BUFSIZ] = { 0 };
	char *buf, *start;
	ssize_t len = 0;
	ssize_t mer_len = 0;

	int save_size = 0;
	int cpy = 0;

	unsigned long long pos = 0;
	unsigned long long cpy_size = 0; 

	int i = 0;

	if(argc != 3)  {
		fprintf(stderr, "usage: strstream merlist.txt\n");
		exit(EXIT_FAILURE);
	}
	//
	FILE *infh = fopen(argv[2], "r");
	if(infh == NULL) {
		fprintf(stderr, "could not open %s\n", argv[2]);
		exit(EXIT_FAILURE);
	}
	// load mers
	FILE *fh = fopen(argv[1], "r");
	if(fh == NULL) {
		fprintf(stderr, "could not open %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	char **mers = load_mers_from_file(fh, &mer_len);
	if(mers == NULL) {
		fprintf(stderr, "could not load mers from %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}

	// get max argument length
	for(i = 0; i < mer_len; i++) {
		int current_len = strlen(mers[i]);
		if( current_len > save_size)
			save_size = current_len;
	}

	cpy = save_size - 1;
	cpy_size = BUFSIZ - cpy;

	buf = buffer;
	start = buf + cpy;

	// read into "start" (buf + cpy) from stdin
	while((len = fread(start, 1, cpy_size, infh)) != 0) {
		for(i = 0; i < mer_len; i++) {
			char *p = buffer;
				while((p = strstr(p, mers[i])) != NULL) {
				printf("%d %llu\n", i, pos + (p - buffer));
				p++;
			}
		}
		memcpy(buffer, buffer + len, cpy);
		pos = pos + len;
	}
	return 0;
}
