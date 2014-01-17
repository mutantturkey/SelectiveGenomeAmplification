// find string in 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv){

	char buffer[BUFSIZ] = { 0 };
	char *buf, *start;
	ssize_t len = 0;

	int cpy = 0;

	unsigned long long pos = 0;
	unsigned long long cpy_size = 0; 

	// get max argument length
	int save_size = strlen(argv[1]);

	cpy = save_size - 1;
	cpy_size = BUFSIZ - cpy;

	buf = buffer;
	start = buf + cpy;

	// copy our first cpy length into the first part of our buffer
	len = fread(buffer, 1, cpy, stdin);
	if(len == 0) 
		exit(EXIT_FAILURE);

	// read into "start" (buf + cpy) from stdin
	while((len = fread(start, 1, cpy_size, stdin)) != 0) {
		char *p = buffer;
		while((p = strstr(p, argv[1])) != NULL) {
			printf("%llu\n", pos + (p - buffer));
			p++;
		}
		memcpy(buffer, buffer + len, cpy);
		pos = pos + len;
	}
	return 0;
}
