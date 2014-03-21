// Copyright 2013 Calvin Morrison
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
int main() {

	size_t len = 0;

	char buffer[4096];
	bool header = false;

	len = fread(&buffer, 1, 1, stdin);

	unsigned long long seq_length = 0;
	if(!errno) {
		if(buffer[0] == '>') {
			header = true;

			while((len = fread(&buffer, 1, 4096, stdin)) != 0) {
				size_t i = 0;                                                                    
				for(i = 0; i < len; i++) {
					if(buffer[i] == '>') {
						printf("%llu\n", seq_length);
						header = true;
						continue;
					}   
					else if(buffer[i] == '\n' && header == true) {
						header = false;
						continue;
					}
					if(header == false && buffer[i] != '\n') { 
						seq_length++;
					}   
				}   
			}   
		} 
		else { 
			fprintf(stderr, "this does not look like a fasta file\n"); 
			return EXIT_FAILURE;
		}
	} 
	else { 
		fprintf(stderr, "could not read file\n"); 
		return EXIT_FAILURE;
	}

	printf("%llu\n", seq_length);

	return EXIT_SUCCESS;
}

