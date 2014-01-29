#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

float melting_temperature(char *mer) {

	float a = 0;
	float c = 0;
	float g = 0;
	float t = 0;
	size_t i = 0;

	for(i = 0; i < strlen(mer); i++) {
		switch(mer[i]) {
			case 'A':
			a++;
			break;
			case 'C':
			c++;
			break;
			case 'G':
			g++;
			break;
			case 'T':
			t++;
			break;
			default:
			break;
		}
	}

	if(strlen(mer) < 13)
		return  ((a+t) * 2) + ((c+g) * 4);
	else
		return 64.9 + 41.0*(g+c-16.4)/(a+t+g+c);
}

int main(int argc, char **argv){

	if(argc < 3) {
		printf("please supply the min and max as stdargs");
		exit(EXIT_FAILURE);
	}
	float min = atof(argv[1]);
	float max = atof(argv[2]);
	
	char mer[64];
	int count = 0;

	while(fscanf(stdin, "%s\t%d\n", mer, &count) == 2) {
		float temp = melting_temperature(mer);
		if( (temp >= min) && (temp <= max) )
			printf("%s\t%d\n", mer, count);
	}

	exit(EXIT_SUCCESS);
}
