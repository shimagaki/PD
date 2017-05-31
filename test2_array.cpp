#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define L 10
#define pi M_PI/180.0
#define rand01 (double)rand()/RAND_MAX

int main(int argc, char* argv[]){
	double t[L][L],M,E;
	int i,j;
	srand(0);
	/*Iitialize of the t.*/	
	std::cout<<pi<<std::endl;
	std::cout<<M_PI<<std::endl;
	for(i=0;i<10;i++){
		std::cout<<rand01<<std::endl;
	}
	return 0;
}
