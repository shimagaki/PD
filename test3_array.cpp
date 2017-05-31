#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#define L 10
#define pi M_PI/180.0
#define rand01 (double)rand()/RAND_MAX
// The size of the system is L
// The interaction is next nearest interaction.
// t[x][y] is an angle, theta_xy, where xy is location on the lattice.

/********
 *=> array should be vector.
 *=> the size of matrix should be changeable.
 * ************/
template <class T>
void print_mat(T mat[L][L]){
	int i,j;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
		printf("%2.1f ",mat[i][j]);
		}
		printf("\n");
	}
}

/**/
double  erg(double t_s,double t_sr,double t_sl,double t_su,double t_sd){
	return 4.0 - (cos(t_s-t_sr)+cos(t_s-t_sl)
			+cos(t_s-t_su)+cos(t_s-t_sd));
}

double seaqHB(double b,double t_s,double t_sr,double t_sl,double t_su,double t_sd){
	double u, dt, dE_prev,dE_prop,r;
	u = rand01;
	//dt is proposal of the theta.	
	dt = 2.0*M_PI*rand01/180.0;	
	dE_prev = erg(t_s,t_sr,t_sl,t_su,t_sd);	
	dE_prop = erg(t_s+dt,t_sr,t_sl,t_su,t_sd);	
	r = exp(-b*(dE_prop-dE_prev));
	r = r/(1.0+r); 
	if(u<r){
	t_s = dt;
	}
	return t_s;
}

/*sequential update*/
double mc(double b,double t[L][L]){
	int i,j,t_ij;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
		t_ij=seaqHB(b,t[i][j],t[(i-1+L)%L][j],t[(i+1)%L][j]
				,t[i][(j-1+L)%L],t[i][(j+1)%L]);	
		t[i][j]=t_ij;	
		}
	}
}

double *tempfunc(int a, double* b){
	for(int i=0; i<a ; ++i){
		b[i]=i;
	}
	return b;	
}

double *tempfunc2(double *b[L], int a){
	for(int i=0; i<a ; ++i){
		for(int j=0; j<a ; ++j){
			b[i][j]  = i*j;
		}
	}
	return b[1];	
}

//using Box-Muller method
 double rand_normal(double mu, double sigma){
	double z = sqrt( -2.0 * log(rand01) )* sin( 2.0 * M_PI * rand01);
	return mu + sigma * z ; 
}	

int main(int argc, char* argv[]){
	double t[L],M,E;
	double* t2[L+1];
	int i,j;
	srand(0);
	/*Iitialize of the t.*/		
	double* m_temp;
	m_temp = tempfunc(L,t);
	//std::cout<<sizeof(tempfunc(L,t))<<std::endl;
//	std::cout<<m_temp[9]<<std::endl;	
	std::cout<<"test output" <<std::endl;
	std::cout<<tempfunc2(t2,L)<<std::endl;	
	return 0;
}
