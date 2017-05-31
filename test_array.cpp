#include <iostream>
#include <stdio.h>
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
	dt = 2.0*M_PI*rand01/18.0;	
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
double (*mc(double b,double t[L][L]))[L]{
	int i,j,t_ij;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
		t_ij=seaqHB(b,t[i][j],t[(i-1+L)%L][j],t[(i+1)%L][j]
				,t[i][(j-1+L)%L],t[i][(j+1)%L]);	
		t[i][j]=t_ij;	
		}
	}
	return t;
}

//using Box-Muller method
 double rand_normal(double mu, double sigma){
	double z = sqrt( -2.0 * log(rand01) )* sin( 2.0 * M_PI * rand01);
	return mu + sigma * z ; 
}	

int main(int argc, char* argv[]){
	double b,t[L][L],M,E;
	int i,j;
	std::cout<<"temperature:"<<std::endl;
	std::cin>>b;	
	srand(0);
	/*Iitialize of the t.*/	
	for(i=0;i<L;i++){
		for(j=i;j<L;j++){
		t[i][j] = 2*M_PI*(rand01);
		t[j][i] = 2*M_PI*(rand01);
		}
	}
	print_mat(t);
	printf("\n\n");	
	/*How can I alocate the mc(*,*) ?*/
	double** m_temp;
	m_temp = mc(b,t);
	print_mat(mc(b,t));	
	return 0;
}
