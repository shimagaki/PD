#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#define L 30 
#define pi M_PI/180.0
#define rand01 (double)rand()/RAND_MAX
// The size of the system is L
// The interaction is next nearest interaction.
// t[x][y] is an angle, theta_xy, where xy is location on the lattice.

/********
 *=> array should be vector.
 *=> the size of matrix should be changeable.
 * ************/
/*
template <class T>
void print_mat(T mat[L][L]){
	int i,j;
	for(auto v1:v){
		for(auto v2:v1){
			std::cout<<v2<<"  ";	
		}
		std::cout<<std::endl;
	}
	return 0;
}
	
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
		printf("%2.1f ",mat[i][j]);
		}
		printf("\n");
	}
}
*/
template <class T>
void print_mat( std::vector<std::vector<T> > v){
	int i,j;
	for(i=0;i<v.size();i++){
		for(j=0;j<v[0].size();j++){
			std::cout<<v[i][j]<<"  ";
		}
		std::cout<<std::endl;
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

void mc(double b, std::vector<std::vector<double> >& t){
	int i=0,j=0;
	double t_ij;
	for(i=0;i<t.size();i++){
		for(j=0;j<t[0].size();j++){
			t_ij=seaqHB(b,t[i][j],t[(i-1+L)%L][j],t[(i+1)%L][j]
					,t[i][(j-1+L)%L],t[i][(j+1)%L]);	
			t[i][j]=t_ij;
		//	printf("t[%d,%d]=%2.3f",i,j,t_ij);	
		}
	}
}

//using Box-Muller method
 double rand_normal(double mu, double sigma){
	double z = sqrt( -2.0 * log(rand01) )* sin( 2.0 * M_PI * rand01);
	return mu + sigma * z ; 
}	

template <class T>
void write_mat(int n,double b, std::vector<std::vector<T> > v){
	std::string fname="xy-beta"+std::to_string(b)+"-tMC"+std::to_string(n)+".dat";
	std::ofstream file(fname);
	int i,j;
	for(i=0;i<v.size();i++){
		for(j=0;j<v[0].size();j++){
			file<<v[i][j]<<"  ";
		}
		file<<std::endl;
	}
	file.close();
}
int main(int argc, char* argv[]){
	double b,M,E;
	std::vector<std::vector<double> > t(L,std::vector<double>(L,0.0)); 
	std::vector<std::vector<double> > t_temp(L,std::vector<double>(L,0.0)); 
	int i,j,n,t_mc;
	std::cout<<"1/T="<<std::endl;
	std::cin>>b;	
	srand(0);
	std::cout<<"t_mc="<<std::endl;
	std::cin>>t_mc;	
	print_mat(t);
	for(n=0;n<t_mc;n++){
		write_mat(n,b,t);
		mc(b,t);
	std::cout<<"t="<<n<<std::endl;	
	print_mat(t);	
	}
		write_mat(n+1,b,t);
	return 0;
}
