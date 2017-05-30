#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#define L  64 
#define rand01 (double)rand()/RAND_MAX
// The size of the system is L
// The interaction is next nearest interaction.
// t[x][y] is an angle, theta_xy, where xy is location on the lattice.

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

template <class T>
void random_vect( std::vector<std::vector<T> >& v){
	int i,j;
	for(i=0;i<v.size();i++){
		for(j=0;j<v[0].size();j++){
			v[i][j] = rand01*2.0*M_PI;
		}
	}
}

template <class T>
void cavity_vector( std::vector<std::vector<T> >& v){
	int i,j,Lx,Ly;
	double a,b;
	Lx = v.size();
	Ly = v[0].size();
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			a = pow((j-0.25*Ly)*(j-0.25*Ly)+(i-0.25*Lx)*(i-0.25*Lx),0.5);
			b = pow((j-0.75*Ly)*(j-0.75*Ly)+(i-0.75*Lx)*(i-0.75*Lx),0.5);
			v[i][j]= -3;//((j-0.25*Ly)/(i-0.25*Lx))/a
				//	+ ((j-0.75*Ly)/(i-0.75*Lx))/b;
		}
	}
}

double  erg(double t_s,double t_sr,double t_sl,double t_su,double t_sd){
	return  - (cos(t_s-t_sr)+cos(t_s-t_sl)
			+cos(t_s-t_su)+cos(t_s-t_sd));
}

double seaqHB(double b,double t_s,double t_sr,double t_sl,double t_su,double t_sd){
	double u, dt, dE_prev,dE_prop,r;
	u = rand01;
	//dt is proposal of the theta.	
	dt = 2.0*M_PI*rand01/8.0;	
	dE_prev = erg(t_s,t_sr,t_sl,t_su,t_sd);	
	dE_prop = erg(t_s+dt,t_sr,t_sl,t_su,t_sd);	
	r = exp(-b*(dE_prop-dE_prev));
	r = r; 
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
			//t[i][j]=t_ij;
			if(t_ij>M_PI){
			t_ij -=2*M_PI;
			}else if(t_ij<-M_PI){
			t_ij +=2*M_PI;
			}else if(t_ij<=M_PI && t_ij>=-M_PI){
			t[i][j]=t_ij;
			}
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
	std::string fname="xy-tMC"+std::to_string(n)+".dat";
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
	cavity_vector(t);
	int i,j,n,t_mc, t_th=1000;
	std::cout<<"1/T="<<std::endl;
	std::cin>>b;	
	srand(0);
	std::cout<<"t_mc="<<std::endl;
	std::cin>>t_mc;	
	
	for(n=0;n<t_mc+t_th;n++){
		if(n>t_th || n<50){
		write_mat(n,b,t);
		}
		mc(b,t);
	std::cout<<"t="<<n<<std::endl;	
	//print_mat(t);	
	}
		write_mat(n+1,b,t);
	return 0;
}
