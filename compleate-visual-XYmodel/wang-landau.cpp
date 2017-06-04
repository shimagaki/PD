#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;

const int N = 1000;
const int D=2;
vector <double> z(2);
vector <double> E_vec(N,0);
vector <double> S_vec(N,0);
vector <double> g_vec(N,1); // It is stable to be log of g.
vector <double> Prob_hist(N,0);
double T;
double E_max=100.0,E_min=0;
double dz=(E_max-E_min)/N;
double Prob_max, Prob_min;
inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

double calc_E(vector<double> v){
	double erg=0;
	for(int d=0;d<D;++d){
	erg+=v[d]*v[d];
	}
	return 0.5 * erg;
}

void init_erg(){
	double dE=(E_max-E_min)/N;
	for(int n=0;n<N;++n){
	E_vec[n] = dE*n;	
	}
}

int index_of_E(double erg){
	int index=0;
	if(erg < E_min || erg > E_max){
		cout<<"Not found"<<endl;
	}else{
	//index = it-E_vec.begin()-1;
	//auto pos = distance(E_vec.begin(), it);
	vector<double>::iterator it = upper_bound(E_vec.begin(), E_vec.end(),erg);
	index = (distance(E_vec.begin(), it)-1); // "-1" is crucial point.
	}
	return index;	
}
/*
void  MH(){
	double erg,erg_prop;
	vector<double> z_prop(D);
	double u ,r;
	erg = calc_E(z);
	for(int d=0;d<D;++d){
	z_prop[d] = z[d]+ (std_rand()-0.5)*std_rand();}
	erg_prop = calc_E(z_prop);
	u = std_rand();
	r = exp(-(erg_prop-erg)/T); 
	if(u<r){
	z = z_prop;	
	}
}
*/
void  Wang_Landau(double eps,double flatness){
	double f = exp(1);
	double erg,erg_prop;
	vector<double> z_prop(D);
	double g,g_prop;
	double u ,r;
	double dS;
	int index_g,index_g_prop;
	double Prob_max, Prob_min;
	int stop_f=1000,stop_hist=1000;
        int count_f=0,count_hist=0;
	int min_hist;double ave_hist;
	int n_update=0;
	while(f>(1.0+eps)){
		for(int d=0;d<D;++d){
		z_prop[d] = z[d]+ (std_rand()-0.5)*std_rand();
		}
		erg=calc_E(z); 
		erg_prop=calc_E(z_prop);	
		if(E_min<erg_prop && erg_prop<E_max){
			n_update+=1;
			index_g=index_of_E(erg);	
			index_g_prop=index_of_E(erg_prop);
			u = std_rand();
			//This is a heart of Wang_Landau algo.
			//r = g_vec[index_g]/g_vec[index_g_prop]; 
			dS = S_vec[index_g]-S_vec[index_g_prop];
			if(dS>0){
			r=1;	
			}else{
			r = exp(dS); 
			}
			if(u<r){
			z = z_prop;
			Prob_hist[index_g_prop]+=1;	
			S_vec[index_g_prop]+=log(f);
			//g_vec[index_g_prop]*=f;	
			}else {
			Prob_hist[index_g]+=1;	
			S_vec[index_g]+=log(f);
			//g_vec[index_g]*=f;	
			}
			//Check the criteria at very 100-th step.	
			if(count_hist%100==0){
			ave_hist = accumulate(Prob_hist.begin(), Prob_hist.end(),0)/Prob_hist.size();
			min_hist = *min_element(Prob_hist.begin(), Prob_hist.end());
			if(min_hist>(ave_hist*flatness)){
				//RESET
				for(int i=0;i<N;i++){Prob_hist[i]=0;}
				f = sqrt(f);
				cout<<"level-f="<<f<<endl;
			}
			}
			count_hist+=1;	
		}
	}
	cout <<"update="<<n_update <<endl;
	//for(int i=0;i<N;i++){S_vec[i]=1.0/n_update;}
}

int main(){
	init_erg();
	double eps = 0.000001;
	double flatness = 0.8;
	Wang_Landau(eps,flatness);
	cout << "*****n,S_vec, g_vec, Prob_hist*******" <<endl;
	for(int n=0;n<N;++n){
	//printf("%d %2.2f %2.2f %2.2f\n",n,S_vec[n],g_vec[n],Prob_hist[n]);
	printf("%d %8.7f\n",n,S_vec[n]);
	}

}
