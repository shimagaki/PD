#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

template <class T>
void print1(std::vector<T> v){
	for(auto v1:v){
	std::cout<<v1<<"  ";
	}
}

template <class T>
void print2(std::vector<std::vector<T> > v){
	int i,j;
	for(i=0; i<v.size();i++){
		j=0;
		for(j=0; j<v[0].size();j++){
			std::cout<<v[i][j]<<"  ";	
		}
		std::cout<<std::endl;
	}}

template <class T>
void seq_init(std::vector<std::vector<T> >& v){
	int i=0, j=0;
	for(i=0; i<v.size();++i){
		for(j=0; j<v[0].size();++j){
			v[i][j] =i*j;	
		}
	}
}


template <class T>
void f1( std::vector<T>& v){
	int it,i;
	for(auto& v1:v){
	v1 *=-1;
	}
	//return v;
}

template <class T>
void f2( std::vector<std::vector<T> >& v){
	for(auto& v1:v){
		for(auto& v2:v1){
			v2 *=-1;	
		}
		std::cout<<std::endl;
	}
}

int main(int argc, char* argv[]){
	double b;
	int i,j;
	int n=5,m=3;
	std::vector<int> v_1d(m,1);	
	std::vector<std::vector<int> >v_2d(m,std::vector<int>(n,1));	
	print2(v_2d);	
	std::cout<<"############"<<std::endl;
	seq_init(v_2d);
	std::cout<<"############"<<std::endl;


	std::cout<<"\n This is f1"<<std::endl;
	print1(v_1d);
	//v_1d=f1(v_1d);
	f1(v_1d);
	print1(v_1d);
	
	std::cout<<"\nv=" <<std::endl;	
	print2(v_2d);	
	
	std::cout<<"\nv =-1" <<std::endl;	
	f2(v_2d);
	print2(v_2d);	
	
	return 0;
}
