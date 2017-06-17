#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;
int N = 8;

vector <int> set(N,-1); // Each element remember their parent's index.
//	Root rememter their "(-1) x depth", remember every root have negative.

int get_root(int x){
	int root_x;
	if(set[x]<0){ root_x = x; }
	if(set[x]>=0){ root_x = get_root(set[x]); }
	return root_x;
}

void myunion(int root_x, int root_y){
		/* Recal that set[root] remember their depth.*/
		/* Longer root is more robast.*/
		if( (-set[root_x]) > (-set[root_y]) ){
			set[root_y] = root_x;
			// "-set[root_x]" does not change.
		}
		if( (-set[root_x]) < (-set[root_y]) ){
			set[root_x] = root_y;	
			// "-set[root_x]" does not change.
		}	
		if( (-set[root_x]) == (-set[root_y]) ){
			set[root_x] -= 1; // Added 1 to the depth of new root.
			set[root_y] = root_x;
		}
}

bool myfind(int x, int y){
	bool flag=false;
	int root_x, root_y;
	root_x = get_root(x);	
	root_y = get_root(y);
	if(root_x != root_y){
		myunion(root_x, root_y);	
		flag = true;
	}
	return flag;
}

int main(){
	myfind(0,1);	
	myfind(2,3);	
	myfind(4,5);	
	myfind(6,7);
	myfind(1,0);
	myfind(1,3);
	myfind(4,3);
	myfind(5,3);
	
	for(int j=0;j<N;++j){
	cout<<j<<"\t"<<set[j]<<endl;
	}	

	return 0;
}
