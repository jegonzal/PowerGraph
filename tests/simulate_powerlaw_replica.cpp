#include <iostream>
#include <math.h>
using namespace std;

double h(int num,double alpha){
	double sum=0;
	for(int i=1;i<=num;i++){
		sum+=pow(i,-alpha);
	}
	return sum;
}
double random_replication(int V,double alpha,int p){
	double E=h(V,alpha-1)/h(V,alpha)*V;
	double tmp=V/h(V,alpha);
	double sum=0;
	for(int i=1;i<=V;i++){
		double num_replica=p*(1-pow((1-1.0/p),i+E/V));
		num_replica+=(p-num_replica)/p;
		sum+=tmp * pow(i,-alpha)* num_replica;
	}
	return sum/V;
}

double grid_replication(int V,double alpha,int p){
	double E=h(V,alpha-1)/h(V,alpha)*V;
	double fp=2*pow(p,0.5)-1;
	double tmp=V/h(V,alpha);
	double sum=0;
	for(int i=1;i<=V;i++){
		double num_replica=fp*(1-pow((1-1.0/fp),i+E/V));
		num_replica+=(fp-num_replica)/fp;
		sum+=tmp * pow(i,-alpha)* num_replica;
	}
	return sum/V;
}
double hybrid_replication(int V,double alpha,int p,int threshold){
	double E=h(V,alpha-1)/h(V,alpha)*V;
	double tmp=V/h(V,alpha);
	double sum=0;
	double R_E_H=1-h(threshold,alpha-1)/h(V,alpha-1);
	for(int i=1;i<=V;i++){
		double num_replica;
		if(i<=threshold){
			num_replica=p*(1-pow((1-1.0/p),E/V*(1-R_E_H)));
		} else {
			num_replica=p*(1-pow((1-1.0/p),i+E/V*(1-R_E_H)));
		}
		num_replica+=(p-num_replica)/p;
		sum+=tmp * pow(i,-alpha)* num_replica;
	}
	return sum/V;
}

int main(){
	int V=10000000;//10m
	cout<<"alpha\trandom\tgrid\thybrid\tp=48"<<endl;
	for(double alpha=1.8;alpha<=2.2;alpha+=0.1){
		cout<<alpha<<"\t"
			<<random_replication(V,alpha,48)<<"\t"
			<<grid_replication(V,alpha,48)<<"\t"
			<<hybrid_replication(V,alpha,48,100)<<endl;
	}

	cout<<"p\t";
	for(double alpha=1.8;alpha<=2.2;alpha+=0.1){
		cout<<alpha<<"\t";
	}
	cout<<"grid/hybrid"<<endl;
	for(int i=1;i<=15;i++){
		cout<<i*i<<"\t";
		for(double alpha=1.8;alpha<=2.2;alpha+=0.1){
			cout<<grid_replication(V,alpha,i*i)/hybrid_replication(V,alpha,i*i,100)<<"\t";
		}
		cout<<endl;
	}

	cout<<"p\t";
	for(double alpha=1.8;alpha<=2.2;alpha+=0.1){
		cout<<alpha<<"\t";
	}
	cout<<"random/hybrid"<<endl;
	for(int i=1;i<=15;i++){
		cout<<i*i<<"\t";
		for(double alpha=1.8;alpha<=2.2;alpha+=0.1){
			cout<<random_replication(V,alpha,i*i)/hybrid_replication(V,alpha,i*i,100)<<"\t";
		}
		cout<<endl;
	}
}
