#include <cstdlib>
#include </opt/local/include/armadillo.h>
#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include "ising.h"

using namespace std;
using namespace arma;

void initialize(int n, mat& lattice, double T, double& E,double& M,
        vector<double>&w,default_random_engine&e,string order){
    double de=-8;
    int indx;
    int spin;
    for(int i=0;i<=16;i++){
        if(i%4==0){w[i]=(exp(-de/T));}
        else{w[i]=(0);}
        de++;
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(order=="rand"){
                indx=rand_index(e,2);
                spin=rand_spin(indx);
            }
            else{spin=1;}
            lattice(i,j)=spin;
            M+=spin;
        }
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(j+1<n){
               E-=lattice(i,j)*lattice(i,j+1);
               //cout<<"("<<i<<","<<j<<")("<<i<<","<<j+1<<")"<<endl;
            }
            else{
                E-=lattice(i,j)*lattice(i,0);
                //cout<<"("<<i<<","<<j<<")("<<i<<","<<0<<")"<<endl;
            }
        }
    }
    for(int j=0;j<n;j++){
        for(int i=0;i<n;i++){
            if(i+1<n){
                E-=lattice(i,j)*lattice(i+1,j);
               // cout<<"("<<i<<","<<j<<")("<<i+1<<","<<j<<")"<<endl;
            }
            else{
                E-=lattice(i,j)*lattice(0,j);
                //cout<<"("<<i<<","<<j<<")("<<0<<","<<j<<")"<<endl;
            }
        }
    }
}

void metropolis(int n, mat& lattice, default_random_engine& e,double& E,
        double& M,vector<double> w,long& acc){
    int x,y,delta_E;
    vector<pair<int,int> > nn;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            x=rand_index(e,n);
            y=rand_index(e,n);
            nn = nearest_neighbors(x,y,n);
            delta_E=2*lattice(x,y)*(lattice(nn[0].first,nn[0].second)+
                    lattice(nn[1].first,nn[1].second)+lattice(nn[2].first,nn[2].second)+
                    lattice(nn[3].first,nn[3].second));
            if(rand_num(e)<=w[delta_E+8]){
                lattice(x,y)*=-1;
                M+=(double)2*lattice(x,y);
                E+=(double)delta_E;
                acc++;
            }
        }
    }
    
}

vector<pair<int,int> > nearest_neighbors(int i, int j, int n){
    vector<pair<int,int> > neighbors;
    if(i+1==n){
        neighbors.push_back(make_pair(0,j));
    }
    else{
        neighbors.push_back(make_pair(i+1,j));
    }
    if(j+1==n){
        neighbors.push_back(make_pair(i,0));
    }
    else{
        neighbors.push_back(make_pair(i,j+1));
    }
    if(i-1<0){
        neighbors.push_back(make_pair(n-1,j));
    }
    else{
        neighbors.push_back(make_pair(i-1,j));
    }
    if(j-1<0){
        neighbors.push_back(make_pair(i,n-1));
    }
    else{
        neighbors.push_back(make_pair(i,j-1));
    }
    return neighbors;
}

int rand_index(default_random_engine&e,int n){
    uniform_int_distribution<int> dist(0,n-1);
    return dist(e);
}

double rand_num(default_random_engine&e){
    uniform_real_distribution<double> dist(0.0,1.0);
    return dist(e);
}

int rand_spin(int indx){
    vector<int> spins(2);
    spins[0]=-1;
    spins[1]=1;
    return spins[indx];
}