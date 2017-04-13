/* 
 * File:   main.cpp
 * Author: alaina
 *
 * Created on April 13, 2017, 12:06 PM
 */

#include <cstdlib>
#include </opt/local/include/armadillo.h>
#include <iostream>
#include <random>
#include "ising.h"


using namespace std;
using namespace arma;

/*
 * 
 */
int main(int argc, char** argv) {
    int n=4;
    mat lattice = zeros<mat>(n,n);
    double E=0,M=0,T=1;
    vector<double> w;
    initialize(n,lattice,T,E,M,w);
    for(auto element: w){
        cout<<element<<" ";
    }
    
    default_random_engine engine;
    vector<pair<int,int> > neigh= nearest_neighbors(0,3,n);

    
    
    return 0;
}



