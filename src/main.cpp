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
    int n=40, cmax=5000;
    double Ti=1.5,Tf=3.,Tstep=0.1;
    mat lattice = zeros<mat>(n,n);
    default_random_engine e;
    double E,M,Eavg,EEavg,Mavg,MMavg;
    vector<double> w(17);
    ofstream Efile("avg_energy.txt");
    ofstream EEfile("avg_energy2.txt");
    ofstream Mfile("avg_mag.txt");
    ofstream MMfile("avg_mag2.txt");
    for(double T = Ti;T<=Tf;T+=Tstep){
        E=0;M=0;Eavg=0;EEavg=0;Mavg=0;MMavg=0;
        initialize(n,lattice,T,E,M,w);
        for(int cycles = 1;cycles<=cmax;cycles++){
            metropolis(n,lattice,e,E,M,w);
            Eavg+=E;EEavg+=E*E;Mavg+=abs(M);MMavg+=M*M;
        }
        Eavg=Eavg/((double)cmax);
        EEavg=EEavg/((double)cmax);
        Mavg=Mavg/((double)cmax);
        MMavg=MMavg/((double)cmax);
        Efile<<T<<" "<<Eavg/((double)n*n)<<endl;
        EEfile<<T<<" "<<EEavg/((double)n*n)<<endl;
        Mfile<<T<<" "<<Mavg/((double)n*n)<<endl;
        MMfile<<T<<" "<<MMavg/((double)n*n)<<endl;
    }

    
    default_random_engine engine;
    vector<pair<int,int> > neigh= nearest_neighbors(0,3,n);

    
    
    return 0;
}



