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
    int n=20, cmax=5000;
    double Ti=1.8,Tf=3.0,Tstep=0.1;
    mat lattice = zeros<mat>(n,n);
    default_random_engine e;
    double E,M,Eavg,EEavg,Mavg,MMavg;
    vector<double> w(17);
    ofstream Efile("avg_energy.txt");
    ofstream EEfile("avg_energy2.txt");
    ofstream Evar("var_energy.txt");
    ofstream Mfile("avg_mag.txt");
    ofstream MMfile("avg_mag2.txt");
    ofstream Mvar("var_mag.txt");
//    ofstream config_T("config_T.txt");
//    ofstream E_mc("energy_MC_rand.txt");
//    ofstream M_mc("mag_MC_rand.txt");
//    ofstream config_mc("config_MC_rand.txt");
//    ofstream P_e("Ei.txt");
    for(double T = Ti;T<=Tf;T+=Tstep){
        if(T>=2.2 && T<2.4){
            Tstep=0.02;
            cmax=100000;
        }
        else{Tstep = 0.1;cmax=5000;}
//    double T=1.0;
//    for(int cmax=500;cmax<=cmax_max;cmax+=500){
        long acc=0;
        E=0;M=0;Eavg=0;EEavg=0;Mavg=0;MMavg=0;
        initialize(n,lattice,T,E,M,w,e,"order");
        for(int cycles = 1;cycles<=cmax;cycles++){
            metropolis(n,lattice,e,E,M,w,acc);
            Eavg+=E;EEavg+=E*E;Mavg+=abs(M);MMavg+=M*M;
//            if(cycles>7000){P_e<<cycles<<" "<<E/((double)n*n)<<endl;}
        }
        Eavg=Eavg/((double)cmax);
        EEavg=EEavg/((double)cmax);
        Mavg=Mavg/((double)cmax);
        MMavg=MMavg/((double)cmax);
        //cout<<"Evar = "<<(EEavg-Eavg*Eavg)/((double)n*n);
        Efile<<T<<" "<<Eavg/((double)n*n)<<endl;
        EEfile<<T<<" "<<EEavg/((double)n*n)<<endl;
        Evar<<T<<" "<<(EEavg-Eavg*Eavg)/((double)n*n*T*T)<<endl;
        Mfile<<T<<" "<<Mavg/((double)n*n)<<endl;
        MMfile<<T<<" "<<MMavg/((double)n*n)<<endl;
        Mvar<<T<<" "<<(MMavg-Mavg*Mavg)/((double)n*n*T)<<endl;
//        config_T<<T<<" "<<acc<<endl;
//        E_mc<<cmax<<" "<<Eavg/((double)n*n)<<endl;
//        M_mc<<cmax<<" "<<Mavg/((double)n*n)<<endl;
//        config_mc<<cmax<<" "<<acc<<endl;
    }

    
    default_random_engine engine;
    vector<pair<int,int> > neigh= nearest_neighbors(0,3,n);

    
    
    return 0;
}



