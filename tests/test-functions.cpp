#include "catch.hpp"
#include "ising.h"

using namespace std;

TEST_CASE("Testing initial energy/magnetization"){
    int n=2;
    mat lattice = zeros<mat>(n,n);
    vector<double> w(17);
    double E=0,M=0,T=1;
    initialize(n,lattice,T,E,M,w);
    REQUIRE(E==Approx(-8));
    REQUIRE(M==Approx(4));
}

TEST_CASE("Testing nearest neighbors function"){
    vector<pair<int,int> > neigh= nearest_neighbors(0,3,4);
    REQUIRE(neigh[0].first==Approx(1));
    REQUIRE(neigh[0].second==Approx(3));
    REQUIRE(neigh[1].first==Approx(0));
    REQUIRE(neigh[1].second==Approx(0));
    REQUIRE(neigh[2].first==Approx(3));
    REQUIRE(neigh[2].second==Approx(3));
    REQUIRE(neigh[3].first==Approx(0));
    REQUIRE(neigh[3].second==Approx(2));
}

TEST_CASE("Testing E,E^2,Cv,|M|,M^2,X"){
    int n=4;
    double T=1,cmax=5000;
    mat lattice = zeros<mat>(n,n);
    default_random_engine e;
    double E,M,Eavg,EEavg,Mavg,MMavg;
    vector<double> w(17);
    E=0;M=0;Eavg=0;EEavg=0;Mavg=0;MMavg=0;
    initialize(n,lattice,T,E,M,w);
    for(int cycles = 1;cycles<=cmax;cycles++){
        metropolis(n,lattice,e,E,M,w);
        Eavg+=E;EEavg+=E*E;Mavg+=M;MMavg+=M*M;
    }
    Eavg=Eavg/((double)cmax*n);
    EEavg=EEavg/((double)cmax*n*n);
    double Z = 2*exp(-8/T)+2*exp(8/T)+12;
    double Eavg_exp=(-16*exp(8/T)+16*exp(-8/T))/Z;
    double EEavg_exp=(128*exp(-8/T)+128*exp(8/T))/Z;
    REQUIRE(Eavg==Approx(Eavg_exp).epsilon(0.01));
    REQUIRE(EEavg==Approx(EEavg_exp).epsilon(0.01));
    REQUIRE((EEavg-Eavg*Eavg)/(T*T)==Approx((EEavg_exp-Eavg_exp*Eavg_exp)/(T*T)).epsilon(0.1));
    
    Mavg=Mavg/((double)cmax*n);
    MMavg=MMavg/((double)cmax*n*n);
    double Mavg_exp=(8*exp(8/T)+8)/Z;
    double MMavg_exp=(32*exp(8/T)+32)/Z;
    REQUIRE(Mavg==Approx(Mavg_exp).epsilon(0.01));
    REQUIRE(MMavg==Approx(MMavg_exp).epsilon(0.01));
    REQUIRE((MMavg-Mavg*Mavg)/T==Approx((MMavg_exp-Mavg_exp*Mavg_exp)/T).epsilon(0.05));
    
}

