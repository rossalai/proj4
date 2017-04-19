#include "catch.hpp"
#include "ising.h"

using namespace std;

TEST_CASE("Testing initial energy/magnetization"){
    int n=2;
    mat lattice = zeros<mat>(n,n);
    vector<double> w(17);
    default_random_engine e;
    double E=0,M=0,T=1;
    initialize(n,lattice,T,E,M,w,e,"order");
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
    int n=2;
    long acc=0;
    double T=1.0,cmax=1000000;
    mat lattice = zeros<mat>(n,n);
    default_random_engine e;
    double E,M,Eavg,EEavg,Mavg,MMavg;
    vector<double> w(17);
    E=0;M=0;Eavg=0;EEavg=0;Mavg=0;MMavg=0;
    initialize(n,lattice,T,E,M,w,e,"order");
    for(int cycles = 1;cycles<=cmax;cycles++){
        metropolis(n,lattice,e,E,M,w,acc);
        Eavg+=E;EEavg+=E*E;Mavg+=abs(M);MMavg+=M*M;
    }
    Eavg=Eavg/((double)cmax);
    EEavg=EEavg/((double)cmax);
    double Z = 2*exp(-8/T)+2*exp(8/T)+12;
    double Eavg_exp=(-16*exp(8/T)+16*exp(-8/T))/Z;
    double EEavg_exp=(128*exp(-8/T)+128*exp(8/T))/Z;
    double Cv = (EEavg-Eavg*Eavg)/(T*T);
    double Cv_exp = (EEavg_exp-Eavg_exp*Eavg_exp)/(T*T);
    REQUIRE(Eavg==Approx(Eavg_exp).epsilon(0.001));
    REQUIRE(EEavg==Approx(EEavg_exp).epsilon(0.001));
    REQUIRE(Cv==Approx(Cv_exp).epsilon(0.01));
    
    Mavg=Mavg/((double)cmax);
    MMavg=MMavg/((double)cmax);
    double Mavg_exp=(8*exp(8/T)+8)/Z;
    double MMavg_exp=(32*exp(8/T)+32)/Z;
    double X = (MMavg - abs(Mavg*Mavg))/T;
    double X_exp = (MMavg_exp - abs(Mavg_exp*Mavg_exp))/T;
    REQUIRE(Mavg==Approx(Mavg_exp).epsilon(0.001));
    REQUIRE(MMavg==Approx(MMavg_exp).epsilon(0.001));
    REQUIRE(X==Approx(X_exp).epsilon(0.05));
    
//    cout<<"<E>calc = "<<Eavg<<" <E>exp = "<<Eavg_exp<<" % diff = "<<100*(Eavg-Eavg_exp)/Eavg_exp<<endl;
//    cout<<"<E^2>calc = "<<EEavg<<" <E^2>exp = "<<EEavg_exp<<" % diff = "<<100*(EEavg-EEavg_exp)/EEavg_exp<<endl;
//    cout<<"Cvcalc = "<<Cv<<" Cvexp = "<<Cv_exp<<" % diff = "<<100*(Cv-Cv_exp)/Cv_exp<<endl;
//    cout<<"<M>calc = "<<Mavg<<" <M>exp = "<<Mavg_exp<<" % diff = "<<100*(Mavg-Mavg_exp)/Mavg_exp<<endl;
//    cout<<"<M^2>calc = "<<MMavg<<" <M^2>exp = "<<MMavg_exp<<" % diff = "<<100*(MMavg-MMavg_exp)/MMavg_exp<<endl;
//    cout<<"Xcalc = "<<X<<" Xexp = "<<X_exp<<" % diff = "<<100*(X-X_exp)/X_exp<<endl;
    
}



