#include "catch.hpp"
#include "ising.h"

using namespace std;

TEST_CASE("Testing initial energy/magnetization"){
    int n=2;
    mat lattice = zeros<mat>(n,n);
    double E=0,M=0,T=1;
    initialize(n,lattice,T,E,M);
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

