/* 
 * File:   ising.h
 * Author: alaina
 *
 * Created on April 13, 2017, 1:16 PM
 */

#ifndef ISING_H
#define	ISING_H

#include <cstdlib>
#include </opt/local/include/armadillo.h>
#include <iostream>
#include <vector>
#include <random>
#include <utility>

using namespace std;
using namespace arma;

void initialize(int, mat&, double, double&,double&,vector<double>&);
int rand_index(default_random_engine&,int);
vector<pair<int,int> > nearest_neighbors(int , int , int );

#endif	/* ISING_H */

