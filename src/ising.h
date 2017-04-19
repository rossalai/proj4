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
#include <string>

using namespace std;
using namespace arma;

void initialize(int, mat&, double, double&,double&,vector<double>&,default_random_engine&,string);
int rand_index(default_random_engine&,int);
double rand_num(default_random_engine&e);
int rand_spin(int );
vector<pair<int,int> > nearest_neighbors(int , int , int );
void metropolis(int, mat& , default_random_engine&,double&,double& ,vector<double>,long&);


#endif	/* ISING_H */

