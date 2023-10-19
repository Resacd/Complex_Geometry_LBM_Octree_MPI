#pragma once
#include <vector>
#include <iostream>
#include "Struct.h"


int gridx = 25;
double U = 0.1;
double Re = 100;
double niu = U * gridx / Re;
double tau = 3.0 * niu + 0.5;
double spacing = 0.1;
int mdepth = 0;
std::vector<double> inVel{ 1.0,0.0,0.0 };
std::vector<double> InOutCoord{ 0.0,2.0 };
//std::vector<double> InOutCoord{ -192.514,-157.027 };
double pressure = 1.0;