#pragma once
#include <cmath>
#include <vector>
#include"matrix.h"

double foo1(double x);
double dfoo1(double x);
double dfoo1_err(double x);

double foo2(double x);
double dfoo2(double x);
double dfoo2_err(double x);

double foo3(double x);
double dfoo3(double x);
double dfoo3_err(double x);

double foo4(double x);
double dfoo4(double x);
double dfoo4_err(double x);

std::vector<double> system1(std::vector<double> x_y);
matrix dsistem1(std::vector<double> x_y);
matrix dsistem1_err(std::vector<double> x_y);