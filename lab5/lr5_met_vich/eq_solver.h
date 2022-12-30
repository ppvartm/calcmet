#pragma once
#include "file_helper.h"
#include <vector>
#include "matrix.h"
#include <iostream>
class eq_solver
{
public:
	
	std::pair<double, double> root_localization(double (*foo)(double), std::pair<double, double> left_right, int N);

	double chord_method(double (*foo)(double), std::pair<double, double> left_right);

	double newton_method(double (*foo)(double),double (*dfoo)(double), std::pair<double, double> left_right);
	double method_sek(double (*foo)(double), std::pair<double, double> left_right);

	std::vector<double> newton_matrix_method(std::vector<double> (*system1)(std::vector<double> ), matrix (*dsistem1)(std::vector<double> ), std::vector<double> x0);
};

