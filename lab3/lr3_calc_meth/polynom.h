#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <initializer_list>
#include"file_helper.h"
#include "C:/Users/Артем/source/repos/calc_method_1/calc_method_1/matrix.h"

class polynom
{
private:
	size_t n;
	std::vector<double> coef;
	
public:
	polynom(size_t n_ = 10);
	polynom(std::vector<double> vec);
	polynom(std::initializer_list<double> list);

	void lagrange_interpol(std::vector<double> nodes, double(*foo)(double), std::vector<double> args);
	double lagrange_interpol(std::vector<double> nodes, double(*foo)(double), double arg);
	void spline_interpol(std::vector<double>& nodes, double(*foo)(double),const std::vector<double>& args);
	double spline_interpol(const std::vector<double>& nodes, double(*foo)(double), double arg);
	double get_prec(const std::vector<double>& nodes, double(*foo)(double), const std::string& name_of_method);
	std::vector<double> progonka(std::vector<std::vector<double>> A);
};

std::vector<std::vector<double>> create_coef_for_spline(const std::vector<double>& vec, double(*foo)(double));