#pragma once
#include <vector>
class matrix
{
public:
	std::vector<double> A;
public:
	matrix();
	matrix inverse();
	std::vector<double> operator*(std::vector<double> vec);
	static std::vector <double> sum(std::vector<double> A, std::vector<double> B);
	static std::vector <double> dif(std::vector<double> A, std::vector<double> B);
	static std::vector <double> pr(std::vector<double> A, double b);
	static double norm(std::vector<double> B);
};