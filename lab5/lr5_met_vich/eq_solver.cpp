#include "eq_solver.h"

std::pair<double, double> eq_solver::root_localization(double (*foo)(double), std::pair<double, double> left_right, int N)
{
	double step = (left_right.second - left_right.first)/N;
	for (int i = 0; i < N; ++i)
	{
		if (foo(left_right.first + step * i) * foo(left_right.first + step * (i + 1)) <= 0)
			return {left_right.first + step * i, left_right.first + step * (i + 1) };
	}
	return left_right;
	throw("na otrezke net kornei");
}



double eq_solver::chord_method(double (*foo)(double), std::pair<double, double> left_right)
{
	return (foo(left_right.first) * left_right.second - foo(left_right.second) * left_right.first)
		/ (foo(left_right.first) - foo(left_right.second));   
}

double eq_solver::newton_method(double (*foo)(double),double (*dfoo)(double), std::pair<double, double> left_right)
{
	double x0, x = chord_method(foo, left_right);
	double eps = 0.000000001;
	int i = 0;
	std::vector<double> temp;	
	do {
//		std::cout << "x = " << x << std::endl;
		x0 = x;
		temp.push_back(x0);
		++i;
		x = x0 - foo(x0) / dfoo(x0);
		if(i>1)
		//std::cout << "p =  " << log(fabs(x - 1.) / fabs(x0 - 1.)) / log(fabs(x0 - 1.) / fabs(temp[i-2] - 1.)) << std::endl;
		if (x > left_right.second)
		{
			x = chord_method(foo, { left_right.first,x });
		}
		if (x < left_right.first)
		{
			x = chord_method(foo, { x, left_right.second });
		}
	} while (fabs(x0 - x) > eps);
	return x0;
}

double eq_solver::method_sek(double (*foo)(double), std::pair<double, double> left_right)
{
	double eps = 0.0000001;
	//double a = left_right.first, b = left_right.second;
	//int i = -1;
	//std::vector<double> temp;
	//while (fabs(b - a) > eps) {
	//	a = a - (b - a) * foo(a) / (foo(b) - foo(a));
	//	b = b - (a - b) * foo(b) / (foo(a) - foo(b));
	//	temp.push_back(b); ++i;
	//	if(i>1) 
	//		std::cout << "p =  " << log(fabs(temp[i] - 0.22) / fabs(temp[i-1] - 0.22)) / log(fabs(temp[i-1] - 0.22) / fabs(temp[i - 2] - 0.22)) << std::endl;
	//}
	double x_last = left_right.first;
	double x_last_last = left_right.second;
	double x=DBL_MAX;
	int i = -1;
	std::vector<double> temp;
	while (true)
	{
		x = x_last - foo(x_last) / (foo(x_last) - foo(x_last_last)) * (x_last - x_last_last);
		temp.push_back(x); ++i;
		if (i > 1)
		std::cout << "p =  " << log(fabs(temp[i] - 0.) / fabs(temp[i-1] - 0.)) / log(fabs(temp[i-1] - 0.) / fabs(temp[i - 2] - 0.)) << std::endl;
		if (fabs(x - x_last) < eps) break;
		x_last_last = x_last;
		x_last = x;
	}
	return x;
}

std::vector<double> eq_solver::newton_matrix_method(std::vector<double>(*system1)(std::vector<double>), matrix(*dsistem1)(std::vector<double>), std::vector<double> x0)
{
	
	double eps = 0.000000001;
	std::vector<double> X0=x0;
	std::vector<double> X;
	matrix J = dsistem1(X0);
	std::vector<double> temp = matrix::pr(system1(X0), -1);
	std::vector<double> Y;
	Y = J.inverse() * temp;
	X = matrix::sum(Y,X0);
	int k = 0;
	while (matrix::norm(matrix::dif(X, X0)) > eps)
	{
		X0 = X;
    	    J = dsistem1(X0);
	    Y = J.inverse() * matrix::pr(system1(X0), -1);
		X = matrix::sum(Y, X0);
		//std::cout << ++k << std::endl;
		++k;
	}
//	std::cout << k << std::endl;
	X.push_back(k);
	return X;
}