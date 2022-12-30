#include "foo_for_ex.h"

double foo1(double x)
{
	return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}
double dfoo1(double x)
{
	return 0.121495 - 1.5119 * x + 5.9535 * x * x - 9.28 * x * x * x + 5 * x * x * x * x;
}
double dfoo1_err(double x)
{
	double eps = 0.0000000001;
	return(foo1(x + eps) - foo1(x)) / eps;
}

double foo2(double x)
{
	return (sqrt(x+1)-1);
}
double dfoo2(double x)
{
	return 1/(2*sqrt(x+1));
}
double dfoo2_err(double x)
{
	double eps = 0.0000000001;
	return(foo2(x + eps) - foo2(x)) / eps;
}

double foo3(double x)
{
	return 35*x*x*x-67*x*x-3*x+3;
}
double dfoo3(double x)
{
	return 105*x*x-134*x-3;
}
double dfoo3_err(double x)
{
	double eps = 0.0000000001;
	return(foo3(x + eps) - foo3(x)) / eps;
}

double foo4(double x)
{
	return (x - 1) * (x - 1)* (x - 1) * (x - 1)* (x - 1) * (x - 1) * (x - 1) * (x - 1);
}
double dfoo4(double x)
{
	return 8 * (x - 1) * (x - 1) * (x - 1) * (x - 1)* (x - 1) * (x - 1) * (x - 1);
}
double dfoo4_err(double x)
{
	double eps = 0.00001;
	return(foo4(x + eps) - foo4(x)) / eps;
}

std::vector<double> system1(std::vector<double> x_y)
{
	return { x_y[0] * x_y[0] - x_y[1] * x_y[1] - 15,x_y[0] * x_y[1] + 4 };
}

matrix dsistem1(std::vector<double> x_y)
{
	matrix res;
	res.A[0] = 2 * x_y[0];
	res.A[1] = -2 * x_y[1];
	res.A[2] = x_y[1];
	res.A[3] = x_y[0];
	return res;
}
matrix dsistem1_err(std::vector<double> x_y)
{
	double eps = 0.01;
	matrix res;
	res.A[0] = ((x_y[0]+eps) * (x_y[0]+eps) - x_y[1] * x_y[1] - 15- (x_y[0] * x_y[0] - x_y[1] * x_y[1] - 15))/eps;
	res.A[1] = (x_y[0] * x_y[0] - (x_y[1]+eps) * (x_y[1]+eps) - 15 - (x_y[0] * x_y[0] - x_y[1] * x_y[1] - 15)) / eps;;
	res.A[2] = ((x_y[0]+eps) * x_y[1] + 4 - (x_y[0] * x_y[1] + 4))/eps;
	res.A[3] = (x_y[0] * (x_y[1]+eps) + 4 - (x_y[0] * x_y[1] + 4)) / eps;
	return res;
}