#include "functions.h"
#include <cmath>
double foo1(double x)
{
	return x * x;
}

double foo2(double x)
{
	return 1 / (1 + x * x);
}

double foo3(double x)
{
	return 1 / (atan(1 + 10 * x * x));
}

double foo4(double x)
{
	return exp(sqrt(2) * log(4 * x * x * x + 2 * x * x - 4 * x + 2)) + asin(1 / (5 + x - x * x)) - 5;
}

double foo5(double x)
{
	return exp(x);
}
double foo6(double x)
{
	return 1;
}

double foo7(double x)
{
	return exp(x*log(2));
}