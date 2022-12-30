#include "matrix.h"


matrix::matrix()
{
	A.resize(4);
}

matrix matrix::inverse()
{
	double det = A[0] * A[3] - A[1] * A[2];
	if (det == 0) throw("det is null");
	matrix res;
	res.A[0] = A[3] / det;
	res.A[3] = A[0] / det;
	res.A[2] = -A[2] / det;
	res.A[1] = -A[1] / det;
	return res;
}

std::vector<double> matrix::operator*(std::vector<double> vec)
{
	std::vector<double> res(2);
	res[0] = A[0] * vec[0] + A[1] * vec[1];
	res[1] = A[2] * vec[0] + A[3] * vec[1];
	return res;
}
 std::vector <double> matrix::sum(std::vector<double> A, std::vector<double> B)
{
	return { A[0] + B[0],A[1] + B[1] };
}

 std::vector <double> matrix::dif(std::vector<double> A, std::vector<double> B)
{
	return { A[0] - B[0],A[1] - B[1] };
}

 std::vector <double> matrix::pr(std::vector<double> A, double b)
{
	return { A[0] * b, A[1] * b };
}
 double matrix::norm(std::vector<double> B)
 {
	 return sqrt(B[0] * B[0] + B[1] * B[1]);
 }