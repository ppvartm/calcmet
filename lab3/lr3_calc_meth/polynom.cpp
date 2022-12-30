#include "polynom.h"


polynom::polynom(size_t n_) :n(n_)
{}

polynom::polynom(std::vector<double> vec) :n(vec.size()), coef(vec)
{}

polynom::polynom(std::initializer_list<double> list) :n(list.size()), coef(list)
{}

void polynom::lagrange_interpol(std::vector<double> nodes, double(*foo)(double), std::vector<double> args)
{
	std::vector<double> values_in_nodes(nodes.size());
	for (size_t i = 0; i < nodes.size(); ++i)
		values_in_nodes[i] = foo(nodes[i]);
	file_helper file;
	std::vector<double> c(values_in_nodes.size());
	for (size_t i = 0; i < args.size(); ++i)
	{
		double values_of_polynom = 0;
		c = std::vector<double>(values_in_nodes.size());
		for (size_t k = 0; k < c.size(); ++k)
		{
			c[k] = 1;
			for (size_t j = 0; j < c.size(); ++j)
			{
				if (k == j) continue;
				c[k] *= ((args[i] - nodes[j]) / (nodes[k] - nodes[j]));
			}
			values_of_polynom += c[k] * values_in_nodes[k];
		}
		file.write({ args[i],values_of_polynom });
	}
}

double polynom::lagrange_interpol(std::vector<double> nodes, double(*foo)(double), double arg)
{

	std::vector<double> values_in_nodes(nodes.size());
	for (size_t i = 0; i < nodes.size(); ++i)
		values_in_nodes[i] = foo(nodes[i]);
	std::vector<double> c(values_in_nodes.size());
	double values_of_polynom = 0;
	c = std::vector<double>(values_in_nodes.size());
	for (size_t k = 0; k < c.size(); ++k)
	{
		c[k] = 1;
		for (size_t j = 0; j < c.size(); ++j)
		{
			if (k == j) continue;
			c[k] *= ((arg - nodes[j]) / (nodes[k] - nodes[j]));
		}
		values_of_polynom += c[k] * values_in_nodes[k];
	}
	return values_of_polynom;
}

void polynom::spline_interpol( std::vector<double>& nodes, double(*foo)(double), const std::vector<double>& args)
{
	std::vector<double> c = matrix<double>(base_type_of_matrix::singular,1).get_decision_zeid_for_special_matrix(create_coef_for_spline(nodes, foo));
	std::vector<double> a(nodes.size() - 1);
	for (size_t i = 0; i < a.size(); ++i)
		a[i] = foo(nodes[i]);
	std::vector<double> h(nodes.size());
	std::vector<double> g(nodes.size());
	h[0] = 0;
	for (size_t i = 0; i < nodes.size()-1; ++i)
		h[i] = nodes[i+1] - nodes[i];
	for (size_t i = 0; i < nodes.size()-1; ++i)
		g[i] = (foo(nodes[i+1]) - foo(nodes[i ])) / h[i];
	std::vector<double> b(nodes.size() - 1);
	for (size_t i = 0; i < b.size(); ++i)
		b[i] = g[i] - (c[i + 1] + 2*c[i]) * h[i] / 3;
	std::vector<double> d(nodes.size() - 1);
	for (size_t i = 0; i < d.size(); ++i)
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);

	file_helper file;
	size_t j = 0;
	std::sort(nodes.begin(), nodes.end());
	for (size_t i = 1; i < nodes.size(); ++i)
	{
		while (args[j] <= nodes[i])
		{
			double answ = a[i-1] + b[i-1] * (args[j] - nodes[i-1]) +
				c[i-1] * (args[j] - nodes[i-1]) * (args[j] - nodes[i-1]) +
				d[i-1] * (args[j] - nodes[i-1]) * (args[j] - nodes[i-1]) * (args[j] - nodes[i-1]);
			file.write({ args[j++], answ });
			if (j == args.size()) break;
		}
	}

}
double polynom::spline_interpol(const std::vector<double>& nodes, double(*foo)(double), double arg)//возращает значения сплайна в точке
{
	//вектор коэфов c
  // std::vector<double> c = matrix<double>(base_type_of_matrix::singular, 1).get_decision_zeid_for_special_matrix(create_coef_for_spline(nodes, foo));
	std::vector<double> c = progonka(create_coef_for_spline(nodes, foo));
	//значения в узлах
	std::vector<double> a(nodes.size() - 1);
	for (size_t i = 0; i < a.size(); ++i)
		a[i] = foo(nodes[i]);
	std::vector<double> h(nodes.size());
	std::vector<double> g(nodes.size());
	h[0] = 0;
	for (size_t i = 0; i < nodes.size() - 1; ++i)
		h[i] = nodes[i + 1] - nodes[i];
	for (size_t i = 0; i < nodes.size() - 1; ++i)
		g[i] = (foo(nodes[i + 1]) - foo(nodes[i])) / h[i];
	std::vector<double> b(nodes.size() - 1);
	for (size_t i = 0; i < b.size(); ++i)
		b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
	std::vector<double> d(nodes.size() - 1);
	for (size_t i = 0; i < d.size(); ++i)
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	size_t j = 0;
	double answ = 0;
	for (size_t i = 1; i < nodes.size(); ++i)
	{
		if ((arg <= nodes[i])&&(arg >= nodes[i-1]))
			 answ = a[i - 1] + b[i - 1] * (arg - nodes[i - 1]) +
				c[i - 1] * (arg - nodes[i - 1]) * (arg - nodes[i - 1]) +
				d[i - 1] * (arg - nodes[i - 1]) * (arg - nodes[i - 1]) * (arg - nodes[i - 1]);
	}
	return answ;
}

double polynom::get_prec(const std::vector<double>& nodes, double(*foo)(double), const std::string& name_of_method)
{
	int N = 10000;
	if (name_of_method == "spline") {
		double max = DBL_MIN;
		for (size_t i = 0; i < N + 1; i++)
		{
			double temp = fabs(spline_interpol(nodes, foo, nodes[0] + (nodes[nodes.size() - 1] - nodes[0]) / N * i) - foo(nodes[0] + (nodes[nodes.size() - 1] - nodes[0]) / N * i));
			if (temp > max)
				max = temp;
		}
		return max;
	}
	else if (name_of_method == "lagrange")
	{
		double max = DBL_MIN;
		for (size_t i = 0; i < N + 1; i++)
		{
			double temp = fabs(lagrange_interpol(nodes, foo, nodes[0] + (nodes[nodes.size() - 1] - nodes[0]) / N * i) - foo(nodes[0] + (nodes[nodes.size() - 1] - nodes[0]) / N * i));
			if (temp > max)
				max = temp;
		}
		return max;
	} return 10000000;
}

std::vector<std::vector<double>> create_coef_for_spline(const std::vector<double>& x, double(*foo)(double))
{
	std::vector<double> a(x.size() - 1);
	std::vector<double> b(x.size() - 1);
	std::vector<double> c(x.size() - 1);
	std::vector<double> d(x.size() - 1);
	std::vector<double> h(x.size());
	std::vector<double> g(x.size());
	h[0] = 0;
	for (size_t i = 1; i < x.size(); ++i)
		h[i] = x[i] - x[i - 1];
	for (size_t i = 1; i < x.size(); ++i)
		g[i] = (foo(x[i]) - foo(x[i - 1])) / h[i];
	//a[0] = b[0] = c[0] = d[0] = 0.;
	for (size_t i = 2; i < x.size(); ++i)
	{
		a[i - 1] = h[i - 1];
		b[i - 1] = 2 * (h[i - 1] + h[i]);
		c[i - 1] = h[i];
		d[i - 1] = 3 * (g[i] - g[i - 1]);
	}
	return { a,b,c,d };
}
//std::vector<double> polynom::progonka(std::vector<std::vector<double>> A)
//{
//int m = A[0].size();
//vector<vector<double>> koef_progonki(m - 1);
//koef_progonki[m - 2].resize(2);
//koef_progonki[m - 2][0] = A[2][m - 1] / (-A[1][m - 1]);
//koef_progonki[m - 2][1] = -A[3][m - 1] / (-A[1][m - 1]);
//double temp = 0;
//for (int i = m - 3; i > 0; --i)
//{
//	koef_progonki[i].resize(2);
//	temp = -A[1][i + 1] - A[0][i + 1] * koef_progonki[i + 1][0];
//	koef_progonki[i][0] = A[2][i + 1] / temp;
//	koef_progonki[i][1] = (-A[3][i + 1] + A[0][i + 1] * koef_progonki[i + 1][1]) / temp;
//}
//vector<double> X(m - 1);
//X[0] = (-A[3][1] + A[0][1] * koef_progonki[1][1]) / (-A[1][1] - A[0][1] * koef_progonki[1][0]);
//for (int i = 1; i < m - 2; ++i)
//	X[i] = koef_progonki[i][0] * X[i - 1] + koef_progonki[i][1];
//X[m - 2] = koef_progonki[m - 2][0] * X[m - 3] + koef_progonki[m - 2][1];
//
//	return X;
//}

std::vector<double> polynom::progonka(std::vector<std::vector<double>> A)
{
	std::vector<double> a = A[0];
	std::vector<double> b = A[1];
	std::vector<double> c = A[2];
	std::vector<double> d = A[3];

	std::for_each(b.begin(), b.end(), [](double& n) { n*=-1; });
	std::for_each(d.begin(), d.end(), [](double& n) { n *= -1; });

	size_t n = a.size()-1;

	std::vector<double> alpha(n+2);
	std::vector<double> beta(n+2);

	alpha[2] =  c[1] / b[1];
	beta[2] = d[1] / b[1];

	for (size_t i = 2; i <= n-1; ++i)
	{
		alpha[i+1] = c[i] / (b[i] - a[i] * alpha[i]);
		beta[i+1] = (d[i]+a[i]*beta[i]) / (b[i] - a[i] * alpha[i]);
	}
	beta[n + 1] = (d[n] + a[n] * beta[n]) / (b[n] - a[n] * alpha[n]);
	std::vector<double> x(n+1);
	x[n] = beta[n+1];

	for (size_t i = n - 1; i >= 1; --i)
	{
		x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
	}

	x.push_back(0.);
	return x;
	/*int m = A[0].size();
	std::vector<std::vector<double>> koef_progonki(m - 1);
	koef_progonki[m - 2].resize(2);
	koef_progonki[m - 2][0] = A[2][m - 1] / (-A[1][m - 1]);
	koef_progonki[m - 2][1] = -A[3][m - 1] / (-A[1][m - 1]);
	double temp = 0;
	for (int i = m - 3; i >= 0; --i)
	{
		koef_progonki[i].resize(2);
		temp = -A[1][i + 1] - A[0][i + 1] * koef_progonki[i + 1][0];
		koef_progonki[i][0] = A[2][i + 1] / temp;
		koef_progonki[i][1] = (-A[3][i + 1] + A[0][i + 1] * koef_progonki[i + 1][1]) / temp;
	}
	std::vector<double> X(m);
	X[0] = (-A[3][0] + A[0][0] * koef_progonki[0][1]) / (-A[1][0] - A[0][0] * koef_progonki[0][0]);
	for (int i = 1; i < m - 1; ++i)
		X[i] = koef_progonki[i - 1][0] * X[i - 1] + koef_progonki[i - 1][1];
	X[m - 1] = koef_progonki[m - 2][0] * X[m - 2] + koef_progonki[m - 2][1];*/
	/*return x;*/
}