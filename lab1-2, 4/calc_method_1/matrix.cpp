#include "matrix.h"

extern const string connection_string_for_matrix_;
extern const string connection_string_for_res_;
template class matrix<double>;
template class matrix<float>;
extern int counter;
extern int iter;

template<typename T>
matrix<T>::matrix(base_type_of_matrix type, int size_)
{
	switch (type)
	{
	case base_type_of_matrix::singular:
		matrix::n = size_;
		A.resize(n);
		for (int i = 0; i < n; ++i)
			A[i].resize(n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (i == j) A[i][j] = 1; else A[i][j] = 0;
		break;
	case base_type_of_matrix::zerro:
		matrix::n = size_;
		A.resize(n);
		for (int i = 0; i < n; ++i)
			A[i].resize(n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				A[i][j] = 0;
		break;
	default:
		matrix::n = size_;
		break;
	}
}
template<typename T>
matrix<T>::matrix(type_of_matrix_from_file type, string connection_string_for_matrix_)
{

	switch (type)
	{
	case type_of_matrix_from_file::square:
	{
		std::ifstream file;
		file.open(connection_string_for_matrix_);
		file >> n;
		A.resize(n);
		for (int i = 0; i < n; ++i)
			A[i].resize(n + 1);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n + 1; ++j)
				file >> A[i][j];
		for (int i = 0; i < n; ++i)
			A[i].resize(A[i].size() - 1);
		file.close();
	}
	break;
	case type_of_matrix_from_file::n_slau:
	{
		std::ifstream file;
		file.open(connection_string_for_matrix_);
		file >> n;
		A.resize(n);
		for (int i = 0; i < n; ++i)
			A[i].resize(n + 1);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n + 1; ++j)
				file >> A[i][j];
		file.close();
	}
	break;
	default:
		matrix::n = 0;
		break;
	}
}


template<typename T>
vector<T> matrix<T>::get_decision_gaus()
{
	matrix B = *this;
	up_tringle();
// if (is_zerro_det()) throw invalid_argument("determinant is zero");
	vector<T> answ;
	answ.resize(n);
	T temp = 0;
	for (int i = n - 1; i >= 0; --i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			temp += A[i][j] * answ[j];
			counter++;
		}
		answ[i] = (A[i][n] - temp) / A[i][i];
		counter++;
		temp = 0;
	}
	ofstream write;
	write.open(connection_string_for_res_, std::ios::app);
	write << std::endl;
	for (int i = 0; i < n; ++i)
		write << setprecision(16) << answ[i] << std::endl;
	write.close();
	*this = B;
	return answ;
}
template<typename T>
vector<T> matrix<T>::get_decision_QR()
{
	matrix<T> Q = get_Q_R().first;
	matrix<T> R = get_Q_R().second;
	if (R.is_zerro_det()) throw invalid_argument("determinant is zero");

	vector<T> b(n);
	vector<T> b_zvezd(n);
	b = get_b();

	b_zvezd = Q.transponate() * b;
	for (int i = 0; i < n; ++i)
		R.A[i].push_back(b_zvezd[i]);

	vector<T> answ(n);
	T temp = 0;
	for (int i = n - 1; i >= 0; --i)
	{
		for (int j = i + 1; j < n; ++j)
			temp += R.A[i][j] * answ[j];
		answ[i] = (R.A[i][n] - temp) / R.A[i][i];
		temp = 0;
	}

	ofstream write;
	write.open(connection_string_for_res_);
	write << std::endl;
	for (int i = 0; i < n; ++i)
		write << setprecision(8) << answ[i] << endl;
	write.close();
	return answ;
}
template<typename T>
vector<T> matrix<T>::get_decision_simple_iteration()
{
	T tao = 0.003;
	int k = 0;
	double eps = 0.0000001;
	vec<T> x(n);
	vec<T> x0(n);
	vec<T> x1(n);
	vec<T> y(n);
	matrix<T> E(base_type_of_matrix::singular, n);
	matrix<T> C(base_type_of_matrix::singular, n);
	matrix<T> A = this->get_matrix();
	C = E - A;
	while (C.norm_kub() > 1)
	{
		tao /= 2;
		A = A * tao;
		C = E - A;
		A = this->get_matrix();
	}
	
	for (int i = 0; i < n; ++i)
		y[i] = get_b()[i] * tao;
	vec<T> xk(n);
	//for(int i = 0; i < 9249;++i)
		while (true)
	{
		++k;
		xk = vec<T>(C * x) + y;
		if (k == 1) x1 = xk;
	//	if (norm_kub(xk - x) < (1 - C.norm_kub()) / C.norm_kub() * eps) break;
		if (norm_kub(xk - x) <  eps) break;
	//	if (norm_kub(vec<T>(*this * xk) - vec<T>(get_b())) < eps) break;//cout<<"K_OSH= "<<k;
		x = xk;
	}
	cout << "tao = " << tao << endl;
	cout << "C.norm = " << C.norm_kub() << endl;
	cout << "k = " << k << endl;
	cout << "norm oshibki = " << norm_kub(x - vector<T> {5, -7, 12, 4}) << endl;
	cout << "k_est = " << ceil(log((1 - C.norm_kub()) * eps / norm_kub(x0 - x1)) / log(C.norm_kub()))<<endl;
	ofstream write;
	write.open(connection_string_for_res_, std::ios::app);
	write << std::endl;
	for (int i = 0; i < n; ++i)
		write << xk[i] << std::endl;
	write.close();
	return xk;
}
template<typename T>
vector<T> matrix<T>::get_decision_yakobi()
{
	T tao = 1;
	vec<T> x(n);
	vec<T> x0(n);
	vec<T> x1(n);
	vec<T> y(n);
	double eps = 0.0000001;
	matrix C(base_type_of_matrix::singular, n);
	C = get_matrix_for_yacobi();
	y = C.get_b();
	vec<T> xk(n);
	int k = 0;
	while (true)
	{
		++k;
		xk = vec<T>(C * x) + y;
		if (k == 1) x1 = xk;
	//	if (norm_kub(xk - x) < (1 - C.norm_kub()) / C.norm_kub() * eps) break;
		if (norm_kub(xk - x) <  eps) break;
	//	if (norm_kub(vec<T>(*this * xk) - vec<T>(get_b())) < eps) break;
		x = xk;
	}
	cout << "tao = " << tao << endl;
	cout << "C.norm = " << C.norm_kub() << endl;
	cout << "k = " << k << endl;
	cout << "norm oshibki = " << norm_kub(x - vector<T> {5, -7, 12, 4}) << endl;
	cout << "k_est = " << ceil(log((1 - C.norm_kub()) * eps / norm_kub(x0 - x1)) / log(C.norm_kub())) << endl;
	ofstream write;
	write.open(connection_string_for_res_, std::ios::app);
	write << std::endl;
	for (int i = 0; i < n; ++i)
		write << xk[i] << std::endl;
	write.close();
	return xk;

}
template<typename T>
vector<T> matrix<T>::get_decision_zeid()
{
	T omega = 1.1;
	vec<T> x(n);
	vec<T> x0(n);
	vec<T> x1(n);
	T sum1 = 0; T sum2 = 0;
	vec<T> xk(n);
	int k = 0;
	double eps = 0.0000001;
	matrix<T> C = (get_matrix_D() + get_matrix_L() * omega).inverse() * (get_matrix_D() * (1 - omega) - get_matrix_U() * omega);
	while (true)
	{
		++k;
		for (int i = 0; i < n; ++i)
		{
			sum1 = 0;
			for (int j = 0; j < i; ++j)
				sum1 += A[i][j] / A[i][i] * xk[j];
			sum2 = 0;
			for (int j = i+1; j < n; ++j)
				sum2 += A[i][j] / A[i][i] * x[j];
			xk[i] = -omega * sum1 + (1 - omega) * x[i] - omega * sum2 + omega * get_b()[i] / A[i][i];
		}
		if (k == 1) x1 = xk;
		if ((norm_kub(xk - x) < (1 - C.norm_kub()) / C.norm_kub() * eps)) break;
	//	if ((norm_kub(xk - x) <  eps)) break;
	//	if (norm_kub(vec<T>( * this * xk) - vec<T>(get_b())) < eps) break;
		x = xk;
	}

	//cout << "tao = " << tao << endl;
	cout << "C.norm = " << C.norm_kub() << endl;
	cout << "k = " << k << endl;
	cout << "norm oshibki = " << norm_kub(x - vector<T> {5, -7, 12, 4}) << endl;
	cout << "k_est = " << ceil(log((1 - C.norm_kub()) * eps / norm_kub(x0 - x1)) / log(C.norm_kub())) << endl;
	ofstream write;
	write.open(connection_string_for_res_, std::ios::app);
	write << std::endl;
	for (int i = 0; i < n; ++i)
		write << x[i] << std::endl;
	write.close();
	return x;
}
template<typename T>
vector<T> matrix<T>::get_decision_zeid_for_special_matrix(const vector<vector<T>>& G_)
{
	int n_ = 217;
	double eps = 0.0000001;
	vector<vector<T>> G = G_;
	vector<T> a = G[0];
	vector<T> b = G[1];
	vector<T> c = G[2];
	vector<T> d = G[3];
	T omega = 1;
	vec<T> x(n_);
	T sum1 = 0; T sum2 = 0;
	vec<T> xk(n_);
	int k = 0;
	while (true)
	{
		++k;
		xk[0] = d[0] / b[0] - c[0] / b[0] * x[1];
		for (int i = 1; i < n_ - 1; ++i)
		{
			xk[i] = d[i] / b[i] - a[i] / b[i] * xk[i - 1] - c[i] / b[i] * x[i + 1];
		}
		xk[n_ - 1] = d[n_ - 1] / b[n_ - 1] - a[n_ - 1] / b[n_ - 1] * xk[n_ - 2];
		if (fabs(norm_kub(xk - x)) < eps) break;
		x = xk;
	}
	cout << k << endl;
	ofstream write;
	write.open(connection_string_for_res_, std::ios::app);
	write << std::endl;
	for (int i = 0; i < n_; ++i)
		write << x[i] << std::endl;
	write.close();
	return x;
}
template<typename M>
vector<M>  get_decision_gaus_(matrix<M> B, vector<M> b)
{
	
	
		for (int i = 0; i < B.n; ++i)
		{
			B.A[i].resize(B.n);
			B.A[i].push_back(b[i]);
		}

	return B.get_decision_gaus();
}

template<typename T>
vector<T> matrix<T>::get_eigenvalues()
{
	T eps = 0.00001;
	//matrix<T> A_ = *this;
	matrix<T> A_ = (* this).get_hesinberg_matrix();
	matrix<T> E(base_type_of_matrix::singular,n);
	matrix<T> Q(base_type_of_matrix::singular,n);
	matrix<T> R(base_type_of_matrix::singular,n);
	vector<T> res(n);
	
	for (int k = 0; k < n; ++k)
	{
		

		while (true)
		//for (int m = 0; m < 200; ++m)
		{
			double omega = A_.A[n-k-1][n-k-1];
			A_ = A_ - E * omega;//+4 умножений
			Q = A_.get_Q_R().first;
			R = A_.get_Q_R().second;
			A_ = R * Q +E * omega;
			counter += 4;
			if (A_.norm_last_podstring() < eps) break;
			iter++;
		}
		res[n-k-1] = A_.A[n - k-1][n - k-1];
		A_ = A_.get_minor();
	}

    for(int i = 0;i < n; ++i)
	cout << setprecision(16) << res[i] << std::endl;
	return res;
}
template<typename T>
vector<vector<T>> matrix<T>::get_eigenvectors()
{
	T eps = 0.000001;
	matrix<T> A_ = *this;
	matrix E(base_type_of_matrix::singular, n);
	vector <T> values = get_eigenvalues();
	//vec<T> x(n);
	//x[0] = 1;
	vec<T> x_(  {0.,-0.7,0.6,-0.35});
	x_ = x_ * (1 / norm_defoalt(x_));
	vec<T> y(n);
	vector<vector<T>> res(n);
	for (int i = 0; i < n; ++i)
	{
		vec<T> x(n);
		x = x_;
		while (true)
		{
			y = get_decision_gaus_(A_ - E * values[i], static_cast<vector<T>>(x));
			x = y * (1/norm_defoalt(y));
			if (fabs(norm_kub(A_ * x) - norm_kub(x * values[i])) < eps) break;
		}
		res[i] = x;
	}
	cout << "*************" << endl;
	for (int i = 0; i < n; ++i)
		print(res[i]);
	return res;
}
template<typename T>
pair<vector<vector<T>>, vector<T>> matrix<T>::relay_algorithm()
{
	T eps = 0.00001;
	matrix<T> A_ = *this;
	matrix E(base_type_of_matrix::singular, n);
	vector<T> res_values(n);
	vector<vector<T>> res_vectors(n);
	vec<T> x_({ 0.,-0.7,0.6,-0.35 });
	x_ = x_ * (1 / norm_defoalt(x_));
	vec<T> x(n);
	x = x_;
	vec<T> y(n);
	//x[1] = 1.;
    for (size_t i = 0; i < n; ++i)
	{
		x = vec<T>({ { 0.,-0.7,0.6,-0.35 } });
		while (true)
		{
		    res_values[i] = dot(A_ * x, x);
			res_vectors[i] = static_cast<vector<T>>(x);
			y = get_decision_gaus_(A_-E*res_values[i], static_cast<vector<T>>(x));
			x = y * (1 / norm_defoalt(y));
			if (fabs(norm_kub(A_ * x) - norm_kub(x * res_values[i])) < eps) break;
		}
	}
	for (int i = 0; i < n; ++i)
		cout << setprecision(16) << res_values[i] << std::endl;
	cout << "*************" << endl;
	for (int i = 0; i < n; ++i)
		print(res_vectors[i]);
	return {res_vectors, res_values};
}

template<typename T>
matrix<T> matrix<T>::get_minor()
{
	matrix A_ = *this;
	A_.n = n - 1;
	A_.A.resize(A_.n);
	for (int i = 0; i < A_.n; ++i)
		A_.A[i].resize(A_.n);
	return A_;
}
template<typename T>
vector<T> matrix<T>::get_b()
{
	std::vector<T> b(n);
	for (int i = 0; i < n; ++i)
		b[i] = A[i][n];
	return b;
}
template<typename T>
void matrix<T>::set_b(const vector<T>& b)
{

	for (int i = 0; i < n; ++i)
	{
		A[i].resize(n + 1);
		A[i][n] = b[i];
	}
}
template<typename T>
matrix<T> matrix<T>::get_matrix()
{
	
	matrix <T> m = *this;
	if (A[1].size() > n) {
		for (int i = 0; i < n; ++i)
			m.A[i].resize(m.A[i].size() - 1);
		return m;
	}
	else return *this;
}
template<typename T>
void matrix<T>::up_tringle()
{

	for (int k = 0; k < n - 1; ++k)
	{
		part_choice(k);
		//	if (A[k][k] == 0) return (false);
		for (int i = k + 1; i < n; ++i)
		{
			T c = A[i][k] / A[k][k];
			counter++;
			for (int j = k; j < n + 1; j++)
			{
				A[i][j] = A[i][j] - c * A[k][j];
				counter++;
			}
		}
	}
	std::ofstream write;
	write.open(connection_string_for_res_);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n + 1; ++j)
			write << std::setprecision(6) << A[i][j] << " ";
		write << std::endl;
	}
	write.close();
	return;
}
template<typename T>
bool matrix<T>::is_zerro_det()
{
	T det = 1;
	for (int i = 0; i < n; ++i)
		det *= A[i][i];
	if (fabs(det) < 10e-16) return true;
	return false;
	//bool f = up_tringle();
	//double det = 1;
	//for (int i = 0; i < n; ++i)
	//	det *= A[i][i];
	//if ((!f) || (fabs(det) < 10e-16)) return true;

}
template<typename T>
void matrix<T>::part_choice(int num_of_current_line)
{
	int num_line_with_min = num_of_current_line;
	for (int i = num_of_current_line + 1; i < n; ++i)
	{
		if (fabs(A[num_of_current_line][num_of_current_line]) < fabs(A[i][num_of_current_line]))
			swap(A[num_of_current_line], A[i]);
	}
}

template<typename T>
matrix<T> matrix<T>::get_hesinberg_matrix()
{
	matrix A = *this;
	matrix T_(base_type_of_matrix::singular, n);
	double temp = 0;
	for(int j = 0; j < n-2; ++j)
		for (int i = j + 2; i < n; ++i)
		{
			int k = j + 1;
			int l = i;
			temp = sqrt(A.A[k][k - 1] * A.A[k][k - 1] + A.A[l][k - 1] * A.A[l][k - 1]);
			T_.A[k][k] = A.A[k][k - 1] / temp;
			T_.A[l][l] = A.A[k][k - 1] / temp;
			T_.A[k][l] = A.A[l][k - 1] / temp;
			T_.A[l][k] = -A.A[l][k - 1] / temp;
		//	A = T_ * A * T_.inverse();
			A = smart_dot_inverse(smart_dot(T_, A, k, l), T_.inverse(),k,l);
			counter += 5;
			T_ = matrix(base_type_of_matrix::singular, n);
		}
	return A;
}
template<typename T>
vector<vector<T>> matrix<T>::get_special_zeid_matrix()
{
	int n_ = 217;
	vector<T> under_dig(n_);
	vector<T> dig(n_);
	vector<T> up_dig(n_);
	vector<T> d(n_);
	d[0] = 1;
	for (int i = 0; i < n_-1; ++i)
	{
		under_dig[i+1] = 1;
		up_dig[i] = 6;
		dig[i] = 8;
		d[i+1] = i+1;
	}
	under_dig[n_- 1] = 1;
	dig[n_ - 1] = 8;
	d[n_-1] = n;
	return {under_dig, dig, up_dig, d};
}
template<typename T>
pair<matrix<T>, matrix<T>> matrix<T>::get_Q_R()
{
	matrix<T> B(base_type_of_matrix::singular, n);
	matrix<T> R = *this;
	matrix<T> temp(base_type_of_matrix::singular, n);

	for (int i = 0; i <= n - 2; ++i)
	{
	//	int j = i + 1;
		for (int j = i + 1; j < n; ++j)
		{

			temp.A[i][i] = temp.A[j][j] = (R.A[i][i]) / (sqrt(R.A[i][i] * R.A[i][i] + R.A[j][i] * R.A[j][i]));
			temp.A[i][j] = (R.A[j][i]) / (sqrt(R.A[i][i] * R.A[i][i] + R.A[j][i] * R.A[j][i]));
			temp.A[j][i] = -temp.A[i][j];
			R = smart_dot_for_R(temp, R, i, j);
			B = smart_dot(temp, B, i, j);
			temp = matrix<T>(base_type_of_matrix::singular, n);
			counter += 5;
		}
	}
	for (int i = 0; i < n; ++i)
		R.A[i].resize(n);

	std::ofstream write;
	write.open(connection_string_for_res_);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n ; ++j)
			write << std::setprecision(16) << R.A[i][j] << " ";
		write << std::endl;
	}
	write << std::endl;
	write << std::endl;
	write << std::endl;
	matrix<T> BB = B.transponate();
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			write << std::setprecision(16) << BB.A[i][j] << " ";
		write << std::endl;
	}
	write.close();

	return { BB, R };
}
template<typename T>
T matrix<T>::get_obus(string type_of_norm)
{
	matrix<T> A_(base_type_of_matrix::zerro, n);
	matrix<T> E(base_type_of_matrix::singular, n);


	matrix<T> Q = get_Q_R().first;
	matrix<T> R = get_Q_R().second;

	vector<T> b_zvezd;
	vector<T> x;

	E.transponate();
	for (int i = 0; i < n; ++i)
	{
		b_zvezd = Q.transponate() * E.A[i];
		x = get_decision_gaus_(R, b_zvezd);
		//	for (int j = 0; j < n; ++j)
		A_.A[i] = x;
		Q.transponate();
	}
	A_.transponate();
	if (type_of_norm == "okt")
		return norm_okt() * A_.norm_okt();
	else return norm_kub() * A_.norm_kub();
}
template<typename T>
matrix<T> matrix<T>::transponate()
{
	for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
			std::swap(A[i][j], A[j][i]);
	counter += 16;
	return *this;
}
template<typename T>
matrix<T> matrix<T>::inverse()
{
	matrix<T> res(base_type_of_matrix::singular, n);
	matrix<T> E(base_type_of_matrix::singular, n);
	vector<T> b_zvezd(n);
	for (int i = 0; i < n; ++i)
	{
	//	b_zvezd = *this * E.A[i];
		res.A[i] = get_decision_gaus_(this->get_matrix(), E.A[i]);
	}
	res.transponate();
	return res;
}
template<typename M>
matrix<M> smart_dot(const matrix<M>& A, const matrix<M>& B, int k, int m)
{
	matrix<M> res = B;
	for (int i = 0; i < A.n; ++i)
	{
		counter += 4;
		res.A[k][i] = A.A[k][k] * B.A[k][i] + A.A[k][m] * B.A[m][i];
		res.A[m][i] = A.A[m][k] * B.A[k][i] + A.A[k][k] * B.A[m][i];
	}
	return res;
}
template<typename M>
matrix<M> smart_dot_inverse(const matrix<M>& A, const matrix<M>& B, int k, int m)
{
	matrix<M> res = A;
	for (int i = 0; i < A.n; ++i)
	{
		counter += 4;
		res.A[i][k] = A.A[i][k] * B.A[k][k] + A.A[i][m] * B.A[m][k];
		res.A[i][m] = A.A[i][k] * B.A[k][m] + A.A[i][m] * B.A[m][m];
	}
	return res;
}
template<typename M>
matrix<M> smart_dot_for_R(const matrix<M>& A, const matrix<M>& B, int k, int m)
{
	matrix<M> res = B;
	for (int i = k; i < A.n; ++i)
	{
		counter += 4;
		res.A[k][i] = A.A[k][k] * B.A[k][i] + A.A[k][m] * B.A[m][i];
		res.A[m][i] = A.A[m][k] * B.A[k][i] + A.A[k][k] * B.A[m][i];
	}
	return res;
}
template<typename T>
matrix<T>& matrix<T>::operator=(matrix<T> B)
{
	n = B.n;
	A = B.A;
	return *this;
}
template<typename T>
matrix<T> matrix<T>::operator*(const matrix<T>& B)
{
	matrix<T> C(base_type_of_matrix::zerro, n);
	if (n != B.n) throw invalid_argument("incorrect size");
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				C.A[i][j] += A[i][k] * B.A[k][j];
	return C;
}
template<typename T>
vector<T> matrix<T>::operator*(vector<T> b)
{
	vector<T> res(n);
	if (n != b.size()) throw invalid_argument("incorrect size");
	for (int i = 0; i < n; ++i)
		for (int k = 0; k < n; ++k)
			res[i] += A[i][k] * b[k];
	return res;
}
template<typename T>
matrix<T> matrix<T>::operator+(const matrix<T>& B)
{
	matrix<T> temp = *this;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			temp.A[i][j] += B.A[i][j];
	return temp;
}
template<typename T>
matrix<T> matrix<T>::operator-(const matrix<T>& B)
{
	matrix<T> temp = *this;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			temp.A[i][j] -= B.A[i][j];
	return temp;
}
template<typename T>
matrix<T> matrix<T>::operator*(double a)
{
	matrix<T> temp = *this;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			temp.A[i][j] *= a;
	return temp;
}

template<typename T>
T matrix<T>::dot(vector<T> vec1, vector<T> vec2)
{
	T res = 0;
	for (int i = 0; i < vec1.size(); ++i)
		res += vec1[i] * vec2[i];
	return res;
}

template<typename T>
T matrix<T>::norm_kub()
{
	T max = 0;
	T sum = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			sum += fabs(A[i][j]);
		if (max < sum) max = sum;
		sum = 0;
	}
	return max;
}
template<typename T>
T matrix<T>::norm_okt()
{
	T max = 0;
	T sum = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			sum += fabs(A[j][i]);
		if (max < sum) max = sum;
		sum = 0;
	}
	return max;
}
template<typename T>
T matrix<T>::norm_okt(vector<T> vec)
{
	T res = 0;
	for (int i = 0; i < n; ++i)
		res += fabs(vec[i]);
	return res;
}
template<typename T>
T matrix<T>::norm_kub(vector<T> vec)
{
	T res = 0;
	for (int i = 0; i < n; ++i)
		if (res < fabs(vec[i])) res = fabs(vec[i]);
	return res;
}
template<typename T>
T matrix<T>::norm_defoalt(vector<T> vec)
{
	T res = 0;
	for (int i = 0; i < n; ++i)
		res+=vec[i]*vec[i];
	return sqrt(res);
}
template<typename T>
T matrix<T>::norm_last_podstring()
{
	T temp = 0;
	for (int i = 0; i < n - 1; ++i)
		temp += A[n - 1][i]* A[n - 1][i];
	return sqrt(temp);
}

template<typename T>
T matrix<T>::get_nevyaz()
{
	vector<T> b = get_b();
	vector<T> x = get_decision_QR();
	vector<T> b_zvezd = *this * x;
	vector<T> nevyaz(n);
	T res = 0.;
	for (int i = 0; i < n; ++i)
	{
		nevyaz[i] = b[i] - b_zvezd[i];
		res += nevyaz[i] * nevyaz[i];
	}
	res = sqrt(res);
	return res;
}
template<typename T>
vector<T> matrix<T>::get_outraged_decision()
{
	vector<T> b = get_b();
	for (int i = 0; i < n; ++i)
		b[i] += 0.01;
	matrix<T> B = *this;
	for (int i = 0; i < n; ++i)
		B.A[i].resize(n);
	return(get_decision_gaus_(B, b));
}
template<typename T>
T matrix<T>::get_down_est(string type_of_norm)
{
	vector<T> delt_b(n);
	for (int i = 0; i < n; ++i)
		delt_b[i] = 0.01;
	vector<T> b = get_b();
	vector<T> x = get_decision_gaus();
	matrix B = *this;
	for (int i = 0; i < n; ++i)
		B.A[i].resize(n);
	vector<T> delt_x = get_decision_gaus_(B, delt_b);
	if (type_of_norm == "okt") return (norm_okt(delt_x) / norm_okt(x)) / (norm_okt(delt_b) / norm_okt(b));
	return (norm_kub(delt_x) / norm_kub(x)) / (norm_kub(delt_b) / norm_kub(b));
}
template<typename T>
matrix<T> matrix<T>::get_matrix_for_yacobi()
{
	matrix<T> C(base_type_of_matrix::singular,n);
	vector<T> b(n);
	for (int i = 0; i < n; ++i)
	{
		b[i] = get_b()[i] / A[i][i];
		for (int j = 0; j < n; ++j)
			if (i != j)
				C.A[i][j] = -A[i][j] / A[i][i];
			else C.A[i][j] = 0;
	}
	C.set_b(b);
	return C;
}

template<typename T>
matrix<T> matrix<T>::get_matrix_L()
{
	matrix<T> L(base_type_of_matrix::zerro, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j)
			L.A[i][j] = A[i][j];
	return L;
}
template<typename T>
matrix<T> matrix<T>::get_matrix_D()
{
	matrix<T> D(base_type_of_matrix::zerro, n);
	for (int i = 0; i < n; ++i)
		D.A[i][i] = A[i][i];
	return D;
}
template<typename T>
matrix<T> matrix<T>::get_matrix_U()
{
	matrix<T> U(base_type_of_matrix::zerro, n);
	for (int i = 0; i < n; ++i)
		for (int j = i+1; j < n; ++j)
			U.A[i][j] = A[i][j];
	return U;
}


template<typename T>
void matrix<T>::print(const vector<T>& v)
{
	for (int i = 0; i < v.size(); ++i)
		cout << v[i] << endl;
	cout << "*************" << endl;
}



template<typename T>
vec<T>::vec(int n_)
{
	n = n_;
	b.resize(n_);
	for (int i = 0; i < n; ++i)
		b[i] = 0;
}
template<typename T>
vec<T>::vec(const vector<T>& v)
{
	n = v.size();
	b = v;
}
template<typename T>
vec<T> vec<T>::operator+(vec<T> c)
{
	vec<T> temp = *this;
	for (int i = 0; i < c.b.size(); ++i)
		temp.b[i] += c.b[i];
	return temp;
}
template<typename T>
vec<T> vec<T>::operator-(vec<T> c)
{
	vec<T> temp = *this;
	for (int i = 0; i < c.b.size(); ++i)
		temp.b[i] -= c.b[i];
	return temp;
}
template<typename T>
vec<T> vec<T>::operator*(T c)
{
	vec<T> temp = *this;
	for (int i = 0; i < b.size(); ++i)
		temp.b[i] *= c;
	return temp;
}

template<typename T>
void vec<T>::push_back(T value)
{
	b.push_back(value);
	n = b.size();
}

template<typename T>
T& vec<T>::operator[](int i)
{
	return b[i];
}
template<typename T>
T vec<T>::operator[](int i) const
{
	return b[i];
}
template<typename T>
vec<T>::operator vector<T>()
{
	vector<T> res(n);
	res = b;
	return res;
}