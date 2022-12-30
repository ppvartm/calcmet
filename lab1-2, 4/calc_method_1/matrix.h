#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>


using namespace std;
enum class base_type_of_matrix { singular, zerro };
enum class type_of_matrix_from_file { square, n_slau };
extern const string connection_string_for_matrix_;
template<typename T>
class vec;

template<typename T>
class matrix
{
private:
	int n;
	vector<vector<T>> A;
public:
	matrix(base_type_of_matrix, int n_);
	matrix(type_of_matrix_from_file, string connection_string_for_matrix_ = connection_string_for_matrix_);

	vector<T> get_decision_gaus();
	vector<T> get_decision_QR();
	vector<T> get_decision_simple_iteration();
	vector<T> get_decision_yakobi();
	vector<T> get_decision_zeid();
	vector<T> get_decision_zeid_for_special_matrix(const vector<vector<T>>& G_);

	template<typename M>
	vector<M> friend get_decision_gaus_(matrix<M> B, vector<M> b);


	vector<T> get_eigenvalues();
	vector<vector<T>> get_eigenvectors();
	pair<vector<vector<T>>, vector<T>> relay_algorithm();
		
	vector<T> get_b();
	void set_b(const vector<T>& b);
	matrix<T> get_matrix();
	bool is_zerro_det();
	void up_tringle();
	void part_choice(int num_of_current_line);



	matrix<T> get_minor();
	matrix<T> get_hesinberg_matrix();
	vector<vector<T>> get_special_zeid_matrix();
	pair<matrix<T>, matrix<T>> get_Q_R();
	T get_obus(string type_of_norm);
	matrix<T> transponate();
	matrix<T> inverse();
	template<typename M>
	friend matrix<M> smart_dot(const matrix<M>& A, const matrix<M>& B, int k, int m);
	template<typename M>
	friend matrix<M> smart_dot_inverse(const matrix<M>& A, const matrix<M>& B, int k, int m);
	template<typename M>
	friend matrix<M> smart_dot_for_R(const matrix<M>& A, const matrix<M>& B, int k, int m);
	matrix<T>& operator=(matrix<T> B);
	matrix<T> operator*(const matrix<T>& B);
	vector<T> operator*(vector<T> b);
	matrix<T> operator+(const matrix<T>& B);
	matrix<T> operator-(const matrix<T>& B);
	matrix<T> operator*(double a);

	T dot(vector<T> vec1, vector<T> vec2);

	T norm_okt();
	T norm_kub();
	T norm_okt(vector<T> vec);
	T norm_kub(vector<T> vec);
	T norm_defoalt(vector<T> vec);
	T norm_last_podstring();

	T get_nevyaz();
	vector<T> get_outraged_decision();
	T get_down_est(string type_of_norm);
	matrix<T> get_matrix_for_yacobi();
	matrix<T> get_matrix_L();
	matrix<T> get_matrix_D();
	matrix<T> get_matrix_U();



	void print(const vector<T>& v);
};

template<typename T>
class vec
{

private:
	int n;
	vector<T> b;

public:
	vec(int n_);
	vec(const vector<T>& v);

	void push_back(T value);
	vec<T> operator+(vec<T> c);
	vec<T> operator-(vec<T> c);
	vec<T> operator*(T c);
	T& operator[](int i);
	T operator[](int i) const;

	operator vector<T>();
};