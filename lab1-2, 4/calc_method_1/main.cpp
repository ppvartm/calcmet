#include "matrix.h"
#include <iostream>
extern const string connection_string_for_res_ = "C:/Users/Артем/Desktop/res.txt";
extern const string connection_string_for_matrix_ = "C:/Users/Артем/Desktop/тест4.txt";

int counter;
int iter;

int main()
{
	//ura
	matrix<double> A(type_of_matrix_from_file::n_slau);
	//cout << A.norm_kub() << endl;
	//cout << A.norm_okt() << endl;
//	A.inverse();
//	matrix<double> E(base_type_of_matrix::singular, 4);
//	cout<<E.get_obus("kub");
	//	vector<double> b = A.get_outraged_decision();
//cout << A.get_nevyaz() << endl;
	//int i = 5;
	//i = i++;
	//cout << i << endl;
//	A.relay_algorithm();
//A.get_eigenvectors();
//cout << "counter= " << counter << endl;
//cout << "iter= " << iter << endl;
//A = A.transponate();
//A.get_eigenvectors();
	//A.get_decision_simple_iteration();
//	A.get_decision_yakobi();
	//A.get_decision_zeid();
	A.get_decision_zeid_for_special_matrix(A.get_special_zeid_matrix());


//	cout << A.get_obus("kub") << endl;
//	cout << A.get_down_est("kub") << endl;
//	matrix B(type_of_matrix_from_file::square, "C:/Users/Артем/Desktop/B.txt");
//	matrix C(type_of_matrix_from_file::square, "C:/Users/Артем/Desktop/C.txt");
//	smart_dot(B, C, 0, 2);

}