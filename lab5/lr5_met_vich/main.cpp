#include<iostream>
#include<fstream>
#include<vector>
#include "foo_for_ex.h"
#include "eq_solver.h"

std::ofstream file1;
std::ofstream file2;
std::ofstream file3;
std::ofstream file4;
std::ofstream file5;
std::ofstream file6;
std::ofstream file7;
std::ofstream file8;
std::ofstream file9;
std::ofstream file10;

void foo(double k,double x, double y);


int main()
{

	/*
	file1.open("C:/Users/Артем/Desktop/files/1_4.txt");
	file2.open("C:/Users/Артем/Desktop/files/5_6.txt");
	file3.open("C:/Users/Артем/Desktop/files/7_8.txt");
	file4.open("C:/Users/Артем/Desktop/files/9_10.txt");
	file5.open("C:/Users/Артем/Desktop/files/11_12.txt");
	file6.open("C:/Users/Артем/Desktop/files/13_15.txt");
	file7.open("C:/Users/Артем/Desktop/files/16_20.txt");
	file8.open("C:/Users/Артем/Desktop/files/21_25.txt");
	file9.open("C:/Users/Артем/Desktop/files/26_30.txt");
	file10.open("C:/Users/Артем/Desktop/files/30.txt");

	double i = -10, j = -10;
	double h = 0.2;
	std::vector<double> answ;
	while (i < 10)
	{
		while (j < 10) {
		 answ = solver.newton_matrix_method(system1, dsistem1_err, {i, j});
		 foo(answ[2], i, j);
			j += h;
		}
		j = -10;
		i += h;

	}
	std::cout << answ[0] << "   " << answ[1] << std::endl;
	*/

	eq_solver solver; // методы этой хрени решают уравнения 
	double y = solver.method_sek(foo1, { 0,1 });
    double x = solver.newton_method(foo1, dfoo1 , { 0,1 });

	
}



void foo(double k, double x, double y)
{
	if (k < 5) {
		file1 << x << "  " << y << std::endl; return;
	}
	if (k < 7) {
		file2 << x << "  " << y << std::endl; return;
	}
	if (k < 9) {
		file3 << x << "  " << y << std::endl; return;
	}
	if (k < 11) {
		file4 << x << "  " << y << std::endl; return;
	}
	if (k < 13) {
		file5 << x << "  " << y << std::endl; return;
	}
	if (k < 16) {
		file6 << x << "  " << y << std::endl; return;
	}
	if (k < 21) {
		file7 << x << "  " << y << std::endl; return;
	}
	if (k < 26) {
		file8 << x << "  " << y << std::endl; return;
	}
	if (k < 31) {
		file9 << x << "  " << y << std::endl; return;
	}
	if (k > 30) {
		file10 << x << "  " << y << std::endl; return;
	}
	
}