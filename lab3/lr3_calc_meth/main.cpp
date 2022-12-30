#include <iostream>
#include "polynom.h"
#include "functions.h"

std::vector<double> create_regular_nodes(std::pair<double, double> left_right, int n);
std::vector<double> create_Chebish_nodes(std::pair<double, double> left_right, int n);


int main()
{
   /* std::vector<double> diap1 = create_regular_nodes({ -2,2 }, 3); 
	std::vector<double> diap2 = create_regular_nodes({ -1,1 }, 10);
	std::vector<double> diap3 = create_regular_nodes({ -1,1 }, 100);
	std::vector<double> diap4 = create_Chebish_nodes({ -1,1 }, 3);
	std::vector<double> diap5 = create_Chebish_nodes({ -1,1 }, 10);
	std::vector<double> diap6 = create_Chebish_nodes({ -1,1 }, 50);
	std::vector<double> diap7 = create_regular_nodes({ -2,2 }, 300);*/
	polynom p(10);
	//p.lagrange_intrepol(diap1, &foo2,diap2);
	//p.spline_interpol(diap1, &foo3, diap2);
	std::vector<double> diap3 = create_regular_nodes({ -3,3 }, 8);
	std::vector<double> diap4 = create_regular_nodes({ -3,3 }, 16);
	std::vector<double> diap5 = create_regular_nodes({ -3,3 }, 32);
	std::vector<double> diap6 = create_regular_nodes({ -3,3 }, 64);
	std::vector<double> diap7 = create_regular_nodes({ -3,3 }, 128);
	std::vector<double> diap8 = create_regular_nodes({ -3,3 }, 256);
	std::vector<double> diap9 = create_Chebish_nodes({ -2,2 }, 200);
	std::vector<double> diap10 = create_regular_nodes({ -10,10 }, 2000);
	std::vector<double> diap11 = create_Chebish_nodes({ -2,2 }, 100);
	std::cout <<"lg8:   " << std::setprecision(13) << foo3(1.23452)- p.lagrange_interpol(diap3, &foo3, 1.23452) <<  std::endl;
	std::cout << "sp8:   "<< std::setprecision(13) << foo3(1.23452)- p.spline_interpol(diap3, &foo3, 1.23452) << std::endl;
	std::cout << "lg16:   " << std::setprecision(13) << foo3(1.23452)- p.lagrange_interpol(diap4, &foo3, 1.23452)  << std::endl;
	std::cout << "sp16:   " << std::setprecision(13) << foo3(1.23452)-p.spline_interpol(diap4, &foo3, 1.23452)  << std::endl;
	std::cout << "lg32:   " << std::setprecision(13) << foo3(1.23452)-p.lagrange_interpol(diap5, &foo3, 1.23452)  << std::endl;
	std::cout << "sp32:   " << std::setprecision(13) << foo3(1.23452)-p.spline_interpol(diap5, &foo3, 1.23452)  << std::endl;
	std::cout << "lg64:   " << std::setprecision(13) << foo3(1.23452) - p.lagrange_interpol(diap6, &foo3, 1.23452) << std::endl;
	std::cout << "sp64:   " << std::setprecision(13) << foo3(1.23452) - p.spline_interpol(diap6, &foo3, 1.23452) << std::endl;
	std::cout << "lg128:   " << std::setprecision(13) << foo3(1.23452) - p.lagrange_interpol(diap7, &foo3, 1.23452) << std::endl;
	std::cout << "sp128:   " << std::setprecision(13) << foo3(1.23452) - p.spline_interpol(diap7, &foo3, 1.23452) << std::endl;
	std::cout << "lg256:   " << std::setprecision(13) << foo3(1.23452) - p.lagrange_interpol(diap8, &foo3, 1.23452) << std::endl;
	std::cout << "sp256:   " << std::setprecision(13) << foo3(1.23452) - p.spline_interpol(diap8, &foo3, 1.23452) << std::endl;
	
	//p.progonka({ {0,0,3,-5},{0,1,1,-6},{0,-4,1,0},{0,2,5,-9}});
	//matrix<double> M(base_type_of_matrix::singular, 1);
	//M.get_decision_zeid_for_special_matrix({ {0,0,3,-5},{0,1,1,-6},{0,-4,1,0},{0,2,5,-9} });
//	p.spline_interpol(diap1, &foo1, diap7);
//std::cout << " lagrange, regular, n = 3           " << p.get_prec(diap1, &foo1, "lagrange")<<std::endl;
//std::cout << " lagrange, regular, n = 10         " << p.get_prec(diap2, &foo1, "lagrange") << std::endl;
//std::cout << " lagrange, regular, n = 50         " << p.get_prec(diap3, &foo1, "lagrange") << std::endl;
//std::cout << " lagrange, Chebish, n = 3           " << p.get_prec(diap4, &foo1, "lagrange") << std::endl;
//std::cout << " lagrange, Chebish, n = 10         " << p.get_prec(diap5, &foo1, "lagrange") << std::endl;
//std::cout << " lagrange, Chebish, n = 50         " << p.get_prec(diap6, &foo1, "lagrange") << std::endl;
}

std::vector<double> create_regular_nodes(std::pair<double, double> left_right, int n)
{
	std::vector<double> answ(n + 1);
	for (size_t i = 0; i < answ.size(); ++i)
		answ[i] = left_right.first + i * (left_right.second - left_right.first) / n;
	
	return answ;
}
std::vector<double> create_Chebish_nodes(std::pair<double, double> left_right, int n)
{
	std::vector<double> answ(n + 1);
	for (size_t i = 0; i < answ.size(); ++i)
		answ[i] = (left_right.first + left_right.second) / 2 + (left_right.second - left_right.first) / 2
		* cos(((2 * i + 1) * atan(1) * 4) / (2 * (n + 1)));
	return answ;
}
