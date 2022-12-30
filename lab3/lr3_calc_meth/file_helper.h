#pragma once
#include <vector>
#include <fstream>
class file_helper
{
private:
	std::string connection_string;
	std::ofstream file;
public:
	file_helper(std::string connection_string_="C:/Users/Артем/Desktop/polynom points5.txt");
	void write(std::pair<double, double> point);
	~file_helper();
};

