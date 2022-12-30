#pragma once
#include<fstream>
#include<string>
class file_helper
{
private:
	std::ofstream file;
	std::string connection_string;
public:
	file_helper(std::string connection_string_);
	void write(double x, double y);

	~file_helper();

};

