#include "file_helper.h"


file_helper::file_helper(std::string connection_string_) :connection_string(connection_string_)
{
	file.open(connection_string);
}

void file_helper::write(std::pair<double, double> point)
{
		file <<point.first << "  " << point.second << std::endl;
}
file_helper::~file_helper()
{
	file.close();
}