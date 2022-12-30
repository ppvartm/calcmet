#include "file_helper.h"
file_helper::file_helper(std::string connection_string_) :connection_string(connection_string_) 
{
	file.open(connection_string);
};
void file_helper::write(double x, double y)
{

}

file_helper::~file_helper()
{
	file.close();
}