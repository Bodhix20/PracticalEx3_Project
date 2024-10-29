#include "textobject.h"

// implementations of the friend stream input/output overloaded functions
// In each case they call the appropriate virtual read/write method to do the work
std::ofstream & operator<< (std::ofstream &ofs, const TextObject& r)   // write to a file
{
	r.writefile(ofs);
	return ofs;
}


std::ifstream & operator>> (std::ifstream &ifs, TextObject& r)  // read from a file
{
	r.readfile(ifs);
	return ifs;
}

std::ostream & operator<< (std::ostream &os, const TextObject& r)    // write to the screen
{
	r.write(os);
	return os;
}

std::istream & operator>> (std::istream &is, TextObject& r)   // read from the keyboard
{
	r.read(is);
	return is;
}


