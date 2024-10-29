#ifndef _TEXTOBJECT
#define _TEXTOBJECT
#include <iostream>
#include <fstream>


// abstract base class providing a read/write interface to all classes that inherit it
class TextObject {
private:
	// overloaded friend stream input/ouput functions 
	friend std::ostream & operator<< (std::ostream& os, const TextObject& r);
	friend std::istream & operator>> (std::istream& is, TextObject& r);
	friend std::ifstream & operator>> (std::ifstream& ifs, TextObject& r);
	friend std::ofstream & operator<< (std::ofstream& ofs, const TextObject& r);
public:
	// pure virtual methods for reading and writing, need to be overriden in derived classes
	virtual void write(std::ostream& os) const = 0;
	virtual void read(std::istream &is) = 0;
	virtual void readfile(std::ifstream& ifs) = 0;
	virtual void writefile(std::ofstream& ofs) const  = 0;
};


#endif