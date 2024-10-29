#ifndef VECTOR_H //Include guard
#define VECTOR_H

#include <iostream> //Generic IO operations
#include <fstream>  //File IO operations
#include <stdexcept> //provides exceptions
#include <vector>  // std vector upon which our Vector is based


/**
*  A vector class for data storage of a 1D array of doubles
*  \n The implementation is derived from the standard container vector std::vector
*  \n We use private inheritance to base our vector upon the library version whilst
*  \nallowing usto expose only those base class functions we wish to use - in this
*  \ncase the array access operator []
*
* The Vector class provides:
* \n-basic constructors for creating vector obcjet from other vector object,
or by creating empty vector of a given size,
* \n-input and oput operation via >> and << operators using keyboard or file
* \n-basic operations like access via [] operator, assignment and comparision
*/
template<typename T>
class Vector : private std::vector<T> {
	typedef std::vector<T> vec;
public:
	using vec::operator[];  // elevates the array access operator inherited from std::vector
							// to public access in Vector
	// CONSTRUCTORS
	/**
	* Default constructor.  Intialize an empty Vector object
	* @see Vector(int Num)
	* @see Vector(const Vector& v)
	*/
    Vector(); // default constructor


	/**
	* Explicit alterative constructor takes an intiger.
	* it is explicit since implicit type conversion int -> vector doesn't make sense
	* Intialize Vector object of size Num
	* @see Vector()
	* @see Vector(const Vector& v)
	* @exception invalid_argument ("vector size negative")
	*/
    explicit Vector(int Num /**< int. Size of a vector */);

	/**
	* Copy constructor takes an Vector object reference.
	* Intialize Vector object with another Vector object
	* @see Vector()
	* @see Vector(int Num)
	*/
    Vector(const Vector<T>& v); 


	// OVERLOADED OPERATORS
	/**
	* Overloaded assignment operator
	* @see operator==(const Vector& v)const
	* @param v Vector to assign from
	* @return the object on the left of the assignment
	*/
    Vector<T>& operator=(const Vector<T>& v  /**< Vecto&. Vector to assign from */);

	
	// ACCESSOR METHODS
	/** Normal get method that returns integer, the size of the vector
	* @return int. the size of the vector
	*/
	int getSize() const;  


	// KEYBOARD/SCREEN INPUT AND OUTPUT
	/**
	* Overloaded istream >> operator. Keyboard input
	* if vector has size user will be asked to input only vector values
	* if vector was not initialized user can choose vector size and input it values
	* @see operator>>(std::ifstream& ifs, Vector& v)
	* @see operator<<(std::ostream& os, const Vector& v)
	* @see operator<<(std::ofstream& ofs, const Vector& v)
	* @return std::istream&. the input stream object is
	* @exception std::invalid_argument ("read error - negative vector size");
	*/
    template <typename T> friend std::istream& operator >> (std::istream& is, /**< keyboard input straem. For user input */
		Vector<T>& v /**< Vector&. vector to write to */
		);


	/**
	* Overloaded ifstream << operator. Display output.
	* @see operator>>(std::istream& is, Vector& v)
	* @see operator>>(std::ifstream& ifs, Vector& v)
	* @see operator<<(std::ofstream& ofs, const Vector& v)
	* @return std::ostream&. the output stream object os
	*/
	template <typename T> friend std::ostream& operator<<(std::ostream& os, /**< output file stream */
		const Vector<T>& v /**< vector to read from */
		);

	/**
	* Overloaded ifstream >> operator. File input
	* the file output operator is compatible with file input operator,
	* ie. everything written can be read later.
	* @see operator>>(std::istream& is, Vector& v)
	* @see operator<<(std::ostream& os, const Vector& v)
	* @see operator<<(std::ofstream& ofs, const Vector& v)
	* @return ifstream&. the input ifstream object ifs
	* @exception std::invalid_argument ("file read error - negative vector size");
	*/
	template <typename T> friend std::ifstream& operator >> (std::ifstream& ifs, /**< input file straem. With opened matrix file */
		Vector<T>& v /**< Vector&. vector to write to */
		);


	/**
	* Overloaded ofstream << operator. File output.
	* the file output operator is compatible with file input operator,
	* ie. everything written can be read later.
	* @see operator>>(std::istream& is, Vector& v)
	* @see operator>>(std::ifstream& ifs, Vector& v)
	* @see operator<<(std::ostream& os, const Vector& v)
	* @return std::ofstream&. the output ofstream object ofs
	*/
	 template <typename T> friend std::ofstream& operator<<(std::ofstream& ofs, /**< outputfile stream. With opened file */
		const Vector<T>& v /**< Vector&. vector to read from */
		);
};


// CONSTRUCTORS
/*=
* Default constructor (empty vector)
*/
template <typename T>
Vector<T>::Vector() : std::vector<T>() {}


/*
* Alternate constructor - creates a vector of a given size
*/
template <typename T>
Vector<T>::Vector(int Num) : std::vector<T>()
{
	// set the size
	(*this).resize(Num);

	// initialise with zero
	std::size_t i;
	for (i = 0; i < (*this).size(); i++) (*this)[i] = 0;
}

/*
* Copy constructor
*/
template <typename T>
Vector<T>::Vector(const Vector& copy) : std::vector<T>()
{
	(*this).resize(copy.size());
    // copy the data members (if vector is empty then num==0)
	std::size_t i;
    for (i=0; i<copy.size(); i++) (*this)[i]=copy[i]; 
}

/*
* accessor method - get the size
*/
template <typename T>
int Vector<T>::getSize() const
{
	return (*this).size();
}

// OVERLOADED OPERATORS
/*
* Operator= - assignment
*/
template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& copy)
{
	(*this).resize(copy.size());
	std::size_t i;
    for (i=0; i<copy.size(); i++) (*this)[i] = copy[i]; 
    return *this;
}



// INPUT AND OUTPUT
/*
* keyboard input , user friendly
*/
template <typename T>
std::istream& operator>>(std::istream& is, Vector<T>& v)
{
	if (!v.size()) {
		int n;

		std::cout << "input the size for the vector" << std::endl;
		is >> n;
		//check input sanity
		if(n < 0) throw std::invalid_argument("read error - negative vector size");

		// prepare the vector to hold n elements
		v = Vector<T>(n);
	}
	// input the elements
	std::cout << "input "<< v.size() <<" vector elements" << std::endl;
	std::size_t i;
	for (i=0; i<v.size(); i++) is >> v[i];

    // return the stream object
    return is;
}

/*
* file input - raw data, compatible with file writing operator
*/
template <typename T>
std::ifstream& operator>>(std::ifstream& ifs, Vector<T>& v) 
{
    int n;

	if (!v.getSize()) {

		// read size from the file
		ifs >> n;
		//check input sanity
		if(n < 0) throw std::invalid_argument("file read error - negative vector size");

		// prepare the vector to hold n elements
		v = Vector<T>(n);
	}
    // input the elements
    for (int i=0; i<v.getSize(); i++) ifs >> v[i];

    // return the stream object
    return ifs;
}

/*
* screen output, user friendly
*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
{
    if (v.size() > 0) {
		std::size_t i;
        for (i=0; i<v.size(); i++) os << v[i]  << " ";
        os << std::endl;
    }
    else
    {
        os << "Vector is empty." << std::endl;
    }
    return os;
}

/*
* file output - raw data, comaptible with file reading operator
*/
template <typename T>
std::ofstream& operator<<(std::ofstream& ofs, const Vector<T>& v)
{
    //put vector size in first line (even if it is zero)
    ofs << v.size() << std::endl;
    //put data in second line (if size==zero nothing will be put)
	std::size_t i;
    for (i=0; i<v.size(); i++) ofs << v[i]  <<  "\n";
    ofs << std::endl;
    return ofs;
}



#endif