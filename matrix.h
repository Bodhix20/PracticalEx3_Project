#ifndef MATRIX_H  //include guard
#define MATRIX_H

#include <iostream>  //generic IO
#include <fstream>   //file IO
#include <stdexcept> //provides exceptions
#include "vector.h"  //we use Vector in Matrix code


/**
*  A matrix class for data storage of a 2D array of doubles
*  \n The implementation is derived from the standard container vector std::vector
*  \n We use private inheritance to base our vector upon the library version whilst
*  \nallowing usto expose only those base class functions we wish to use - in this
*  \ncase the array access operator []
*
* The Matrix class provides:
* \n-basic constructors for creating a matrix object from other matrix object,
* \nor by creating empty matrix of a given size,
* \n-input and oput operation via >> and << operators using keyboard or file
* \n-basic operations like access via [] operator, assignment and comparision
*/
template<typename T>
class Matrix : private std::vector<std::vector<T> > {
	typedef std::vector<std::vector<T> > vec;
public:
	using vec::operator[];  // make the array access operator public within Matrix

	// CONSTRUCTORS

	/**
	* Default constructor.  Intialize an empty Matrix object
	* @see Matrix(int Nrows, int Ncols)
	* @see Matrix(const Matrix& m)
	*/
    Matrix();

	/**
	* Alternate constructor.
	* build a matrix Nrows by Ncols
	* @see Matrix()
	* @see Matrix(const Matrix& m)
	* @exception invalid_argument ("matrix size negative or zero")
	*/
    Matrix(int Nrows /**< int. number of rows in matrix */, int Ncols /**< int. number of columns in matrix  */); 

	/**
	* Copy constructor.
	* build a matrix from another matrix
	* @see Matrix()
	* @see Matrix(int Nrows, int Ncols)
	*/
	Matrix(const Matrix<T>& m /**< Matrix&. matrix to copy from  */);

	// ACCESSOR METHODS

	/**
	* Normal public get method.
	* get the number of rows
	* @see int getNcols()const
	* @return int. number of rows in matrix
	*/
	int getNrows() const; // get the number of rows

	/**
	* Normal public get method.
	* get the number of columns
	* @see int getNrows()const
	* @return int. number of columns in matrix
	*/
	int getNcols() const; // get the number of cols

	// OVERLOADED OPERATOR

	/**
	* Overloaded assignment operator
	* @see operator==(const Matrix& m)const
	* @return Matrix&. the matrix on the left of the assignment
	*/
	Matrix<T>& operator=(const Matrix<T>& m /**< Matrix&. Matrix to assign from */); // overloaded assignment operator



	// MULTIPLICATION, COMPARISON METHODS and TRANSPOSE METHODS

	/**
	* Overloaded *operator that returns a Matrix.
	* It Performs matrix by matrix multiplication.
	* @see operator*(const Matrix & a) const
	* @exception out_of_range ("Matrix access error")
	* One or more of the matrix have a zero size
	* @exception std::out_of_range ("uncompatible matrix sizes")
	* Number of columns in first matrix do not match number of columns in second matrix
	* @return Matrix. matrix-matrix product
	*/
	//
	Matrix<T> operator*(const Matrix<T>& a /**< Matrix. matrix to multiply by */
		) const;

	/**
	* Overloaded *operator that returns a Vector.
	* It Performs matrix by vector multiplication.
	* @see operator*(const Matrix & a)const
	* @exception std::out_of_range ("Matrix access error")
	* matrix has a zero size
	* @exception std::out_of_range ("Vector access error")
	* vector has a zero size
	* @exception std::out_of_range ("uncompatible matrix-vector sizes")
	* Number of columns in matrix do not match the vector size
	* @return Vector. matrix-vector product
	*/
	//
	template<typename T1>
	Vector<T1> operator*(const Vector<T1>& v /**< Vector. Vector to multiply by */
		) const
	{
		int nrows = getNrows();
		int ncols = getNcols();

		// catch invalid matrix, vector
		if (nrows <= 0 || ncols <=0 ) { throw std::out_of_range("Matrix access error"); }
		if (v.getSize() <= 0) { throw std::out_of_range("Vector access error"); }

		//if the matrix sizes do not match
		if (ncols != v.getSize()) throw std::out_of_range("matrix sizes do not match");

		// matrix to store the multiplication
		Vector<T1> res(nrows);

		// perform the multiplication
		for (int i = 0;i<nrows;i++)
			for (int j = 0;j<ncols;j++) res[i] = res[i] + ((*this)[i][j] * v[j]);

		// return the result
		return res;
	}



	/**
	* public method that returns the transpose of the matrix.
	* It returns the transpose of matrix
	* @return Matrix.  matrix transpose
	*/
	Matrix transpose() const;


	/**
	* Overloaded istream >> operator.
	* Keyboard input
	* if matrix has size user will be asked to input only matrix values
	* if matrix was not initialized user can choose matrix size and input it values
	* @see operator<<(std::ofstream& ofs, const Matrix& m)
	* @see operator>>(std::istream& is, Matrix& m)
	* @see operator<<(std::ostream& os, const Matrix& m)
	* @exception std::invalid_argument ("read error - negative matrix size");
	* @return std::istream&. The istream object
	*/
	  template <typename T> friend std::istream& operator >> (std::istream& is, /**< Keyboard input stream */
		Matrix<T>& m /**< Matrix to write into */
		);// keyboard input


	/**
	* Overloaded ostream << operator.
	* Display output
	* if matrix has size user will be asked to input only matrix values
	* if matrix was not initialized user can choose matrix size and input it values
	* @see operator>>(std::ifstream& ifs, Matrix& m)
	* @see operator>>(std::istream& is, Matrix& m)
	* @see operator<<(std::ostream& os, const Matrix& m)
	* @return std::ostream&. The ostream object
	*/
	 template <typename T> friend std::ostream& operator<<(std::ostream& os, /**< Display output stream */
		const Matrix<T>& m /**< Matrix to read from*/
		);// screen output


	/**
	* Overloaded ifstream >> operator. File input
	* the file output operator is compatible with file input operator,
	* ie. everything written can be read later.
	* @see operator>>(std::ifstream& ifs, Matrix& m)
	* @see operator<<(std::ofstream& ofs, const Matrix& m)
	* @see operator<<(std::ostream& os, const Matrix& m)
	* @return std::ifstream&. The ifstream object
	*/
	 template <typename T> friend std::ifstream& operator >> (std::ifstream& ifs, /**< Input file stream with opened matrix file */
		Matrix<T>& m /**< Matrix to write into */
		);// file input

	/**
	* Overloaded ofstream << operator. File output
	* the file output operator is compatible with file input operator,
	* ie. everything written can be read later.
	* @see operator>>(std::ifstream& ifs, Matrix& m)
	* @see operator<<(std::ofstream& ofs, const Matrix& m)
	* @see operator>>(std::istream& is, Matrix& m)
	* @exception std::invalid_argument ("file read error - negative matrix size");
	* @return std::ofstream&. The ofstream object
	*/
	 template <typename T> friend std::ofstream& operator<<(std::ofstream& ofs,
		const Matrix<T>& m /**< Matrix to read from*/
		);// file output
};


// CONSTRUCTORS
/*=
 *Default constructor (empty matrix)
 */
template <typename T>
Matrix<T>::Matrix() : std::vector<std::vector<T> >()  {}

/*
 * Alternate constructor - creates a matrix with the given values
 */
template <typename T>
Matrix<T>::Matrix(int Nrows, int Ncols) : std::vector<std::vector<T> >()
{
    //check input
    if(Nrows < 0 || Ncols < 0) throw std::invalid_argument("matrix size negative");
    
	// set the size for the rows
	(*this).resize(Nrows);
	// set the size for the columns
	for (int i = 0; i < Nrows; i++) (*this)[i].resize(Ncols);

	// initialise the matrix to contain zero
	for (int i = 0; i < Nrows; i++)
		for (int j = 0; j < Ncols; j++) (*this)[i][j] = 0.0;
}

/*
 * Copy constructor
 */
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& m) : std::vector<std::vector<double> >()
{
	// set the size of the rows
	(*this).resize(m.size());
	// set the size of the columns
	std::size_t i;
	for (i = 0; i < m.size(); i++) (*this)[i].resize(m[0].size());

	// copy the elements
	for (int i = 0; i < m.getNrows(); i++)
		for (int j = 0; j < m.getNcols(); j++)
			(*this)[i][j] = m[i][j];
}

// ACCESSOR METHODS
/*
* accessor method - get the number of rows
*/
template <typename T>
int Matrix<T> :: getNrows() const
{
	return (*this).size();
}

/*
* accessor method - get the number of columns
*/
template <typename T>
int Matrix<T>::getNcols() const
{
	return (*this)[0].size();
}




// OVERLOADED OPERATORS
/*
* Operator= - assignment
*/
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
	(*this).resize(m.size());
	std::size_t i;
	std::size_t j;
	for (i = 0; i < m.size(); i++) (*this)[i].resize(m[0].size());

	for (i = 0; i<m.size(); i++)
		for (j = 0; j<m[0].size(); j++)
			(*this)[i][j] = m[i][j];
	return *this;
}

/*
* Operator* multiplication of a matrix by a matrix
*/
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& a) const {

	int nrows = getNrows();
	int ncols = getNcols();
	
	// catch invalid matrices
	if (nrows <= 0 || ncols <= 0) { throw std::out_of_range("Matrix access error"); }
	if (a.getNrows() <= 0 || a.getNcols() <= 0) { throw std::out_of_range("Matrix access error"); }

	//if the matrix sizes do not match
	if (ncols != a.getNrows()) throw std::out_of_range("matrix sizes do not match");

	//  matrix to store the result
	Matrix mmult = Matrix(getNrows(), a.getNcols());

	//matrix multiplication
	for (int i = 0;i<getNrows();i++) {
		for (int j = 0;j<a.getNcols();j++) {
			for (int k = 0;k<getNcols();k++) {
				mmult[i][j] += ((*this)[i][k] * a[k][j]);
			}
		}
	}
	return mmult;
}



// OTHER METHODS
/*
 * Transpose of the matrix
 */
template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
	int nrows = getNrows();
	int ncols = getNcols();

	// matrix to stoore the transpose
	Matrix temp(ncols, nrows);

	for (int i=0; i < ncols; i++)
		for (int j=0; j < nrows; j++)
			temp[i][j] = (*this)[j][i];

	return temp;
}


// INPUT AND OUTPUT FUNCTIONS
/*
 * keyboard input , user friendly
 */
template <typename T>
std::istream& operator>>(std::istream& is, Matrix<T>& m) {
	
	int nrows, ncols;

	// test to see whether the matrix m is empty
	if (!m.getNrows()) {
		std::cout << "input the number of rows for the matrix" << std::endl;
		is >> nrows;
		std::cout << "input the number of cols for the matrix" << std::endl;
		is >> ncols;
		//check input 
		if(nrows < 0 || ncols < 0) throw std::invalid_argument("read error - negative matrix size");

		// prepare the matrix to hold n elements
		m = Matrix<T>(nrows, ncols);
	}
    // input the elements
    std::cout << "input "<< m.getNrows() * m.getNcols() << " matrix elements" << std::endl;
	for (int i = 0; i < m.getNrows(); i++)
		for (int j=0; j< m.getNcols(); j++) is >> m[i][j];

    // return the stream object
    return is;
}

/* 
 * screen output, user friendly
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {

	// test to see whether there are any elements
    if (m.getNrows() > 0) {
        os << "The matrix elements are" << std::endl;
        for (int i=0; i<m.getNrows();i++) {
			for (int j=0;j<m.getNcols();j++) {
				os << m[i][j]  << " ";
			}
			os << "\n";
		}
        os << std::endl;
    }
    else
    {
        os << "Matrix is empty." << std::endl;
    }
    return os;
}

/*
* file input - raw data, compatible with file writing operator
*/
template <typename T>
std::ifstream& operator>>(std::ifstream& ifs, Matrix<T>& m) {
	
	int nrows, ncols;

    // read size from the file
    ifs >> nrows; ifs>> ncols;
    //check input sanity
    if(nrows < 0 || ncols < 0) throw std::invalid_argument("file read error - negative matrix size");

    // prepare the vector to hold n elements
    m = Matrix<T>(nrows, ncols);

    // input the elements
    for (int i=0; i<nrows; i++) 
		for (int j=0; j<ncols; j++) ifs >> m[i][j];

    // return the stream object
    return ifs;
}


/*
* file output - raw data, comaptible with file reading operator
*/
template <typename T>
std::ofstream& operator<<(std::ofstream& ofs, const Matrix<T>& m) {
    //put matrix rownumber in first line (even if it is zero)
    ofs << m.getNrows() << std::endl;
	//put matrix columnnumber in second line (even if it is zero)
	ofs << m.getNcols() << std::endl;
    //put data in third line (if size==zero nothing will be put)
    for (int i=0; i<m.getNrows(); i++) {
		for (int j=0; j<m.getNcols(); j++) ofs << m[i][j] <<  " ";
		ofs << std::endl;
	}
    return ofs;
}


#endif