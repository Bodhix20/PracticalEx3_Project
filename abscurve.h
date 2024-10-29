#ifndef _ABS_CURVE
#define _ABS_CURVE

#include "textobject.h"
#include "vector.h"


// abstract base curve template class encapsulating a generic curve defined by two limits
// pure virtual methods are the read/write functions and the overloaded operator for evaluation of the curve at a point
template<class T>
class AbsCurve : public TextObject
{
protected:
	double x1; // left limit
	double x2; // right limit
public:
	// constructors
	AbsCurve();
	AbsCurve(double X1, double X2);

	// pure virtual method for evaluation
	virtual T operator()(double) const = 0;

	// compute an array of points along the curve 
	virtual Vector<T> computePoints(int n) const;

	// accessor methods getting the limits
	virtual double getLeftLimit() const;
	virtual double getRightLimit() const;

	// purer virtual read/write methods
	virtual void readfile(std::ifstream& ifs) = 0;
	virtual void writefile(std::ofstream& ofs) const = 0;
	virtual void write(std::ostream& os) const = 0;
	virtual void read(std::istream &is) = 0;
};


// constructor implemntations
template<class T>
AbsCurve<T>::AbsCurve() : x1(0.0), x2(0.0) {}


template<class T>
AbsCurve<T>::AbsCurve(double X1, double X2): x1(X1), x2(X2) {} 


// compute an array of points using the overloaded operator
template<class T>
Vector<T> AbsCurve<T>::computePoints(int n) const
{
	Vector<T> v(n);

	double step = (x2-x1)/(double)(n-1);
	
	// evaluate
	double val;
	
	for (int i=0; i<n; i++) {
		val = x1+i*step;
		v[i] = (*this)(val);
	}
	return v;
}


// accessor methods
template<class T>
double AbsCurve<T>::getLeftLimit() const
{
	return x1;
}

template<class T>
double AbsCurve<T>::getRightLimit() const
{
	return x2;
}

#endif