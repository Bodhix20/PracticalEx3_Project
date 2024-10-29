#ifndef _CDATA
#define _CDATA

#include "textobject.h"
#include "vector.h"
#include "point.h"


// class encapsulating the data points and parameter values for the curve data fitting
class Cdata : public TextObject
{
private:
	int n;					// number of points
	Vector<Point3D> pts;	// data points
	Vector<double> param;	// parameter values
	Vector<Point3D> tang;	// tangent vectors
	void uniform();			// uniform parameetrisation method
	void chord();			// chord length parameetrisation method
	void centripetal();		// centripetal parameetrisation method
	void computeTangents();	// computes the tangents according to Bessel's method
public:
	// constructors
	Cdata();
	Cdata(int num, const Vector<Point3D>& data);

	// accessor methods
	int getNum() const;
	Vector<Point3D> getPoints() const;
	Vector<double> getParam() const;
	Vector<Point3D> getTang() const;

	// read/write methods
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) const;
	virtual void write(std::ostream& os) const;
	virtual void read(std::istream &is);
};

#endif