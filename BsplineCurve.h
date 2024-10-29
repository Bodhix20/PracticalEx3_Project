#ifndef _BSPLINE_CURVE
#define _BSPLINE_CURVE

#include "abscurve.h"
#include "point.h"


// class encapulating a polynomial curve written in Bspline form (i.e using the Bspline basis) 
class BsplineCurve : public AbsCurve<Point3D>
{
private:
	int ord;						// order ot the curve
	int dim;						// dimension (= number of control points)
	Vector<Point3D> d;				// control points
	Vector<double> knots;			// vector of knots
	int findIndex(double t0) const; // helper function to find the senment where t0 lies
public:

	// constructors
	BsplineCurve();
	BsplineCurve(double t1, double t2, int k, int num, const Vector<Point3D>& cpts, const Vector<double>& kts);

	// overloaded operator for evaluating the Bspline curve at t0
	Point3D operator()(double t0) const;

	// accessor methods
	Vector<Point3D> getCpoints() const;
	Vector<double> getKnots() const;

	// read/write methods
	void readfile(std::ifstream& ifs);
	void writefile(std::ofstream& ofs) const;
	void write(std::ostream& os) const;
	void read(std::istream &is);
};

#endif
