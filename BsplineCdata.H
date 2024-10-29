#ifndef _BSPLINECDATA_CURVE
#define _BSPLINECDATA_CURVE

#include "bsplinecurve.h"
#include "cdata.h"
#include "matrix.h"


// class for encapsulating an interpolating or approximating Bspline curve
class BsplineCdata {
private:
	Cdata cdat;			// the Cdata object containing the data points, param,eter values etc
	BsplineCurve b;		// the underlying Bspline curve

	// helper methods to compute the interpolating and least squares approximating curves 
	Matrix<double> computeLeastSquaresMatrix(int nseg, int ord) const;
	Vector<Point3D> computeRHSVectorLeastSq(int nseg, int ord) const;
	Matrix<double> computeInterpolationMatrix() const;
	Vector<Point3D> computeRHSVectorInterp() const;
	Vector<double> createVectorInterp(int num, int ord, const Vector<double>& knots, double val) const;
	int findIndex(int ord, int dim, const Vector<double>& knots, double t0) const;
	Vector<double> computeKnots(const Vector<double>& param, int n, int ord) const;
	void convert(const AbsCurve<Point3D>& c, int num);
public:

	// constructors
	BsplineCdata();
	BsplineCdata(const Cdata& data);

	// accessor method
	BsplineCurve getBspline() const;

	// compute the interpolating and approximating Bsplines
	void computeBsplineInterp();
	void computeBsplineLeastSq(int nseg, int ord);
};

#endif